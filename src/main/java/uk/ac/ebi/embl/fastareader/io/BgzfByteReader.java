/*
 * Copyright 2026 EMBL - European Bioinformatics Institute
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
 * file except in compliance with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software distributed under the
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 */
package uk.ac.ebi.embl.fastareader.io;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * {@link SeekableByteReader} over a BGZF-compressed file, presenting the decompressed
 * (logical) byte stream.
 *
 * <p>At open time it builds a <em>block directory</em> — an ordered list of
 * {@code (compressedOffset, uncompressedStart, uncompressedSize)} per BGZF block — by reading
 * only each block's gzip header (for {@code BSIZE}) and trailer (for {@code ISIZE}), without
 * inflating any payload. Logical reads binary-search the directory, inflate the covering
 * block(s) (cached via a small LRU), and copy out the requested bytes.</p>
 *
 * <p>Logical offsets are flat and additive and equal the offsets of the file's decompressed
 * twin, so the {@code SequenceIndex} built from a BGZF file is identical to its uncompressed
 * equivalent.</p>
 *
 * <p>This class is not thread-safe. Concurrent reads require external synchronization or
 * separate reader instances.</p>
 */
public final class BgzfByteReader implements SeekableByteReader {

    /** Number of decompressed blocks held in the LRU cache (≤ ~1 MiB at 16 × 64 KiB). */
    private static final int CACHE_CAPACITY = 16;

    private static final int GZIP_ID1 = 0x1f;
    private static final int GZIP_ID2 = 0x8b;
    private static final int DEFLATE_CM = 0x08;
    private static final int FEXTRA = 0x04;
    private static final int BGZF_SI1 = 66; // 'B'
    private static final int BGZF_SI2 = 67; // 'C'
    private static final int GZIP_TRAILER_SIZE = 8; // CRC32 (4) + ISIZE (4)

    private final FileChannel channel;

    // Block directory (parallel arrays sorted by uncompressedStart).
    private final long[] compressedOffsets;
    private final long[] uncompressedStarts;
    private final int[] uncompressedSizes;
    private final int[] blockSizes;
    private final int[] dataOffsets; // offset within block where the deflate payload begins
    private final int[] crc32s; // gzip trailer CRC32 for each block
    private final long totalUncompressed;

    // LRU cache of decompressed block payloads keyed by block index.
    private final Map<Integer, byte[]> cache;

    private long position = 0;

    public BgzfByteReader(File file) throws IOException {
        this.channel = FileChannel.open(file.toPath(), StandardOpenOption.READ);
        boolean close = true;
        try {
            List<long[]> directory = buildDirectory();
            int n = directory.size();
            this.compressedOffsets = new long[n];
            this.uncompressedStarts = new long[n];
            this.uncompressedSizes = new int[n];
            this.blockSizes = new int[n];
            this.dataOffsets = new int[n];
            this.crc32s = new int[n];
            long total = 0;
            for (int i = 0; i < n; i++) {
                long[] b = directory.get(i);
                compressedOffsets[i] = b[0];
                uncompressedStarts[i] = b[1];
                uncompressedSizes[i] = (int) b[2];
                blockSizes[i] = (int) b[3];
                dataOffsets[i] = (int) b[4];
                crc32s[i] = (int) b[5];
                total += b[2];
            }
            this.totalUncompressed = total;
            this.cache = new LinkedHashMap<>(CACHE_CAPACITY * 2, 0.75f, true) {
                @Override
                protected boolean removeEldestEntry(Map.Entry<Integer, byte[]> eldest) {
                    return size() > CACHE_CAPACITY;
                }
            };
            close = false;
        } finally {
            if (close) channel.close();
        }
    }

    /** Walks the file block-by-block, reading only headers + trailers (no inflation). */
    private List<long[]> buildDirectory() throws IOException {
        long fileSize = channel.size();
        List<long[]> directory = new ArrayList<>();
        long compressedOffset = 0;
        long uncompressedStart = 0;

        while (compressedOffset < fileSize) {
            byte[] head = readAt(compressedOffset, 12);
            if ((head[0] & 0xFF) != GZIP_ID1 || (head[1] & 0xFF) != GZIP_ID2 || (head[2] & 0xFF) != DEFLATE_CM) {
                throw new IOException("Malformed BGZF: bad gzip header at compressed offset " + compressedOffset);
            }
            if ((head[3] & 0xFF & FEXTRA) == 0) {
                throw new IOException("Malformed BGZF: missing FEXTRA flag at compressed offset " + compressedOffset);
            }
            int xlen = (head[10] & 0xFF) | ((head[11] & 0xFF) << 8);
            byte[] extra = readAt(compressedOffset + 12, xlen);
            int bsize = findBsize(extra, xlen, compressedOffset);
            int blockSize = bsize + 1;
            int dataOffset = 12 + xlen;
            if (blockSize <= dataOffset + GZIP_TRAILER_SIZE) {
                throw new IOException("Malformed BGZF: implausible BSIZE at compressed offset " + compressedOffset);
            }

            byte[] tail = readAt(compressedOffset + blockSize - GZIP_TRAILER_SIZE, GZIP_TRAILER_SIZE);
            int crc32 =
                    (tail[0] & 0xFF) | ((tail[1] & 0xFF) << 8) | ((tail[2] & 0xFF) << 16) | ((tail[3] & 0xFF) << 24);
            long isize = (tail[4] & 0xFFL)
                    | ((tail[5] & 0xFFL) << 8)
                    | ((tail[6] & 0xFFL) << 16)
                    | ((tail[7] & 0xFFL) << 24);
            if (isize > 65536) {
                throw new IOException("Malformed BGZF: ISIZE " + isize + " exceeds 65536 at compressed offset " + compressedOffset);
            }

            if (isize > 0) {
                directory.add(new long[] {compressedOffset, uncompressedStart, isize, blockSize, dataOffset, crc32});
                uncompressedStart += isize;
            }
            compressedOffset += blockSize;
        }
        return directory;
    }

    /** Reads the BGZF {@code BC} subfield's {@code BSIZE} value from the gzip extra field. */
    private int findBsize(byte[] extra, int xlen, long compressedOffset) throws IOException {
        int pos = 0;
        while (pos + 4 <= xlen) {
            int si1 = extra[pos] & 0xFF;
            int si2 = extra[pos + 1] & 0xFF;
            int slen = (extra[pos + 2] & 0xFF) | ((extra[pos + 3] & 0xFF) << 8);
            if (si1 == BGZF_SI1 && si2 == BGZF_SI2 && slen == 2) {
                if (pos + 6 > xlen) {
                    throw new IOException("Malformed BGZF: BC subfield truncated at compressed offset " + compressedOffset);
                }
                return (extra[pos + 4] & 0xFF) | ((extra[pos + 5] & 0xFF) << 8);
            }
            pos += 4 + slen;
        }
        throw new IOException("Malformed BGZF: missing BC subfield at compressed offset " + compressedOffset);
    }

    @Override
    public int read(ByteBuffer dst, long position) throws IOException {
        if (position < 0) throw new IllegalArgumentException("Negative position: " + position);
        if (!channel.isOpen()) throw new IOException("Reader is closed");
        if (position >= totalUncompressed) return -1;

        int total = 0;
        long pos = position;
        while (dst.hasRemaining() && pos < totalUncompressed) {
            int bi = blockIndexFor(pos);
            byte[] block = blockBytes(bi);
            int within = (int) (pos - uncompressedStarts[bi]);
            int n = Math.min(block.length - within, dst.remaining());
            dst.put(block, within, n);
            pos += n;
            total += n;
        }
        return total;
    }

    @Override
    public long size() {
        return totalUncompressed;
    }

    @Override
    public long position() {
        return position;
    }

    @Override
    public SeekableByteReader position(long newPosition) {
        if (newPosition < 0) throw new IllegalArgumentException("Negative position: " + newPosition);
        this.position = newPosition;
        return this;
    }

    @Override
    public boolean isOpen() {
        return channel.isOpen();
    }

    @Override
    public void close() throws IOException {
        cache.clear();
        channel.close();
    }

    // ------------------------------------------------------------------

    /** Largest block index whose {@code uncompressedStart <= pos}. */
    private int blockIndexFor(long pos) {
        int lo = 0;
        int hi = uncompressedStarts.length - 1;
        int ans = 0;
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            if (uncompressedStarts[mid] <= pos) {
                ans = mid;
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return ans;
    }

    private byte[] blockBytes(int index) throws IOException {
        byte[] cached = cache.get(index);
        if (cached != null) return cached;
        byte[] inflated = inflateBlock(index);
        cache.put(index, inflated);
        return inflated;
    }

    private byte[] inflateBlock(int index) throws IOException {
        long offset = compressedOffsets[index];
        int blockSize = blockSizes[index];
        int dataOffset = dataOffsets[index];
        int dataLen = blockSize - dataOffset - GZIP_TRAILER_SIZE;

        byte[] compressed = readAt(offset + dataOffset, dataLen);
        byte[] out = new byte[uncompressedSizes[index]];

        Inflater inflater = new Inflater(true);
        try {
            inflater.setInput(compressed);
            int total = 0;
            while (total < out.length) {
                int n = inflater.inflate(out, total, out.length - total);
                if (n == 0) {
                    if (inflater.finished() || inflater.needsInput() || inflater.needsDictionary()) break;
                }
                total += n;
            }
            if (total != out.length) {
                throw new IOException("Malformed BGZF: inflated " + total + " bytes, expected ISIZE " + out.length
                        + " at compressed offset " + offset);
            }
        } catch (DataFormatException e) {
            throw new IOException("Malformed BGZF: deflate error at compressed offset " + offset, e);
        } finally {
            inflater.end();
        }

        CRC32 crc = new CRC32();
        crc.update(out, 0, out.length);
        if ((int) crc.getValue() != crc32s[index]) {
            throw new IOException("Malformed BGZF: CRC32 mismatch at compressed offset " + offset);
        }

        return out;
    }

    /** Reads exactly {@code len} bytes at {@code position}; throws on premature EOF. */
    private byte[] readAt(long position, int len) throws IOException {
        if (len <= 0) return new byte[0];
        ByteBuffer buf = ByteBuffer.allocate(len);
        long pos = position;
        while (buf.hasRemaining()) {
            int n = channel.read(buf, pos);
            if (n < 0) {
                throw new IOException("Malformed BGZF: unexpected EOF reading " + len + " bytes at offset " + position);
            }
            pos += n;
        }
        return buf.array();
    }
}
