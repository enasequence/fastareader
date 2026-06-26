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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

/**
 * Test-only helper that writes BGZF files using only {@code java.util.zip}, so fixtures can be
 * produced hermetically without requiring {@code bgzip}/{@code samtools} on the build machine.
 *
 * <p>Each call to {@link #toBgzf(byte[], int)} wraps the payload into one or more standard,
 * independent gzip blocks (each carrying the BGZF {@code BC} extra subfield, {@code CRC32} and
 * {@code ISIZE}) terminated by the canonical empty EOF block. The block size is configurable so
 * tests can force tiny blocks and exercise cross-block reads.</p>
 *
 * <p>This class lives under {@code src/test} only and must never leak into the published
 * artifact.</p>
 */
public final class BgzfTestWriter {

    /** Maximum uncompressed payload per BGZF block as mandated by the format. */
    public static final int MAX_BLOCK_SIZE = 65_536;

    /** The canonical 28-byte empty BGZF EOF block (ISIZE 0). */
    private static final byte[] EOF_BLOCK = {
        0x1f,
        (byte) 0x8b,
        0x08,
        0x04,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        (byte) 0xff,
        0x06,
        0x00,
        0x42,
        0x43,
        0x02,
        0x00,
        0x1b,
        0x00,
        0x03,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00
    };

    private BgzfTestWriter() {}

    /**
     * Returns valid BGZF bytes (single block) with the CRC32 of that block corrupted.
     * Triggers the CRC32 check in {@code BgzfByteReader.inflateBlock}.
     */
    public static byte[] toBgzfCorruptCrc(byte[] data) throws IOException {
        byte[] bytes = toBgzf(data, Math.max(data.length, 1));
        int blockEnd = bytes.length - EOF_BLOCK.length;
        bytes[blockEnd - 8] ^= 0xFF;
        return bytes;
    }

    /**
     * Returns valid BGZF bytes (single block) with the ISIZE trailer set to 65537,
     * triggering the oversized-ISIZE guard in {@code BgzfByteReader.buildDirectory}.
     */
    public static byte[] toBgzfOversizedIsize(byte[] data) throws IOException {
        byte[] bytes = toBgzf(data, Math.max(data.length, 1));
        int blockEnd = bytes.length - EOF_BLOCK.length;
        // 65537 = 0x00010001 in little-endian
        bytes[blockEnd - 4] = 0x01;
        bytes[blockEnd - 3] = 0x00;
        bytes[blockEnd - 2] = 0x01;
        bytes[blockEnd - 1] = 0x00;
        return bytes;
    }

    /**
     * Returns a hand-crafted byte array that passes {@link BgzfDetector} (BC subfield header
     * present) but fails {@code BgzfByteReader.findBsize} because XLEN is too small to hold
     * both bytes of the BSIZE value (bug-1 reproducer).
     */
    public static byte[] toBgzfShortBcXlen() {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        // gzip header: ID1 ID2 CM FLG MTIME(4) XFL OS
        out.write(0x1f);
        out.write(0x8b);
        out.write(0x08);
        out.write(0x04);
        out.write(0);
        out.write(0);
        out.write(0);
        out.write(0);
        out.write(0);
        out.write(0xff);
        // XLEN = 5 (only one byte of BSIZE fits)
        out.write(5);
        out.write(0);
        // extra: SI1='B' SI2='C' SLEN=2 + 1 byte of BSIZE (truncated)
        out.write(0x42);
        out.write(0x43);
        out.write(0x02);
        out.write(0x00);
        out.write(0xFF);
        // append EOF block so file has more than the header
        for (byte b : EOF_BLOCK) out.write(b & 0xFF);
        return out.toByteArray();
    }

    /**
     * Returns BGZF bytes truncated to 10 bytes, always landing mid-gzip-header (headers need ≥18
     * bytes), so {@code BgzfByteReader} always hits a premature-EOF error regardless of payload.
     */
    public static byte[] toBgzfTruncated(byte[] data) throws IOException {
        byte[] bytes = toBgzf(data, Math.max(data.length, 1));
        return Arrays.copyOf(bytes, 10);
    }

    /**
     * Returns BGZF bytes with the CRC32 of the {@code blockIndex}-th data block (0-based) flipped,
     * leaving all other blocks intact so only reads touching that block throw.
     */
    public static byte[] toBgzfCorruptCrcOfBlock(byte[] data, int blockSize, int blockIndex) throws IOException {
        byte[] bytes = toBgzf(data, blockSize);
        int offset = 0;
        for (int i = 0; i <= blockIndex; i++) {
            int xlen = (bytes[offset + 10] & 0xFF) | ((bytes[offset + 11] & 0xFF) << 8);
            int bsize = -1;
            for (int pos = 0; pos + 4 <= xlen; ) {
                int si1 = bytes[offset + 12 + pos] & 0xFF;
                int si2 = bytes[offset + 12 + pos + 1] & 0xFF;
                int slen = (bytes[offset + 12 + pos + 2] & 0xFF) | ((bytes[offset + 12 + pos + 3] & 0xFF) << 8);
                if (si1 == 0x42 && si2 == 0x43 && slen == 2) {
                    bsize = (bytes[offset + 12 + pos + 4] & 0xFF) | ((bytes[offset + 12 + pos + 5] & 0xFF) << 8);
                    break;
                }
                pos += 4 + slen;
            }
            if (bsize < 0) throw new IllegalStateException("no BC subfield at block " + i);
            int thisBlockSize = bsize + 1;
            if (i == blockIndex) {
                bytes[offset + thisBlockSize - 8] ^= 0xFF;
                return bytes;
            }
            offset += thisBlockSize;
        }
        throw new IllegalArgumentException("not enough blocks");
    }

    /** Compresses {@code data} into BGZF bytes using the given uncompressed block size. */
    public static byte[] toBgzf(byte[] data, int blockSize) throws IOException {
        if (blockSize <= 0 || blockSize > MAX_BLOCK_SIZE) {
            throw new IllegalArgumentException("blockSize must be in 1.." + MAX_BLOCK_SIZE);
        }
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        int offset = 0;
        while (offset < data.length) {
            int len = Math.min(blockSize, data.length - offset);
            out.write(makeBlock(data, offset, len));
            offset += len;
        }
        out.write(EOF_BLOCK);
        return out.toByteArray();
    }

    /** Writes {@code data} as a BGZF file at {@code target}, returning it as a {@link File}. */
    public static File writeBgzfFile(Path target, byte[] data, int blockSize) throws IOException {
        Files.write(target, toBgzf(data, blockSize));
        return target.toFile();
    }

    /** Builds a single BGZF block for {@code data[offset..offset+len)}. */
    private static byte[] makeBlock(byte[] data, int offset, int len) throws IOException {
        Deflater deflater = new Deflater(Deflater.DEFAULT_COMPRESSION, true);
        deflater.setInput(data, offset, len);
        deflater.finish();
        byte[] compressed = new byte[len + 64];
        int compressedLen = 0;
        while (!deflater.finished()) {
            if (compressedLen == compressed.length) {
                byte[] bigger = new byte[compressed.length * 2];
                System.arraycopy(compressed, 0, bigger, 0, compressedLen);
                compressed = bigger;
            }
            compressedLen += deflater.deflate(compressed, compressedLen, compressed.length - compressedLen);
        }
        deflater.end();

        CRC32 crc = new CRC32();
        crc.update(data, offset, len);

        // total block size = 18 (header incl. 6-byte extra field) + compressed + 8 (trailer)
        int blockSize = 18 + compressedLen + 8;
        int bsize = blockSize - 1;

        ByteArrayOutputStream block = new ByteArrayOutputStream(blockSize);
        // gzip header (10 bytes)
        block.write(0x1f); // ID1
        block.write(0x8b); // ID2
        block.write(0x08); // CM = deflate
        block.write(0x04); // FLG = FEXTRA
        block.write(0x00); // MTIME (4 bytes)
        block.write(0x00);
        block.write(0x00);
        block.write(0x00);
        block.write(0x00); // XFL
        block.write(0xff); // OS = unknown
        // XLEN = 6 (LE)
        block.write(0x06);
        block.write(0x00);
        // BGZF extra subfield: SI1='B', SI2='C', SLEN=2, BSIZE (LE)
        block.write(0x42); // 'B'
        block.write(0x43); // 'C'
        block.write(0x02);
        block.write(0x00);
        block.write(bsize & 0xff);
        block.write((bsize >> 8) & 0xff);
        // compressed payload
        block.write(compressed, 0, compressedLen);
        // trailer: CRC32 (LE) + ISIZE (LE)
        long crcValue = crc.getValue();
        block.write((int) (crcValue & 0xff));
        block.write((int) ((crcValue >> 8) & 0xff));
        block.write((int) ((crcValue >> 16) & 0xff));
        block.write((int) ((crcValue >> 24) & 0xff));
        block.write(len & 0xff);
        block.write((len >> 8) & 0xff);
        block.write((len >> 16) & 0xff);
        block.write((len >> 24) & 0xff);

        return block.toByteArray();
    }
}
