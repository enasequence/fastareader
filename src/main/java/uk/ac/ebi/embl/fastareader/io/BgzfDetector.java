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

/**
 * Detects whether a file is BGZF-compressed, plain gzip, or uncompressed by inspecting its
 * leading bytes (gzip magic + the BGZF {@code BC} extra subfield).
 *
 * <p>BGZF files are gzip files whose header carries a BGZF extra subfield
 * ({@code SI1='B' (66)}, {@code SI2='C' (67)}, {@code SLEN=2}). Plain gzip files have the
 * gzip magic but lack that subfield; they are detected separately so callers can fail fast
 * with a helpful message (plain gzip is not seekable).</p>
 */
public final class BgzfDetector {

    /** Result of inspecting a file's leading bytes. */
    public enum Compression {
        /** Not a gzip file at all — read it as a plain (uncompressed) file. */
        UNCOMPRESSED,
        /** A BGZF-compressed file (gzip with the BGZF {@code BC} extra subfield). */
        BGZF,
        /** A plain gzip file (gzip magic but no BGZF {@code BC} subfield) — not seekable. */
        PLAIN_GZIP
    }

    private static final int GZIP_ID1 = 0x1f;
    private static final int GZIP_ID2 = 0x8b;
    private static final int DEFLATE_CM = 0x08;
    private static final int FEXTRA = 0x04;
    private static final int BGZF_SI1 = 66; // 'B'
    private static final int BGZF_SI2 = 67; // 'C'
    private static final int BGZF_SLEN = 2;

    private BgzfDetector() {}

    /** Detects the compression of {@code file} by reading only its leading header bytes. */
    public static Compression detect(File file) throws IOException {
        try (FileChannel channel = FileChannel.open(file.toPath(), StandardOpenOption.READ)) {
            byte[] header = readUpTo(channel, 0, 12);
            if (header.length < 2) return Compression.UNCOMPRESSED;
            if ((header[0] & 0xFF) != GZIP_ID1 || (header[1] & 0xFF) != GZIP_ID2) {
                return Compression.UNCOMPRESSED;
            }
            // From here on it is gzip; decide BGZF vs plain gzip.
            if (header.length < 12) return Compression.PLAIN_GZIP;
            if ((header[2] & 0xFF) != DEFLATE_CM) return Compression.PLAIN_GZIP;
            if ((header[3] & 0xFF & FEXTRA) == 0) return Compression.PLAIN_GZIP;

            int xlen = (header[10] & 0xFF) | ((header[11] & 0xFF) << 8);
            byte[] extra = readUpTo(channel, 12, xlen);
            if (extra.length < xlen) return Compression.PLAIN_GZIP;

            return hasBgzfSubfield(extra, xlen) ? Compression.BGZF : Compression.PLAIN_GZIP;
        }
    }

    /** Scans the gzip extra field for the BGZF {@code BC} subfield. */
    private static boolean hasBgzfSubfield(byte[] extra, int xlen) {
        int pos = 0;
        while (pos + 4 <= xlen) {
            int si1 = extra[pos] & 0xFF;
            int si2 = extra[pos + 1] & 0xFF;
            int slen = (extra[pos + 2] & 0xFF) | ((extra[pos + 3] & 0xFF) << 8);
            if (si1 == BGZF_SI1 && si2 == BGZF_SI2 && slen == BGZF_SLEN) return true;
            pos += 4 + slen;
        }
        return false;
    }

    /** Reads up to {@code len} bytes at {@code position}, returning however many were available. */
    private static byte[] readUpTo(FileChannel channel, long position, int len) throws IOException {
        if (len <= 0) return new byte[0];
        ByteBuffer buf = ByteBuffer.allocate(len);
        long pos = position;
        while (buf.hasRemaining()) {
            int n = channel.read(buf, pos);
            if (n < 0) break;
            pos += n;
        }
        if (buf.position() == len) return buf.array();
        byte[] out = new byte[buf.position()];
        System.arraycopy(buf.array(), 0, out, 0, out.length);
        return out;
    }
}
