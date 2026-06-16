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
package uk.ac.ebi.embl.fastareader.encoding;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;

public final class Utf8Detector {

    public static final int DEFAULT_MAX_BYTES = 4 * 1024 * 1024;

    private Utf8Detector() {}

    public static boolean isProbablyUtf8(Path path) throws IOException {
        return isProbablyUtf8(path, DEFAULT_MAX_BYTES);
    }

    public static boolean isProbablyUtf8(Path path, int maxBytes) throws IOException {
        if (maxBytes <= 0) throw new IllegalArgumentException("maxBytes must be > 0");

        try (InputStream in = Files.newInputStream(path)) {
            return isProbablyUtf8(in, maxBytes);
        }
    }

    /**
     * Validates that up to {@code maxBytes} bytes drawn from {@code in} look like UTF-8.
     *
     * <p>This overload validates an arbitrary byte stream rather than a file path, which lets
     * callers feed a <em>decompressed</em> prefix (e.g. of a BGZF file) into the same state
     * machine. The caller owns the stream and is responsible for closing it.</p>
     */
    public static boolean isProbablyUtf8(InputStream in, int maxBytes) throws IOException {
        if (maxBytes <= 0) throw new IllegalArgumentException("maxBytes must be > 0");

        {
            byte[] buf = new byte[1024 * 1024];

            int totalRead = 0;
            boolean bomChecked = false;

            int needed = 0; // continuation bytes remaining
            int minCode = 0; // smallest allowed code point for this sequence (prevents overlong)
            int code = 0; // current code point being assembled

            while (totalRead < maxBytes) {
                int toRead = Math.min(buf.length, maxBytes - totalRead);
                int n = in.read(buf, 0, toRead);
                if (n < 0) break;
                if (n == 0) continue;

                int i = 0;

                if (!bomChecked) {
                    // Check UTF-8 BOM if present at the very start of the file.
                    bomChecked = true;
                    if (n >= 3 && (buf[0] & 0xFF) == 0xEF && (buf[1] & 0xFF) == 0xBB && (buf[2] & 0xFF) == 0xBF) {
                        i = 3;
                    }
                }

                totalRead += n;

                for (; i < n; i++) {
                    int b = buf[i] & 0xFF;

                    if (needed == 0) {
                        // ASCII fast path
                        if (b < 0x80) continue;

                        // Determine sequence length & initial code bits
                        if (b >= 0xC2 && b <= 0xDF) { // 2-byte
                            needed = 1;
                            code = b & 0x1F;
                            minCode = 0x80;
                        } else if (b >= 0xE0 && b <= 0xEF) { // 3-byte
                            needed = 2;
                            code = b & 0x0F;
                            minCode = 0x800;
                        } else if (b >= 0xF0 && b <= 0xF4) { // 4-byte (UTF-8 valid max is U+10FFFF)
                            needed = 3;
                            code = b & 0x07;
                            minCode = 0x10000;
                        } else {
                            // Includes 0x80..0xBF continuation, 0xC0..0xC1 overlong starts, 0xF5..0xFF out of range
                            return false;
                        }
                    } else {
                        // Must be continuation byte 10xxxxxx
                        if ((b & 0xC0) != 0x80) return false;

                        code = (code << 6) | (b & 0x3F);
                        needed--;

                        if (needed == 0) {
                            if (code < minCode) return false; // Reject overlong encodings
                            if (code >= 0xD800 && code <= 0xDFFF) return false; // Reject UTF-16 surrogate halves
                            if (code > 0x10FFFF) return false; // Reject > U+10FFFF

                            // reset for next char
                            code = 0;
                            minCode = 0;
                        }
                    }
                }
            }

            // If we end sampling in the middle of a character -> invalid
            return needed == 0;
        }
    }
}
