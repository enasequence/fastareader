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
package uk.ac.ebi.embl.fastareader;

import java.io.IOException;
import java.nio.ByteBuffer;
import uk.ac.ebi.embl.fastareader.io.SeekableByteReader;

class BgzfReaderUtils {

    private BgzfReaderUtils() {}

    /** Reads up to {@code maxBytes} of the logical (decompressed) content from {@code reader}. */
    static byte[] decompressedPrefix(SeekableByteReader reader, int maxBytes) throws IOException {
        int prefix = (int) Math.min(reader.size(), (long) maxBytes);
        ByteBuffer buf = ByteBuffer.allocate(prefix);
        long pos = 0;
        while (buf.hasRemaining()) {
            int n = reader.read(buf, pos);
            if (n <= 0) break;
            pos += n;
        }
        byte[] out = new byte[buf.position()];
        System.arraycopy(buf.array(), 0, out, 0, out.length);
        return out;
    }
}
