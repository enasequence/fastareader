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

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Random;
import java.util.zip.GZIPOutputStream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class BgzfByteReaderTest {

    @TempDir
    Path tmp;

    @Test
    void detectsBgzfPlainGzipAndUncompressed() throws IOException {
        byte[] payload = "ACGTACGTACGTNNNN".getBytes(StandardCharsets.US_ASCII);

        File bgzf = BgzfTestWriter.writeBgzfFile(tmp.resolve("a.bgzf"), payload, 8);
        assertEquals(BgzfDetector.Compression.BGZF, BgzfDetector.detect(bgzf));

        Path plain = tmp.resolve("a.txt");
        Files.write(plain, payload);
        assertEquals(BgzfDetector.Compression.UNCOMPRESSED, BgzfDetector.detect(plain.toFile()));

        Path gz = tmp.resolve("a.gz");
        try (GZIPOutputStream gos = new GZIPOutputStream(Files.newOutputStream(gz))) {
            gos.write(payload);
        }
        assertEquals(BgzfDetector.Compression.PLAIN_GZIP, BgzfDetector.detect(gz.toFile()));
    }

    @Test
    void factoryRejectsPlainGzip() throws IOException {
        Path gz = tmp.resolve("b.gz");
        try (GZIPOutputStream gos = new GZIPOutputStream(Files.newOutputStream(gz))) {
            gos.write("hello".getBytes(StandardCharsets.US_ASCII));
        }
        IOException e = assertThrows(IOException.class, () -> SeekableByteReaderFactory.open(gz.toFile()));
        assertTrue(e.getMessage().toLowerCase().contains("bgzip"));
    }

    @Test
    void readsBackExactBytesAcrossBlocks() throws IOException {
        Random rnd = new Random(42);
        byte[] payload = new byte[5000];
        for (int i = 0; i < payload.length; i++) {
            payload[i] = (byte) ("ACGTN".charAt(rnd.nextInt(5)));
        }
        // force tiny blocks so the data spans many blocks
        File bgzf = BgzfTestWriter.writeBgzfFile(tmp.resolve("c.bgzf"), payload, 100);

        try (BgzfByteReader reader = new BgzfByteReader(bgzf)) {
            assertEquals(payload.length, reader.size());

            // read the whole stream in a 64-byte buffer that straddles block boundaries
            byte[] got = new byte[payload.length];
            int total = 0;
            ByteBuffer buf = ByteBuffer.allocate(64);
            while (total < payload.length) {
                buf.clear();
                int n = reader.read(buf, total);
                if (n <= 0) break;
                buf.flip();
                buf.get(got, total, n);
                total += n;
            }
            assertEquals(payload.length, total);
            assertArrayEquals(payload, got);

            // EOF behaviour
            ByteBuffer one = ByteBuffer.allocate(1);
            assertEquals(-1, reader.read(one, payload.length));
        }
    }

    @Test
    void positionalReadDoesNotAffectCursor() throws IOException {
        byte[] payload = "0123456789ABCDEF".getBytes(StandardCharsets.US_ASCII);
        File bgzf = BgzfTestWriter.writeBgzfFile(tmp.resolve("d.bgzf"), payload, 4);

        try (BgzfByteReader reader = new BgzfByteReader(bgzf)) {
            reader.position(5);
            ByteBuffer buf = ByteBuffer.allocate(3);
            int n = reader.read(buf, 10);
            buf.flip();
            assertEquals(3, n);
            assertEquals("ABC", StandardCharsets.US_ASCII.decode(buf).toString());
            assertEquals(5, reader.position(), "positional read must not move the cursor");
        }
    }
}
