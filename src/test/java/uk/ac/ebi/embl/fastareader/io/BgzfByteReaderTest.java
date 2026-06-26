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

    // ---- BgzfByteReader contract tests ----

    @Test
    void emptyBgzfReturnsMinusOne() throws IOException {
        File bgzf = BgzfTestWriter.writeBgzfFile(tmp.resolve("empty.bgzf"), new byte[0], BgzfTestWriter.MAX_BLOCK_SIZE);
        try (BgzfByteReader reader = new BgzfByteReader(bgzf)) {
            assertEquals(0, reader.size());
            assertEquals(-1, reader.read(ByteBuffer.allocate(1), 0));
        }
    }

    @Test
    void readAtEofReturnsMinusOne() throws IOException {
        byte[] payload = "ACGT".getBytes(StandardCharsets.US_ASCII);
        File bgzf = BgzfTestWriter.writeBgzfFile(tmp.resolve("eof.bgzf"), payload, 64);
        try (BgzfByteReader reader = new BgzfByteReader(bgzf)) {
            assertEquals(-1, reader.read(ByteBuffer.allocate(1), payload.length));
            assertEquals(-1, reader.read(ByteBuffer.allocate(1), payload.length + 100));
        }
    }

    @Test
    void readAfterCloseThrows() throws IOException {
        byte[] payload = "ACGT".getBytes(StandardCharsets.US_ASCII);
        File bgzf = BgzfTestWriter.writeBgzfFile(tmp.resolve("closed.bgzf"), payload, 64);
        BgzfByteReader reader = new BgzfByteReader(bgzf);
        reader.close();
        assertFalse(reader.isOpen());
        assertThrows(IOException.class, () -> reader.read(ByteBuffer.allocate(1), 0));
    }

    @Test
    void negativePositionThrows() throws IOException {
        byte[] payload = "ACGT".getBytes(StandardCharsets.US_ASCII);
        File bgzf = BgzfTestWriter.writeBgzfFile(tmp.resolve("neg.bgzf"), payload, 64);
        try (BgzfByteReader reader = new BgzfByteReader(bgzf)) {
            assertThrows(IllegalArgumentException.class, () -> reader.read(ByteBuffer.allocate(1), -1));
            assertThrows(IllegalArgumentException.class, () -> reader.position(-1));
        }
    }

    @Test
    void lruEvictionWithMoreThan16Blocks() throws IOException {
        // 20 blocks of 10 bytes each = 200 bytes total; cache holds 16, so eviction occurs
        Random rnd = new Random(7);
        byte[] payload = new byte[200];
        for (int i = 0; i < payload.length; i++) payload[i] = (byte) ('A' + rnd.nextInt(5));
        File bgzf = BgzfTestWriter.writeBgzfFile(tmp.resolve("lru.bgzf"), payload, 10);

        try (BgzfByteReader reader = new BgzfByteReader(bgzf)) {
            // forward pass
            byte[] fwd = new byte[payload.length];
            ByteBuffer buf = ByteBuffer.wrap(fwd);
            int n = reader.read(buf, 0);
            assertEquals(payload.length, n);
            assertArrayEquals(payload, fwd);

            // backward random-access (forces cache misses / re-inflation after eviction)
            for (int i = payload.length - 1; i >= 0; i--) {
                ByteBuffer one = ByteBuffer.allocate(1);
                reader.read(one, i);
                assertEquals(payload[i], one.get(0));
            }
        }
    }

    // ---- BgzfDetector edge cases ----

    @Test
    void detectEmptyFile() throws IOException {
        File empty = tmp.resolve("empty.txt").toFile();
        empty.createNewFile();
        assertEquals(BgzfDetector.Compression.UNCOMPRESSED, BgzfDetector.detect(empty));
    }

    @Test
    void detectShortFile() throws IOException {
        Path p = tmp.resolve("one.txt");
        Files.write(p, new byte[]{0x1f});
        assertEquals(BgzfDetector.Compression.UNCOMPRESSED, BgzfDetector.detect(p.toFile()));
    }

    @Test
    void detectFextraWithoutBcIsPlainGzip() throws IOException {
        // Write a minimal gzip with FEXTRA but an extra field that contains no BC subfield
        Path p = tmp.resolve("fextra_no_bc.gz");
        byte[] header = {
            0x1f, (byte)0x8b,  // ID1 ID2
            0x08,              // CM = deflate
            0x04,              // FLG = FEXTRA
            0,0,0,0,           // MTIME
            0, (byte)0xff,     // XFL, OS
            0x04, 0x00,        // XLEN = 4
            0x41, 0x42, 0x02, 0x00  // extra: SI1='A' SI2='B' SLEN=0 (not BC)
        };
        Files.write(p, header);
        assertEquals(BgzfDetector.Compression.PLAIN_GZIP, BgzfDetector.detect(p.toFile()));
    }

    @Test
    void detectBcNotFirstSubfieldIsStillBgzf() throws IOException {
        // Put a non-BC subfield before the BC subfield
        Path p = tmp.resolve("bc_second.gz");
        byte[] header = {
            0x1f, (byte)0x8b, 0x08, 0x04,
            0,0,0,0, 0, (byte)0xff,
            0x0a, 0x00,  // XLEN = 10
            // subfield 1: SI1='X' SI2='Y' SLEN=2 DATA=0000
            0x58, 0x59, 0x02, 0x00, 0x00, 0x00,
            // subfield 2: BC SLEN=2 BSIZE=27
            0x42, 0x43, 0x02, 0x00
            // (BSIZE bytes omitted — detector only checks for subfield presence, not reads BSIZE)
        };
        // detector reads extra and calls hasBgzfSubfield, which scans all subfields
        Files.write(p, header);
        assertEquals(BgzfDetector.Compression.BGZF, BgzfDetector.detect(p.toFile()));
    }

    // ---- Corrupt BGZF error paths ----

    @Test
    void corruptCrcThrowsIoException() throws IOException {
        byte[] payload = "ACGTACGT".getBytes(StandardCharsets.US_ASCII);
        byte[] corrupt = BgzfTestWriter.toBgzfCorruptCrc(payload);
        File f = tmp.resolve("crc.bgzf").toFile();
        Files.write(f.toPath(), corrupt);
        try (BgzfByteReader reader = new BgzfByteReader(f)) {
            IOException e = assertThrows(IOException.class,
                    () -> reader.read(ByteBuffer.allocate(payload.length), 0));
            assertTrue(e.getMessage().toLowerCase().contains("crc32"));
        }
    }

    @Test
    void oversizedIsizeThrowsIoException() throws IOException {
        byte[] payload = "ACGTACGT".getBytes(StandardCharsets.US_ASCII);
        byte[] corrupt = BgzfTestWriter.toBgzfOversizedIsize(payload);
        File f = tmp.resolve("isize.bgzf").toFile();
        Files.write(f.toPath(), corrupt);
        IOException e = assertThrows(IOException.class, () -> new BgzfByteReader(f));
        assertTrue(e.getMessage().contains("ISIZE"));
    }

    @Test
    void shortBcXlenThrowsIoException() throws IOException {
        byte[] corrupt = BgzfTestWriter.toBgzfShortBcXlen();
        File f = tmp.resolve("shortbc.bgzf").toFile();
        Files.write(f.toPath(), corrupt);
        IOException e = assertThrows(IOException.class, () -> new BgzfByteReader(f));
        assertTrue(e.getMessage().toLowerCase().contains("truncated"));
    }

    @Test
    void truncatedFileThrowsIoException() throws IOException {
        byte[] payload = "ACGTACGTACGT".getBytes(StandardCharsets.US_ASCII);
        byte[] corrupt = BgzfTestWriter.toBgzfTruncated(payload, payload.length);
        File f = tmp.resolve("trunc.bgzf").toFile();
        Files.write(f.toPath(), corrupt);
        assertThrows(IOException.class, () -> new BgzfByteReader(f));
    }

    // ---- SeekableByteReaderFactory branches ----

    @Test
    void factoryOpensBgzfFile() throws IOException {
        byte[] payload = "ACGT".getBytes(StandardCharsets.US_ASCII);
        File bgzf = BgzfTestWriter.writeBgzfFile(tmp.resolve("fac.bgzf"), payload, 64);
        try (var reader = SeekableByteReaderFactory.open(bgzf)) {
            assertTrue(reader instanceof BgzfByteReader);
            assertEquals(payload.length, reader.size());
        }
    }

    @Test
    void factoryOpensUncompressedFile() throws IOException {
        Path p = tmp.resolve("plain.txt");
        Files.write(p, "ACGT".getBytes(StandardCharsets.US_ASCII));
        try (var reader = SeekableByteReaderFactory.open(p.toFile())) {
            assertTrue(reader instanceof FileChannelByteReader);
            assertEquals(4, reader.size());
        }
    }

    // ---- FileChannelByteReader ----

    @Test
    void fileChannelByteReaderBasicRead() throws IOException {
        byte[] content = "ACGTNN".getBytes(StandardCharsets.US_ASCII);
        Path p = tmp.resolve("fc.txt");
        Files.write(p, content);
        try (FileChannelByteReader r = new FileChannelByteReader(p.toFile())) {
            assertEquals(content.length, r.size());
            assertTrue(r.isOpen());
            ByteBuffer buf = ByteBuffer.allocate(content.length);
            int n = r.read(buf, 0);
            assertEquals(content.length, n);
            assertArrayEquals(content, buf.array());
        }
    }

    @Test
    void fileChannelByteReaderPositionAndSeek() throws IOException {
        byte[] content = "0123456789".getBytes(StandardCharsets.US_ASCII);
        Path p = tmp.resolve("fcpos.txt");
        Files.write(p, content);
        try (FileChannelByteReader r = new FileChannelByteReader(p.toFile())) {
            r.position(3);
            assertEquals(3, r.position());
            ByteBuffer buf = ByteBuffer.allocate(4);
            int n = r.read(buf, 3);
            assertEquals(4, n);
            buf.flip();
            assertEquals("3456", StandardCharsets.US_ASCII.decode(buf).toString());
        }
    }

    @Test
    void fileChannelByteReaderCloseIsIdempotent() throws IOException {
        Path p = tmp.resolve("fccl.txt");
        Files.write(p, new byte[]{1, 2, 3});
        FileChannelByteReader r = new FileChannelByteReader(p.toFile());
        assertTrue(r.isOpen());
        r.close();
        assertFalse(r.isOpen());
        r.close(); // must not throw
    }
}
