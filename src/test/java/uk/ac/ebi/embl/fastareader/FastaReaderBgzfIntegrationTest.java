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

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.zip.GZIPOutputStream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import uk.ac.ebi.embl.fastareader.api.SequenceFormatReader;
import uk.ac.ebi.embl.fastareader.api.SequenceFormatReaderFactory;
import uk.ac.ebi.embl.fastareader.api.rereading.SequenceInfoDTO;
import uk.ac.ebi.embl.fastareader.exception.FastaFileException;
import uk.ac.ebi.embl.fastareader.io.BgzfTestWriter;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;

/** Round-trip equivalence: a BGZF copy of each FASTA fixture must read identically. */
class FastaReaderBgzfIntegrationTest {

    @TempDir
    Path tmp;

    @Test
    void bgzfReadsIdenticallyToUncompressed_example() throws Exception {
        assertFastaRoundTrip("example.txt", 16);
    }

    @Test
    void bgzfReadsIdenticallyToUncompressed_carriageReturn() throws Exception {
        assertFastaRoundTrip("example_with_carriage_return_char.txt", 16);
    }

    @Test
    void bgzfReadsIdenticallyToUncompressed_differingHeaders() throws Exception {
        assertFastaRoundTrip("differing_headerline_fasta.txt", 32);
    }

    @Test
    void bgzfReadsIdenticallyToUncompressed_onlyNs() throws Exception {
        assertFastaRoundTrip("only_ns_example.txt", 16);
    }

    /**
     * Partial-file reading over BGZF: a BGZF copy of a GFF3 file (leading annotation lines followed by an
     * embedded FASTA section) read from an explicit {@code startByteOffset} must yield the same entries as
     * the uncompressed file read from the same offset. The offset is interpreted against the decompressed
     * byte space, and the tiny block size forces it to land inside / across BGZF blocks.
     */
    @Test
    void bgzfPartialReadFromByteOffsetMatchesUncompressed() throws Exception {
        String fixture = "gff3_with_embedded_fasta.txt";
        byte[] raw = Files.readAllBytes(TestResources.file("fasta", fixture).toPath());
        long fastaOffset = indexOf(raw, ">ID1".getBytes(StandardCharsets.UTF_8));
        assertTrue(fastaOffset > 0, "embedded FASTA must not start at byte 0");

        File bgzf = BgzfFixtures.bgzfCopyOf("fasta", fixture, tmp, 7);

        SequenceAlphabet alphabet = SequenceAlphabet.defaultNucleotideAlphabet();
        try (FastaReader plain = new FastaReader(TestResources.file("fasta", fixture), alphabet, fastaOffset);
                FastaReader compressed = new FastaReader(bgzf, alphabet, fastaOffset)) {

            assertEquals(List.of(0L, 1L), compressed.getOrderedIds(), "ordered ids");
            assertEquals(plain.getOrderedIds(), compressed.getOrderedIds(), "ordered ids parity");

            for (long id : plain.getOrderedIds()) {
                assertEquals(plain.getHeaderline(id), compressed.getHeaderline(id), "header " + id);
                assertEquals(plain.getSequenceIndex(id), compressed.getSequenceIndex(id), "index " + id);
                assertEquals(plain.getStats(id), compressed.getStats(id), "stats " + id);
                assertEquals(plain.getGapRegions(id), compressed.getGapRegions(id), "gaps " + id);

                long total = plain.getStats(id).totalBases();
                String expected = plain.getSequenceSlice(id, 1, total, SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(
                        expected,
                        compressed.getSequenceSlice(id, 1, total, SequenceRangeOption.WHOLE_SEQUENCE),
                        "slice string " + id);
                assertEquals(expected, streamWhole(compressed, id, total), "streamed slice " + id);
            }
        }
    }

    private static long indexOf(byte[] haystack, byte[] needle) {
        outer:
        for (int i = 0; i <= haystack.length - needle.length; i++) {
            for (int j = 0; j < needle.length; j++) {
                if (haystack[i + j] != needle[j]) continue outer;
            }
            return i;
        }
        return -1;
    }

    private void assertFastaRoundTrip(String fixture, int blockSize) throws Exception {
        File bgzf = BgzfFixtures.bgzfCopyOf("fasta", fixture, tmp, blockSize);

        try (FastaReader plain = new FastaReader(TestResources.file("fasta", fixture));
                FastaReader compressed = new FastaReader(bgzf)) {

            assertEquals(plain.getOrderedIds(), compressed.getOrderedIds(), "ordered ids");

            for (long id : plain.getOrderedIds()) {
                assertEquals(plain.getHeaderline(id), compressed.getHeaderline(id), "header " + id);
                assertEquals(plain.getSequenceIndex(id), compressed.getSequenceIndex(id), "index " + id);
                assertEquals(plain.getStats(id), compressed.getStats(id), "stats " + id);
                assertEquals(plain.getGapRegions(id), compressed.getGapRegions(id), "gaps " + id);

                long total = plain.getStats(id).totalBases();
                String expected = plain.getSequenceSlice(id, 1, total, SequenceRangeOption.WHOLE_SEQUENCE);
                String actual = compressed.getSequenceSlice(id, 1, total, SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(expected, actual, "slice string " + id);

                // streamed slice parity
                String streamed = streamWhole(compressed, id, total);
                assertEquals(expected, streamed, "streamed slice " + id);
            }
        }
    }

    private static String streamWhole(FastaReader reader, long id, long total) throws Exception {
        try (java.io.Reader r = reader.getSequenceSliceReader(id, 1, total, SequenceRangeOption.WHOLE_SEQUENCE)) {
            StringBuilder sb = new StringBuilder();
            char[] cbuf = new char[8192];
            int n;
            while ((n = r.read(cbuf)) != -1) {
                sb.append(cbuf, 0, n);
            }
            return sb.toString();
        }
    }

    // ---- UTF-8 gate ----

    /**
     * A BGZF file whose decompressed content contains non-UTF-8 bytes must be rejected by
     * {@code FastaReader} with a {@link FastaFileException} referencing UTF-8.
     */
    @Test
    void bgzfWithNonUtf8DecompressedContentIsRejected() throws Exception {
        // 0x80 is an invalid UTF-8 start byte — will fail Utf8Detector
        byte[] nonUtf8 = {0x3E, 'h', 'd', 'r', 0x0A, 0x41, (byte) 0x80, 0x43, 0x0A};
        byte[] bgzfBytes = BgzfTestWriter.toBgzf(nonUtf8, 16);
        Path target = tmp.resolve("non_utf8_fasta.fa.bgzf");
        Files.write(target, bgzfBytes);

        FastaFileException ex = assertThrows(FastaFileException.class, () -> new FastaReader(target.toFile()).close());
        assertTrue(
                ex.getMessage().contains("UTF-8"), "exception message should mention UTF-8, got: " + ex.getMessage());
    }

    @Test
    void crossBlockSliceStringMatchesStream() throws Exception {
        // tiny blocks force a single sequence line and slice to span multiple BGZF blocks
        File bgzf = BgzfFixtures.bgzfCopyOf("fasta", "example.txt", tmp, 4);
        try (FastaReader compressed = new FastaReader(bgzf)) {
            List<Long> ids = compressed.getOrderedIds();
            for (long id : ids) {
                long total = compressed.getStats(id).totalBases();
                for (long end = 1; end <= total; end++) {
                    String s = compressed.getSequenceSlice(id, 1, end, SequenceRangeOption.WHOLE_SEQUENCE);
                    String streamed = streamWhole(compressed, id, end);
                    assertEquals(s, streamed, "id " + id + " end " + end);
                }
            }
        }
    }

    // ---- plain-gzip rejection ----

    /**
     * A plain (non-BGZF) gzip file must be rejected by {@code FastaReader} with a
     * {@link FastaFileException} that mentions "bgzip" or "gzip" in its message.
     */
    @Test
    void plainGzipIsRejectedByFastaReader() throws Exception {
        Path gz = tmp.resolve("plain.fasta.gz");
        try (GZIPOutputStream gos = new GZIPOutputStream(Files.newOutputStream(gz))) {
            gos.write(">hdr\nACGT\n".getBytes(StandardCharsets.US_ASCII));
        }
        FastaFileException ex = assertThrows(FastaFileException.class, () -> new FastaReader(gz.toFile()).close());
        String msg = ex.getMessage().toLowerCase();
        assertTrue(
                msg.contains("bgzip") || msg.contains("gzip"),
                "exception message should mention gzip, got: " + ex.getMessage());
    }

    // ---- SequenceInfoDTO re-reading round-trip over BGZF ----

    /**
     * Export {@link SequenceInfoDTO} from a BGZF-backed {@code FastaReader}, reload via
     * {@code SequenceFormatReaderFactory.readBySequenceInfo}, and assert that slices, stats,
     * gap regions, and header lines are identical.
     */
    @Test
    void rereadingViaSequenceInfoDtoMatchesOriginal() throws Exception {
        File bgzf = BgzfFixtures.bgzfCopyOf("fasta", "example.txt", tmp, 16);

        SequenceInfoDTO dto;
        try (SequenceFormatReader service = SequenceFormatReaderFactory.readFasta(bgzf)) {
            dto = service.exportReaderSettings();
        }

        try (SequenceFormatReader original = SequenceFormatReaderFactory.readFasta(bgzf);
                SequenceFormatReader reread = SequenceFormatReaderFactory.readBySequenceInfo(dto)) {

            assertEquals(original.getOrderedIds(), reread.getOrderedIds(), "ordered ids");
            for (long id : original.getOrderedIds()) {
                assertEquals(original.getSequenceIndex(id), reread.getSequenceIndex(id), "index " + id);
                assertEquals(original.getStats(id), reread.getStats(id), "stats " + id);
                assertEquals(original.getGapRegions(id), reread.getGapRegions(id), "gaps " + id);
                assertEquals(original.getHeaderline(id), reread.getHeaderline(id), "header " + id);

                long total = original.getStats(id).totalBases();
                String expectedSlice = original.getSequenceSlice(id, 1, total, SequenceRangeOption.WHOLE_SEQUENCE);
                String rereadSlice = reread.getSequenceSlice(id, 1, total, SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(expectedSlice, rereadSlice, "slice " + id);
            }
        }
    }

    /**
     * Same as {@link #rereadingViaSequenceInfoDtoMatchesOriginal()} but with a multi-entry FASTA
     * and tiny BGZF blocks so header seeking spans block boundaries.
     */
    @Test
    void rereadingMultiEntryBgzfViaSequenceInfoDtoMatchesOriginal() throws Exception {
        // Use tiny blocks so the header line of the second entry spans a block boundary
        File bgzf = BgzfFixtures.bgzfCopyOf("fasta", "differing_headerline_fasta.txt", tmp, 8);

        SequenceInfoDTO dto;
        try (SequenceFormatReader service = SequenceFormatReaderFactory.readFasta(bgzf)) {
            dto = service.exportReaderSettings();
        }

        try (SequenceFormatReader original = SequenceFormatReaderFactory.readFasta(bgzf);
                SequenceFormatReader reread = SequenceFormatReaderFactory.readBySequenceInfo(dto)) {

            assertEquals(original.getOrderedIds(), reread.getOrderedIds(), "ordered ids");
            for (long id : original.getOrderedIds()) {
                assertEquals(original.getSequenceIndex(id), reread.getSequenceIndex(id), "index " + id);
                assertEquals(original.getStats(id), reread.getStats(id), "stats " + id);
                assertEquals(original.getGapRegions(id), reread.getGapRegions(id), "gaps " + id);
                assertEquals(original.getHeaderline(id), reread.getHeaderline(id), "header " + id);

                long total = original.getStats(id).totalBases();
                String expectedSlice = original.getSequenceSlice(id, 1, total, SequenceRangeOption.WHOLE_SEQUENCE);
                String rereadSlice = reread.getSequenceSlice(id, 1, total, SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(expectedSlice, rereadSlice, "slice " + id);
            }
        }
    }
}
