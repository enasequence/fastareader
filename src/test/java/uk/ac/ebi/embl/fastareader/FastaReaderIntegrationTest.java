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
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;
import org.junit.jupiter.api.Test;
import uk.ac.ebi.embl.fastareader.exception.FastaFileException;
import uk.ac.ebi.embl.fastareader.sequenceutils.GapRegion;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;

class FastaReaderIntegrationTest {

    @Test
    void doesNotTolerateImproperHeaders() throws IOException {
        // by improper headers, i mean ones that do not start in > and end with a newline as per the general spec */
        File fasta = TestResources.file("fasta", "malformed_headerline_fasta.txt");

        assertThrows(FastaFileException.class, () -> new FastaReader(fasta));
    }

    @Test
    void toleratesDifferentFastaHeaders() throws IOException {
        File fasta =
                TestResources.file("fasta", "differing_headerline_fasta.txt"); // this example has varying fasta headers

        assertDoesNotThrow(() -> new FastaReader(fasta));
    }

    @Test
    void readsUnicodeHeadersCorrectly() throws IOException, FastaFileException {
        File fasta =
                TestResources.file("fasta", "differing_headerline_fasta.txt"); // this example has varying fasta headers

        try (FastaReader service = new FastaReader(fasta)) {

            List<Long> ids = service.getOrderedIds();
            assertEquals(2, ids.size(), "should parse 2 FASTA entries");
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));

            SequenceStats entry1 = service.getStats(0L);
            SequenceStats entry2 = service.getStats(1L);

            SequenceStats imaginaryEntry = service.getStats(123L);
            assertNotNull(entry1, "index for ID1 must exist");
            assertNotNull(entry2, "index for ID2 must exist");
            assertNull(imaginaryEntry, "index for ID3 must not exist");

            // check header
            var headerLine = service.getHeaderline(0L);
            var headerLine1 = service.getHeaderline(1L);
            assertTrue(headerLine.isPresent(), "header line for ID1 must exist");
            assertTrue(headerLine1.isPresent(), "header line for ID1 must exist");
            assertEquals(
                    ">MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken Zażółć gęślą jaźń — 日本語",
                    headerLine.get());
            assertEquals(
                    "> ID2 | {\"description\":\"这绝对是一段\u200B\u200B描述。\", \"molecule_type\":\"脱氧核糖核酸\", \"topology\":\"linear\"}",
                    headerLine1.get());
        }
    }

    @Test
    void proccessingEntriesWithCarriageReturnsCorrectly() throws Exception {
        File fasta = TestResources.file("fasta", "example_with_carriage_return_char.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<Long> ids = service.getOrderedIds();
            assertEquals(2, ids.size(), "should parse 2 FASTA entries");
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));

            SequenceStats entry1 = service.getStats(0L);
            SequenceStats entry2 = service.getStats(1L);

            SequenceStats imaginaryEntry = service.getStats(123L);
            assertNotNull(entry1, "index for AF123456.1 must exist");
            assertNotNull(entry2, "index for AF123455.2 must exist");
            assertNull(imaginaryEntry, "index for ID3 must not exist");

            // check headerlines
            // check header
            var headerLine = service.getHeaderline(0L);
            var headerLine1 = service.getHeaderline(1L);
            assertTrue(headerLine.isPresent(), "header line for ID1 must exist");
            assertTrue(headerLine1.isPresent(), "header line for ID1 must exist");
            assertEquals(
                    ">AF123456.1 |{\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"circular\"}",
                    headerLine.get());
            assertEquals(
                    ">AF123455.2 |{\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"circular\"}",
                    headerLine1.get());

            // From the sample file above:
            assertEquals(9, entry1.leadingNsCount(), "AF123456.1 leading Ns");
            assertEquals(1, entry1.trailingNsCount(), "AF123456.1 trailing Ns");
            assertEquals(0, entry2.leadingNsCount(), "AF123455.2 leading Ns");
            assertEquals(0, entry2.trailingNsCount(), "AF123455.2 trailing Ns");

            String sequence1StartSlice = service.getSequenceSlice(0L, 1, 11, SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals("CCCGGCGCGGG", sequence1StartSlice);

            String sequence1EndSlice = service.getSequenceSlice(
                    0L,
                    entry1.totalBasesWithoutNBases() - 9,
                    entry1.totalBasesWithoutNBases(),
                    SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals("AAAAAAAAAA", sequence1EndSlice);

            String sequence2withoutNbases = service.getSequenceSlice(
                    1L, 1, entry2.totalBasesWithoutNBases(), SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals(
                    "CCCGGCGCGGGCAAGAAGCTGCCGCGTCTGCCCAAGTGTGCCCGCTGCCGCAACCACGGC"
                            + "TACTCCTCGCCGCTGAAGGGGCACAAGCGGTTCTGCATGTGGCGGGACTGCCAGTGCAAG"
                            + "AAGTGCAGCCTGATCCGCCGAGCGGCAGGGGTGATGGCCGTGCAGGTTGCACTGAGGAGG"
                            + "ATGTGTTTGTAGTGGTTCCTCGTAGGCTCCAGACGTTTTCTCCTCGTATCGCCAAATTAA"
                            + "CGCGTTTTGTAGTGGTTCCTCGTAGGCTCCAGACGTTTTCTCCTCAGACGTGGCCAGCAA"
                            + "ACAAGTCTCAAAAAAAAGTTACGTGCGTTTCTGCGAGTGTTATTTTGTTAAGAACGGCTC"
                            + "ACAGTGTCCTCTTCCTGTGTTACAGAAGCCAACCTGAAATGAAACTAGTCTGGAAAAATT"
                            + "CATTGTTCTCTGTAGTTGCAGCTGTACCTGAAATAAAAATGTTATTGATGACTGAAAAAA"
                            + "AAAAAAAAAAAA",
                    sequence2withoutNbases);
        }
    }

    @Test
    void readsEmbeddedFastaFromExplicitByteOffset() throws Exception {
        // A GFF3 file: leading annotation lines followed by an embedded FASTA section whose bytes are
        // identical to fasta/example.txt. Starting the scan at the offset of the first '>' header must
        // yield exactly the same entries as reading example.txt from byte 0.
        File gff3 = TestResources.file("fasta", "gff3_with_embedded_fasta.txt");
        byte[] bytes = Files.readAllBytes(gff3.toPath());
        long fastaOffset = indexOf(bytes, ">ID1".getBytes(StandardCharsets.UTF_8));
        assertTrue(fastaOffset > 0, "embedded FASTA must not start at byte 0");

        try (FastaReader service = new FastaReader(gff3, SequenceAlphabet.defaultNucleotideAlphabet(), fastaOffset)) {
            assertEquals(List.of(0L, 1L), service.getOrderedIds());

            assertEquals(
                    ">ID1 | {\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"linear\"}",
                    service.getHeaderline(0L).orElseThrow());
            assertEquals(
                    ">ID2 | {\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"circular\"}",
                    service.getHeaderline(1L).orElseThrow());

            SequenceStats entry1 = service.getStats(0L);
            assertEquals(2, entry1.leadingNsCount(), "ID1 leading Ns");
            assertEquals(2, entry1.trailingNsCount(), "ID1 trailing Ns");

            assertEquals(
                    "NNACACGTTTNN",
                    service.getSequenceSlice(0L, 1, entry1.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE));
            assertEquals(List.of(new GapRegion(1, 2), new GapRegion(11, 12)), service.getGapRegions(0L));
            assertEquals(
                    "ACGTGGGG",
                    service.getSequenceSlice(
                            1L, 1, service.getStats(1L).totalBases(), SequenceRangeOption.WHOLE_SEQUENCE));
        }
    }

    @Test
    void tolerantOfOffsetLandingOnNonHeaderPrefixLine() throws Exception {
        // An offset that lands on the '##FASTA' directive line (before the first '>' header) is tolerated:
        // preceding non-header lines are skipped and scanning still finds the FASTA entries.
        File gff3 = TestResources.file("fasta", "gff3_with_embedded_fasta.txt");
        byte[] bytes = Files.readAllBytes(gff3.toPath());
        long directiveOffset = indexOf(bytes, "##FASTA".getBytes(StandardCharsets.UTF_8));
        assertTrue(directiveOffset > 0, "##FASTA directive must exist past byte 0");

        try (FastaReader service =
                new FastaReader(gff3, SequenceAlphabet.defaultNucleotideAlphabet(), directiveOffset)) {
            assertEquals(List.of(0L, 1L), service.getOrderedIds());
            assertEquals(
                    "NNACACGTTTNN",
                    service.getSequenceSlice(
                            0L, 1, service.getStats(0L).totalBases(), SequenceRangeOption.WHOLE_SEQUENCE));
        }
    }

    @Test
    void rejectsOffsetBeyondFileSize() throws IOException {
        File gff3 = TestResources.file("fasta", "gff3_with_embedded_fasta.txt");
        long tooBig = gff3.length() + 1;
        assertThrows(
                FastaFileException.class,
                () -> new FastaReader(gff3, SequenceAlphabet.defaultNucleotideAlphabet(), tooBig));
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

    @Test
    void gettingSequenceSliceAsStringReturnsCorrectly() throws Exception {
        File fasta = TestResources.file("fasta", "example.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<Long> ids = service.getOrderedIds();
            assertEquals(2, ids.size(), "should parse 2 FASTA entries");
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));

            SequenceStats entry1 = service.getStats(0L);
            SequenceStats entry2 = service.getStats(1L);

            SequenceStats imaginaryEntry = service.getStats(123L);
            assertNotNull(entry1, "index for ID1 must exist");
            assertNotNull(entry2, "index for ID2 must exist");
            assertNull(imaginaryEntry, "index for ID3 must not exist");

            // check header
            var headerLine = service.getHeaderline(0L);
            var headerLine1 = service.getHeaderline(1L);
            assertTrue(headerLine.isPresent(), "header line for ID1 must exist");
            assertTrue(headerLine1.isPresent(), "header line for ID1 must exist");
            assertEquals(
                    ">ID1 | {\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"linear\"}",
                    headerLine.get());
            assertEquals(
                    ">ID2 | {\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"circular\"}",
                    headerLine1.get());

            // From the sample file above:
            assertEquals(2, entry1.leadingNsCount(), "ID1 leading Ns");
            assertEquals(2, entry1.trailingNsCount(), "ID1 trailing Ns");
            assertEquals(0, entry2.leadingNsCount(), "ID2 leading Ns");
            assertEquals(0, entry2.trailingNsCount(), "ID2 trailing Ns");

            String sequence1 = service.getSequenceSlice(0L, 1, entry1.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals("NNACACGTTTNN", sequence1);
            assertEquals(List.of(new GapRegion(1, 2), new GapRegion(11, 12)), service.getGapRegions(0L));
            assertEquals(List.of(new GapRegion(11, 12)), service.getGapRegions(0L, 3, 12));

            String sequence2 = service.getSequenceSlice(1L, 1, entry2.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals("ACGTGGGG", sequence2);

            String sequence1withoutNbases = service.getSequenceSlice(
                    0L, 1, entry1.totalBasesWithoutNBases(), SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals("ACACGTTT", sequence1withoutNbases);
        }
    }

    @Test
    void gettingSequenceViaReaderGivesCorrectResult() throws Exception {
        File fasta = TestResources.file("fasta", "example.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<Long> ids = service.getOrderedIds();
            assertEquals(2, ids.size(), "should parse 2 FASTA entries");
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));

            SequenceStats entry1 = service.getStats(0L);
            SequenceStats entry2 = service.getStats(1L);

            // stream whole sequence with the reader
            String streamedSequence;
            try (java.io.Reader r =
                    service.getSequenceSliceReader(0L, 1, entry1.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE)) {
                StringBuilder sb = new StringBuilder();
                char[] cbuf = new char[8192];
                int n;
                while ((n = r.read(cbuf)) != -1) {
                    sb.append(cbuf, 0, n);
                }
                streamedSequence = sb.toString();
            }
            // compare
            assertEquals("NNACACGTTTNN", streamedSequence);

            // stream whole sequence with the reader
            String streamedSequenceWithoutNbases;
            try (java.io.Reader r = service.getSequenceSliceReader(
                    0L, 1, entry1.totalBasesWithoutNBases(), SequenceRangeOption.WITHOUT_EDGE_N_BASES)) {
                StringBuilder sb = new StringBuilder();
                char[] cbuf = new char[8192];
                int n;
                while ((n = r.read(cbuf)) != -1) {
                    sb.append(cbuf, 0, n);
                }
                streamedSequenceWithoutNbases = sb.toString();
            }
            // compare
            assertEquals("ACACGTTT", streamedSequenceWithoutNbases);

            // stream sequence with the reader
            String streamedSequence2;
            try (java.io.Reader r =
                    service.getSequenceSliceReader(1L, 1, entry2.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE)) {
                StringBuilder sb = new StringBuilder();
                char[] cbuf = new char[8192];
                int n;
                while ((n = r.read(cbuf)) != -1) {
                    sb.append(cbuf, 0, n);
                }
                streamedSequence2 = sb.toString();
            }
            // compare
            assertEquals("ACGTGGGG", streamedSequence2);
        }
    }

    @Test
    void gettingStringAsAStringVsStreamProducesSameResultSlices() throws Exception {
        File fasta = TestResources.file("fasta", "example.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<Long> ids = service.getOrderedIds();
            assertEquals(2, ids.size(), "should parse 2 FASTA entries");
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));

            SequenceStats entry1 = service.getStats(0L);
            SequenceStats entry2 = service.getStats(1L);

            for (long end = 2; end <= entry1.totalBases(); end++) {
                // get slice as string
                String sequence = service.getSequenceSlice(0L, 1, end, SequenceRangeOption.WHOLE_SEQUENCE);
                // stream sequence with the reader
                String streamedSequence;
                try (java.io.Reader r =
                        service.getSequenceSliceReader(0L, 1, end, SequenceRangeOption.WHOLE_SEQUENCE)) {
                    StringBuilder sb = new StringBuilder();
                    char[] cbuf = new char[8192];
                    int n;
                    while ((n = r.read(cbuf)) != -1) {
                        sb.append(cbuf, 0, n);
                    }
                    streamedSequence = sb.toString();
                }
                // compare
                assertEquals(sequence, streamedSequence);
            }

            for (long end = 2; end <= entry2.totalBases(); end++) {
                // get slice as string
                String sequence2 = service.getSequenceSlice(1L, 1, end, SequenceRangeOption.WHOLE_SEQUENCE);
                // stream sequence with the reader
                String streamedSequence2;
                try (java.io.Reader r =
                        service.getSequenceSliceReader(1L, 1, end, SequenceRangeOption.WHOLE_SEQUENCE)) {
                    StringBuilder sb = new StringBuilder();
                    char[] cbuf = new char[8192];
                    int n;
                    while ((n = r.read(cbuf)) != -1) {
                        sb.append(cbuf, 0, n);
                    }
                    streamedSequence2 = sb.toString();
                }
                // compare
                assertEquals(sequence2, streamedSequence2);
            }
        }
    }

    // to run this, curl the sequence with: curl -o single_fasta_large_sequence.txt
    // https://www.ebi.ac.uk/ena/cram/md5/11398cc4b68f2cceb4fd50b742d4b1ec
    // then to add the fasta header run something like :
    //
    // tmp="$(mktemp "${TMPDIR:-/tmp}/prepend.XXXXXX")" &&
    // { printf '%s\n' '>ID1 | {"description":"x", "molecule_type":"dna", "topology":"linear"}'; cat --
    // single_fasta_large_sequence.txt; } >"$tmp" &&
    // mv -f -- "$tmp" single_fasta_large_sequence.txt
    //
    // then just move the fasta into whatever/gff3tools/src/test/resources/fasta/
    // and run the test
    // @Test
    void readBigSequenceSuccessfully() throws Exception {
        File fasta = TestResources.file("fasta", "single_fasta_large_sequence.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<Long> ids = service.getOrderedIds();
            assertEquals(1, ids.size(), "should parse 1 FASTA entry");
            assertTrue(ids.contains(0L));
            assertTrue(!ids.contains(1L));

            SequenceStats entry1 = service.getStats(0L);

            // get first 16 chars
            String sequenceStart = service.getSequenceSlice(0L, 1, 16, SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals(sequenceStart, "GGGCTTTAAATGGCTC");

            // get last 16 chars
            String sequenceEnd = service.getSequenceSlice(
                    1L, entry1.totalBases() - 15, entry1.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals(sequenceEnd, "GAATTCTGATGGCTGT");
        }
    }

    @Test
    void parsesMultiEntryFastaWithCrOnlyLineEndings() throws Exception {
        File fasta = TestResources.file("fasta", "two_entries_cr_only.txt");
        try (FastaReader service = new FastaReader(fasta)) {
            List<Long> ids = service.getOrderedIds();
            assertEquals(2, ids.size(), "CR-only FASTA must yield 2 entries");
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));

            assertTrue(service.getHeaderline(0L).isPresent());
            assertTrue(service.getHeaderline(1L).isPresent());

            String seq1 = service.getSequenceSlice(0L, 1, 4, SequenceRangeOption.WHOLE_SEQUENCE);
            String seq2 = service.getSequenceSlice(1L, 1, 4, SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals("ACGT", seq1);
            assertEquals("GGGG", seq2);
        }
    }
}
