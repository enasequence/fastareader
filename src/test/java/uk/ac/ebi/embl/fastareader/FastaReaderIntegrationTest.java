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
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import org.junit.jupiter.api.Test;
import uk.ac.ebi.embl.fastareader.exception.FastaFileException;

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
            List<FastaEntry> entries = service.getFastaEntries();
            assertEquals(2, entries.size(), "should parse 2 FASTA entries");

            Optional<FastaEntry> entry1 = service.getFastaWithId(0L);
            Optional<FastaEntry> entry2 = service.getFastaWithId(1L);
            assertTrue(entry1.isPresent(), "index for 0L fastaReaderId must exist");
            assertTrue(entry2.isPresent(), "index for 1 fastaReaderId must exist");

            // check header
            assertTrue(entry1.get()
                    .headerLine
                    .equals(">MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken Zażółć gęślą jaźń — 日本語"));
            assertTrue(
                    entry2.get()
                            .headerLine
                            .equals(
                                    "> ID2 | {\"description\":\"这绝对是一段\u200B\u200B描述。\", \"molecule_type\":\"脱氧核糖核酸\", \"topology\":\"linear\"}"));
        }
    }

    @Test
    void proccessingEntriesWithCarriageReturnsCorrectly() throws IOException, FastaFileException {
        File fasta = TestResources.file("fasta", "example_with_carriage_return_char.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<FastaEntry> entries = service.getFastaEntries();
            assertEquals(2, entries.size(), "should parse 2 FASTA entries");

            Set<Long> ids = entries.stream().map(FastaEntry::getFastaReaderId).collect(Collectors.toSet());
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));

            Optional<FastaEntry> entry1 = service.getFastaWithId(0L);
            Optional<FastaEntry> entry2 = service.getFastaWithId(1L);
            Optional<FastaEntry> imaginaryEntry = service.getFastaWithId(2L);
            assertTrue(entry1.isPresent(), "index for AF123456.1 must exist");
            assertTrue(entry2.isPresent(), "index for AF123455.2 must exist");
            assertTrue(imaginaryEntry.isEmpty(), "index for ID3 must not exist");

            // check header
            assertTrue(
                    entry1.get()
                            .headerLine
                            .equals(
                                    ">AF123456.1 |{\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"circular\"}"));
            assertTrue(
                    entry2.get()
                            .headerLine
                            .equals(
                                    ">AF123455.2 |{\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"circular\"}"));

            // From the sample file above:
            assertEquals(9, entry1.get().leadingNsCount, "AF123456.1 leading Ns");
            assertEquals(1, entry1.get().trailingNsCount, "AF123456.1 trailing Ns");
            assertEquals(0, entry2.get().leadingNsCount, "AF123455.2 leading Ns");
            assertEquals(0, entry2.get().trailingNsCount, "AF123455.2 trailing Ns");

            String sequence1StartSlice = service.getSequenceSliceString(
                    entry1.get().fastaReaderId, 1, 11, SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals("CCCGGCGCGGG", sequence1StartSlice);

            String sequence1EndSlice = service.getSequenceSliceString(
                    entry1.get().fastaReaderId,
                    entry1.get().totalBasesWithoutNBases - 9,
                    entry1.get().totalBasesWithoutNBases,
                    SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals("AAAAAAAAAA", sequence1EndSlice);

            String sequence2withoutNbases = service.getSequenceSliceString(
                    entry2.get().fastaReaderId,
                    1,
                    entry2.get().totalBasesWithoutNBases,
                    SequenceRangeOption.WITHOUT_EDGE_N_BASES);
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
    void gettingSequenceSliceAsStringReturnsCorrectly() throws IOException, FastaFileException {
        File fasta = TestResources.file("fasta", "example.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<FastaEntry> entries = service.getFastaEntries();

            assertEquals(2, entries.size(), "should parse 2 FASTA entries");

            Set<Long> ids = entries.stream().map(FastaEntry::getFastaReaderId).collect(Collectors.toSet());
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));

            Optional<FastaEntry> entry1 = service.getFastaWithId(0L);
            Optional<FastaEntry> entry2 = service.getFastaWithId(1L);
            Optional<FastaEntry> imaginaryEntry = service.getFastaWithId(123L);
            assertTrue(entry1.isPresent(), "index for ID1 must exist");
            assertTrue(entry2.isPresent(), "index for ID2 must exist");
            assertTrue(imaginaryEntry.isEmpty(), "index for ID3 must not exist");

            // check header
            assertTrue(entry1.get()
                    .headerLine
                    .equals(">ID1 | {\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"linear\"}"));
            assertTrue(entry2.get()
                    .headerLine
                    .equals(">ID2 | {\"description\":\"x\", \"molecule_type\":\"dna\", \"topology\":\"circular\"}"));

            // From the sample file above:
            assertEquals(2, entry1.get().leadingNsCount, "ID1 leading Ns");
            assertEquals(2, entry1.get().trailingNsCount, "ID1 trailing Ns");
            assertEquals(0, entry2.get().leadingNsCount, "ID2 leading Ns");
            assertEquals(0, entry2.get().trailingNsCount, "ID2 trailing Ns");

            String sequence1 = service.getSequenceSliceString(
                    entry1.get().fastaReaderId, 1, entry1.get().totalBases, SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals("NNACACGTTTNN", sequence1);

            String sequence2 = service.getSequenceSliceString(
                    entry2.get().fastaReaderId, 1, entry2.get().totalBases, SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals("ACGTGGGG", sequence2);

            String sequence1withoutNbases = service.getSequenceSliceString(
                    entry1.get().fastaReaderId,
                    1,
                    entry1.get().totalBasesWithoutNBases,
                    SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals("ACACGTTT", sequence1withoutNbases);
        }
    }

    @Test
    void gettingSequenceViaReaderGivesCorrectResult() throws IOException, FastaFileException {
        File fasta = TestResources.file("fasta", "example.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<FastaEntry> entries = service.getFastaEntries();
            assertEquals(2, entries.size(), "should parse 2 FASTA entries");

            Set<Long> ids = entries.stream().map(FastaEntry::getFastaReaderId).collect(Collectors.toSet());
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));
            Optional<FastaEntry> entry1 = service.getFastaWithId(0L);
            Optional<FastaEntry> entry2 = service.getFastaWithId(1L);

            // stream whole sequence with the reader
            String streamedSequence;
            try (java.io.Reader r = service.getSequenceSliceReader(
                    entry1.get().fastaReaderId, 1, entry1.get().totalBases, SequenceRangeOption.WHOLE_SEQUENCE)) {
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
                    entry1.get().fastaReaderId,
                    1,
                    entry1.get().totalBasesWithoutNBases,
                    SequenceRangeOption.WITHOUT_EDGE_N_BASES)) {
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
            try (java.io.Reader r = service.getSequenceSliceReader(
                    entry2.get().fastaReaderId, 1, entry2.get().totalBases, SequenceRangeOption.WHOLE_SEQUENCE)) {
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
    void gettingStringAsAStringVsStreamProducesSameResultSlices() throws IOException, FastaFileException {
        File fasta = TestResources.file("fasta", "example.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<FastaEntry> entries = service.getFastaEntries();
            assertEquals(2, entries.size(), "should parse 2 FASTA entries");

            Set<Long> ids = entries.stream().map(e -> e.getFastaReaderId()).collect(Collectors.toSet());
            assertTrue(ids.contains(0L));
            assertTrue(ids.contains(1L));
            Optional<FastaEntry> entry1 = service.getFastaWithId(0L);
            Optional<FastaEntry> entry2 = service.getFastaWithId(1L);

            for (long end = 2; end <= entry1.get().totalBases; end++) {
                // get slice as string
                String sequence = service.getSequenceSliceString(
                        entry1.get().fastaReaderId, 1, end, SequenceRangeOption.WHOLE_SEQUENCE);
                // stream sequence with the reader
                String streamedSequence;
                try (java.io.Reader r = service.getSequenceSliceReader(
                        entry1.get().fastaReaderId, 1, end, SequenceRangeOption.WHOLE_SEQUENCE)) {
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

            for (long end = 2; end <= entry2.get().totalBases; end++) {
                // get slice as string
                String sequence2 = service.getSequenceSliceString(
                        entry2.get().fastaReaderId, 1, end, SequenceRangeOption.WHOLE_SEQUENCE);
                // stream sequence with the reader
                String streamedSequence2;
                try (java.io.Reader r = service.getSequenceSliceReader(
                        entry2.get().fastaReaderId, 1, end, SequenceRangeOption.WHOLE_SEQUENCE)) {
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
    void readBigSequenceSuccessfully() throws IOException, FastaFileException {
        File fasta = TestResources.file("fasta", "single_fasta_large_sequence.txt");
        try (FastaReader service = new FastaReader(fasta)) {

            List<FastaEntry> entries = service.getFastaEntries();
            assertEquals(1, entries.size(), "should parse 1 FASTA entry");

            Set<Long> ids = entries.stream().map(FastaEntry::getFastaReaderId).collect(Collectors.toSet());
            assertTrue(ids.contains(0));
            Optional<FastaEntry> entry1 = service.getFastaWithId(0L);

            // get first 16 chars
            String sequenceStart = service.getSequenceSliceString(
                    entry1.get().fastaReaderId, 1, 16, SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals(sequenceStart, "GGGCTTTAAATGGCTC");

            // get last 16 chars
            String sequenceEnd = service.getSequenceSliceString(
                    entry1.get().fastaReaderId,
                    entry1.get().totalBases - 15,
                    entry1.get().totalBases,
                    SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals(sequenceEnd, "GAATTCTGATGGCTGT");
        }
    }
}
