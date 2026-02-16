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
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.File;
import java.io.IOException;
import org.junit.jupiter.api.Test;
import uk.ac.ebi.embl.fastareader.exception.SequenceFileException;

public class PlainSequenceReaderIntegrationTest {

    @Test
    void doesNotTolerateAnythingExceptSequenceAlphabetInFile() throws IOException {
        File sequenceFile =
                TestResources.file("sequence", "malformed_example.txt"); // this example has varying fasta headers

        assertThrows(SequenceFileException.class, () -> new SequenceReader(sequenceFile));
    }

    @Test
    void countsOnlyNsSequenceCorrectly() throws IOException, SequenceFileException {
        File sequenceFile = TestResources.file("sequence", "only_ns_example.txt");
        try (SequenceReader service = new SequenceReader(sequenceFile)) {

            SequenceEntry sequence = service.getSequenceInfo();
            assertEquals(43, sequence.leadingNsCount, "leading Ns");
            assertEquals(43, sequence.trailingNsCount, "trailing Ns");
            assertEquals(43, sequence.baseCount.get('N'), "total Ns");
        }
    }

    @Test
    void proccessingSequenceWithCarriageReturnsCorrectly() throws IOException, SequenceFileException {
        File sequenceFile = TestResources.file("sequence", "example_with_carriage_return_char.txt");
        try (SequenceReader service = new SequenceReader(sequenceFile)) {

            SequenceEntry sequence = service.getSequenceInfo();

            // From the sample file above:
            assertEquals(9, sequence.leadingNsCount, "leading Ns");
            assertEquals(1, sequence.trailingNsCount, "trailing Ns");

            String sequence1StartSlice =
                    service.getSequenceSliceString(1, 11, SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals("CCCGGCGCGGG", sequence1StartSlice);

            String sequence1EndSlice = service.getSequenceSliceString(
                    sequence.totalBasesWithoutNBases - 9,
                    sequence.totalBasesWithoutNBases,
                    SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals("AAAAAAAAAA", sequence1EndSlice);

            String sequence1withoutNbases = service.getSequenceSliceString(
                    1, sequence.totalBasesWithoutNBases, SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals(
                    "CCCGGCGCGGGCAAGAAGCTGCCGCGTCTGCCCAAGTGTGCCCGCTGCCGCAACCACGGC"
                            + "TACTCCTCGCCGCTGAAGGGGCACAAGCGGTTCTGCATGTGGCGGGACTGCCAGTGCAAG"
                            + "AAGTGCAGCCTGATCCGCCGAGCGGCAGGGGTGATGGCCGTGCAGGTTGCACTGAGGAGG"
                            + "ATGTGTTTGTAGTGGTTCCTCGTAGGCTCCAGACGTTTTCTCCTCGTATCGCCAAATTAA"
                            + "CGCGTTTTGCATATTACAGTTGAGTGCCTCGACTTAGATTGCAATATAAGCGGCCAGCAA"
                            + "ACAAGTCTCAAAAAAAAGTTACGTGCGTTTCTGCGAGTGTTATTTTGTTAAGAACGGCTC"
                            + "ACAGTGTCCTCTTCCTGTGTTACAGAAGCCAACCTGAAATGAAACTAGTCTGGAAAAATT"
                            + "CATTGTTCTCTGTAGTTGCAGCTGTACCTGAAATAAAAATGTTATTGATGACTGAAAAAA"
                            + "AAAAAAAAAAAA",
                    sequence1withoutNbases);
        }
    }

    @Test
    void gettingSequenceSliceAsStringReturnsCorrectly() throws IOException, SequenceFileException {
        File sequenceFile = TestResources.file("sequence", "example.txt");
        try (SequenceReader service = new SequenceReader(sequenceFile)) {

            SequenceEntry sequence = service.getSequenceInfo();

            // From the sample file above:
            assertEquals(2, sequence.leadingNsCount, "ID1 leading Ns");
            assertEquals(2, sequence.trailingNsCount, "ID1 trailing Ns");

            String sequence1 =
                    service.getSequenceSliceString(1, sequence.totalBases, SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals("NNACACGTTTNN", sequence1);

            String sequence1withoutNbases = service.getSequenceSliceString(
                    1, sequence.totalBasesWithoutNBases, SequenceRangeOption.WITHOUT_EDGE_N_BASES);
            assertEquals("ACACGTTT", sequence1withoutNbases);
        }
    }

    @Test
    void gettingSequenceViaReaderGivesCorrectResult() throws IOException, SequenceFileException {
        File sequenceFile = TestResources.file("sequence", "example.txt");
        try (SequenceReader service = new SequenceReader(sequenceFile)) {

            SequenceEntry sequence = service.getSequenceInfo();

            // stream whole sequence with the reader
            String streamedSequence;
            try (java.io.Reader r =
                    service.getSequenceSliceReader(1, sequence.totalBases, SequenceRangeOption.WHOLE_SEQUENCE)) {
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
                    1, sequence.totalBasesWithoutNBases, SequenceRangeOption.WITHOUT_EDGE_N_BASES)) {
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
        }
    }

    @Test
    void gettingStringAsAStringVsStreamProducesSameResultSlices() throws IOException, SequenceFileException {
        File sequenceFile = TestResources.file("sequence", "example.txt");
        try (SequenceReader service = new SequenceReader(sequenceFile)) {

            SequenceEntry sequenceInfo = service.getSequenceInfo();

            for (long end = 2; end <= sequenceInfo.totalBases; end++) {
                // get slice as string
                String sequence = service.getSequenceSliceString(1, end, SequenceRangeOption.WHOLE_SEQUENCE);
                // stream sequence with the reader
                String streamedSequence;
                try (java.io.Reader r = service.getSequenceSliceReader(1, end, SequenceRangeOption.WHOLE_SEQUENCE)) {
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
        }
    }
}
