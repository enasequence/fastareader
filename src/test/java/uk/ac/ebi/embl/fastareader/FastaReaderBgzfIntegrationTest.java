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
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

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
}
