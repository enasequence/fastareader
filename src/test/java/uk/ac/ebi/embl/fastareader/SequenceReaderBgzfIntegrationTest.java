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
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import uk.ac.ebi.embl.fastareader.exception.SequenceFileException;
import uk.ac.ebi.embl.fastareader.io.BgzfTestWriter;

/** Round-trip equivalence: a BGZF copy of each plain-sequence fixture must read identically. */
class SequenceReaderBgzfIntegrationTest {

    @TempDir
    Path tmp;

    @Test
    void bgzfReadsIdenticallyToUncompressed_example() throws Exception {
        assertSequenceRoundTrip("example.txt", 16);
    }

    @Test
    void bgzfReadsIdenticallyToUncompressed_carriageReturn() throws Exception {
        assertSequenceRoundTrip("example_with_carriage_return_char.txt", 16);
    }

    @Test
    void bgzfReadsIdenticallyToUncompressed_onlyNs() throws Exception {
        assertSequenceRoundTrip("only_ns_example.txt", 16);
    }

    private void assertSequenceRoundTrip(String fixture, int blockSize) throws Exception {
        File bgzf = BgzfFixtures.bgzfCopyOf("sequence", fixture, tmp, blockSize);

        try (SequenceReader plain = new SequenceReader(TestResources.file("sequence", fixture));
                SequenceReader compressed = new SequenceReader(bgzf)) {

            assertEquals(plain.getSequenceIndex(), compressed.getSequenceIndex(), "index");
            assertEquals(plain.getStats(), compressed.getStats(), "stats");
            assertEquals(plain.getGapRegions(), compressed.getGapRegions(), "gaps");

            long total = plain.getStats().totalBases();
            String expected = plain.getSequenceSliceString(1, total, SequenceRangeOption.WHOLE_SEQUENCE);
            String actual = compressed.getSequenceSliceString(1, total, SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals(expected, actual, "slice string");

            String streamed = streamWhole(compressed, total);
            assertEquals(expected, streamed, "streamed slice");
        }
    }

    private static String streamWhole(SequenceReader reader, long total) throws Exception {
        try (java.io.Reader r = reader.getSequenceSliceReader(1, total, SequenceRangeOption.WHOLE_SEQUENCE)) {
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
     * {@code SequenceReader} with a {@link SequenceFileException} referencing UTF-8.
     */
    @Test
    void bgzfWithNonUtf8DecompressedContentIsRejected() throws Exception {
        // 0x80 is an invalid UTF-8 start byte — will fail Utf8Detector
        byte[] nonUtf8 = {0x41, 0x43, (byte) 0x80, 0x54, 0x0A};
        byte[] bgzfBytes = BgzfTestWriter.toBgzf(nonUtf8, 16);
        Path target = tmp.resolve("non_utf8_seq.bgzf");
        Files.write(target, bgzfBytes);

        SequenceFileException ex =
                assertThrows(SequenceFileException.class, () -> new SequenceReader(target.toFile()).close());
        assertTrue(
                ex.getMessage().contains("UTF-8"), "exception message should mention UTF-8, got: " + ex.getMessage());
    }

    @Test
    void crossBlockSliceStringMatchesStream() throws Exception {
        File bgzf = BgzfFixtures.bgzfCopyOf("sequence", "example.txt", tmp, 4);
        try (SequenceReader compressed = new SequenceReader(bgzf)) {
            long total = compressed.getStats().totalBases();
            for (long end = 1; end <= total; end++) {
                String s = compressed.getSequenceSliceString(1, end, SequenceRangeOption.WHOLE_SEQUENCE);
                String streamed = streamWhole(compressed, end);
                assertEquals(s, streamed, "end " + end);
            }
        }
    }
}
