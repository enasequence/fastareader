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
package uk.ac.ebi.embl.fastareader.sequenceutils;

import static org.junit.jupiter.api.Assertions.*;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

public class LineIndexTest {

    /** Explicit layout of a regular run: fixed base width, fixed byte stride. */
    private static List<LineEntry> regularLines(long lineCount, long width, long stride, long firstByteStart) {
        List<LineEntry> lines = new ArrayList<>();
        for (long k = 0; k < lineCount; k++) {
            long baseStart = 1 + k * width;
            long byteStart = firstByteStart + k * stride;
            lines.add(new LineEntry(baseStart, baseStart + width - 1, byteStart, byteStart + width));
        }
        return lines;
    }

    private static LineEntry lineContaining(List<LineEntry> lines, long base) {
        for (LineEntry line : lines) {
            if (base >= line.baseStart && base <= line.baseEnd) return line;
        }
        throw new IllegalArgumentException("no line contains base " + base);
    }

    private static void assertAgreesWithLines(LineIndex index, List<LineEntry> lines) {
        long total = lines.get(lines.size() - 1).baseEnd;
        assertEquals(total, index.totalBases());
        for (long base = 1; base <= total; base++) {
            LineEntry expected = lineContaining(lines, base);
            LineEntry actual = index.lineContainingBase(base);
            assertEquals(expected.baseStart, actual.baseStart, "baseStart @" + base);
            assertEquals(expected.baseEnd, actual.baseEnd, "baseEnd @" + base);
            assertEquals(expected.byteStart, actual.byteStart, "byteStart @" + base);
        }
        assertEquals(lines, index.linesView());
    }

    @Test
    void uniformRecordCollapsesToRegular() {
        List<LineEntry> lines = regularLines(7, 4, 5, 100);
        LineIndex index = LineSegmenter.fromLines(lines);
        assertInstanceOf(RegularLineIndex.class, index, "one segment starting at base 1 -> regular");
        assertAgreesWithLines(index, lines);
    }

    @Test
    void shortFinalLineStaysRegular() {
        List<LineEntry> lines = new ArrayList<>(regularLines(3, 4, 5, 0));
        long tailStart = 3 * 5;
        lines.add(new LineEntry(13, 14, tailStart, tailStart + 2)); // 2-base tail
        LineIndex index = LineSegmenter.fromLines(lines);
        assertInstanceOf(RegularLineIndex.class, index, "a shorter final line is not an irregularity");
        assertAgreesWithLines(index, lines);
    }

    @Test
    void strideChangeProducesTwoSegments() {
        // Third line starts one byte late (as if a blank line were inserted): stride 5 -> 6.
        List<LineEntry> lines = List.of(
                new LineEntry(1, 4, 100, 104), new LineEntry(5, 8, 105, 109), new LineEntry(9, 12, 111, 115));
        LineIndex index = LineSegmenter.fromLines(lines);
        assertInstanceOf(SegmentedLineIndex.class, index);
        assertEquals(2, ((SegmentedLineIndex) index).segments.size());
        assertAgreesWithLines(index, lines);
    }

    @Test
    void widthChangeMidRecordProducesSegments() {
        // Widths 4,4,2,4,4 -> the short middle line splits the run.
        List<LineEntry> lines = List.of(
                new LineEntry(1, 4, 0, 4),
                new LineEntry(5, 8, 5, 9),
                new LineEntry(9, 10, 10, 12), // short, non-final
                new LineEntry(11, 14, 13, 17),
                new LineEntry(15, 18, 18, 22));
        LineIndex index = LineSegmenter.fromLines(lines);
        assertInstanceOf(SegmentedLineIndex.class, index);
        assertEquals(2, ((SegmentedLineIndex) index).segments.size());
        assertAgreesWithLines(index, lines);
    }

    @Test
    void singleLineRecordIsRegular() {
        List<LineEntry> lines = List.of(new LineEntry(1, 6, 50, 56));
        LineIndex index = LineSegmenter.fromLines(lines);
        assertInstanceOf(RegularLineIndex.class, index);
        assertAgreesWithLines(index, lines);
    }

    @Test
    void emptyRecordIsEmptySegmented() {
        LineIndex index = LineSegmenter.fromLines(List.of());
        assertInstanceOf(SegmentedLineIndex.class, index);
        assertTrue(index.isEmpty());
        assertEquals(0, index.totalBases());
        assertTrue(index.linesView().isEmpty());
    }

    @Test
    void bothRepresentationsSurviveJsonRoundTrip() throws Exception {
        ObjectMapper mapper = new ObjectMapper();

        LineIndex regular = new RegularLineIndex(100, 60, 61, 1_000_000_000L);
        LineIndex segmented = LineSegmenter.fromLines(List.of(
                new LineEntry(1, 4, 100, 104), new LineEntry(5, 8, 105, 109), new LineEntry(9, 12, 111, 115)));
        assertInstanceOf(SegmentedLineIndex.class, segmented);

        for (LineIndex original : List.of(regular, segmented)) {
            String json = mapper.writeValueAsString(original);
            LineIndex back = mapper.readValue(json, LineIndex.class);
            assertEquals(original.getClass(), back.getClass(), "concrete type preserved");
            assertEquals(original, back);
            assertEquals(original.totalBases(), back.totalBases());
        }
    }
}
