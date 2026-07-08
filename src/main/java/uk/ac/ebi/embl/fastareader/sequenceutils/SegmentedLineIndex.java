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

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Piecewise-regular line index: a record is stored as a list of {@link LineSegment}s, each a run of
 * lines with constant width and stride. Memory scales with the number of irregularities (segment
 * boundaries), not the number of lines — a uniformly wrapped record is a single segment, while a
 * record with a handful of blank lines or width changes uses a handful of segments. Lookup binary
 * searches the (few) segments, then computes the byte offset arithmetically within one.
 */
@JsonIgnoreProperties(ignoreUnknown = true)
public final class SegmentedLineIndex implements LineIndex {

    public List<LineSegment> segments;

    /** Jackson needs this for JSON->OBJECT conversion. */
    public SegmentedLineIndex() {
        this.segments = new ArrayList<>();
    }

    public SegmentedLineIndex(List<LineSegment> segments) {
        this.segments = new ArrayList<>(segments);
    }

    @Override
    @JsonIgnore
    public boolean isEmpty() {
        return segments.isEmpty();
    }

    @Override
    public long totalBases() {
        if (segments.isEmpty()) return 0;
        return segments.get(segments.size() - 1).baseEnd;
    }

    @Override
    public LineEntry lineContainingBase(long base) {
        return segments.get(findSegment(base)).lineContainingBase(base);
    }

    @Override
    public List<LineEntry> linesView() {
        List<LineEntry> lines = new ArrayList<>();
        for (LineSegment seg : segments) {
            for (long base = seg.firstBase; base <= seg.baseEnd; base += seg.lineWidthBases) {
                lines.add(seg.lineContainingBase(base));
            }
        }
        return List.copyOf(lines);
    }

    @Override
    public LineIndex copy() {
        return new SegmentedLineIndex(segments); // LineSegment is treated as immutable
    }

    private int findSegment(long base) {
        int lo = 0, hi = segments.size() - 1, ans = hi;
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            LineSegment seg = segments.get(mid);
            if (base < seg.firstBase) hi = mid - 1;
            else if (base > seg.baseEnd) lo = mid + 1;
            else return mid;
            ans = lo;
        }
        return Math.max(0, Math.min(ans, segments.size() - 1));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SegmentedLineIndex)) return false;
        SegmentedLineIndex that = (SegmentedLineIndex) o;
        return Objects.equals(segments, that.segments);
    }

    @Override
    public int hashCode() {
        return Objects.hash(segments);
    }
}
