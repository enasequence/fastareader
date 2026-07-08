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

import java.util.ArrayList;
import java.util.List;

/**
 * Streaming builder that groups sequence lines into {@link LineSegment}s greedily, retaining only
 * the currently open segment plus the closed ones. Feeding lines one at a time keeps memory bounded
 * by the number of irregularities rather than the number of lines, so a uniformly wrapped
 * gigabase record collapses to a single segment.
 */
public final class LineSegmenter {

    private final List<LineSegment> segments = new ArrayList<>();

    private boolean hasOpen = false;
    private long firstBase;
    private long firstByte;
    private long width; // W of the open segment
    private long stride = -1; // S of the open segment, set once a second line arrives
    private long lastByteStart;
    private long lastWidth;
    private long baseEnd;
    private long lineCount;

    /** Feed one committed line (1-based inclusive bases, absolute bytes) in ascending base order. */
    public void add(long baseStart, long baseEnd, long byteStart, long byteEndExclusive) {
        long lineWidth = baseEnd - baseStart + 1;
        if (!hasOpen) {
            open(baseStart, baseEnd, byteStart, lineWidth);
            return;
        }
        long lineStride = byteStart - lastByteStart;
        boolean strideOk = (lineCount == 1) || (lineStride == stride);
        boolean prevLineFull = (lastWidth == width); // a non-final line in the run must be full width
        boolean widthFits = (lineWidth <= width); // no line may exceed the run's width
        if (strideOk && prevLineFull && widthFits) {
            if (lineCount == 1) stride = lineStride;
            lastByteStart = byteStart;
            lastWidth = lineWidth;
            this.baseEnd = baseEnd;
            lineCount++;
        } else {
            close();
            open(baseStart, baseEnd, byteStart, lineWidth);
        }
    }

    private void open(long baseStart, long baseEnd, long byteStart, long lineWidth) {
        hasOpen = true;
        firstBase = baseStart;
        firstByte = byteStart;
        width = lineWidth;
        stride = -1;
        lastByteStart = byteStart;
        lastWidth = lineWidth;
        this.baseEnd = baseEnd;
        lineCount = 1;
    }

    private void close() {
        long segStride = (lineCount >= 2) ? stride : Math.max(width, 1);
        segments.add(new LineSegment(firstBase, firstByte, width, segStride, baseEnd));
        hasOpen = false;
    }

    /** Finalize and return the most compact representation for the accumulated lines. */
    public LineIndex build() {
        if (hasOpen) close();
        if (segments.size() == 1) {
            LineSegment seg = segments.get(0);
            if (seg.firstBase == 1) {
                return new RegularLineIndex(seg.firstByte, seg.lineWidthBases, seg.byteStride, seg.baseEnd);
            }
        }
        return new SegmentedLineIndex(segments);
    }

    /** Convenience for building an index from an in-memory list (used outside the file scanner). */
    public static LineIndex fromLines(List<LineEntry> lines) {
        LineSegmenter segmenter = new LineSegmenter();
        for (LineEntry line : lines) {
            segmenter.add(line.baseStart, line.baseEnd, line.byteStart, line.byteEndExclusive);
        }
        return segmenter.build();
    }
}
