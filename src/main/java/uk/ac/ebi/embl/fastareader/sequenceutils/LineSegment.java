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

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import java.util.Objects;

/**
 * A maximal run of consecutive sequence lines that share one base width and one byte stride, with
 * only the final line of the run allowed to be shorter. Any base within the run is located
 * arithmetically, so the whole run costs five longs regardless of how many lines it spans.
 */
@JsonIgnoreProperties(ignoreUnknown = true)
public final class LineSegment {

    public long firstBase; // 1-based base of the first line in this run
    public long firstByte; // absolute byte offset of that first base
    public long lineWidthBases; // bases per full line (W)
    public long byteStride; // bytes between consecutive line starts (W + separator width)
    public long baseEnd; // 1-based base of the last base in this run (inclusive)

    /** Jackson needs this for JSON->OBJECT conversion. */
    public LineSegment() {}

    public LineSegment(long firstBase, long firstByte, long lineWidthBases, long byteStride, long baseEnd) {
        if (lineWidthBases <= 0) {
            throw new IllegalArgumentException("lineWidthBases must be positive: " + lineWidthBases);
        }
        if (baseEnd < firstBase) {
            throw new IllegalArgumentException("bad segment: " + firstBase + ".." + baseEnd);
        }
        this.firstBase = firstBase;
        this.firstByte = firstByte;
        this.lineWidthBases = lineWidthBases;
        this.byteStride = byteStride;
        this.baseEnd = baseEnd;
    }

    public boolean contains(long base) {
        return base >= firstBase && base <= baseEnd;
    }

    public long totalBases() {
        return baseEnd - firstBase + 1;
    }

    public LineEntry lineContainingBase(long base) {
        long k = (base - firstBase) / lineWidthBases; // 0-based line index within this run
        long lineBaseStart = firstBase + k * lineWidthBases;
        long lineBaseEnd = Math.min(lineBaseStart + lineWidthBases - 1, baseEnd);
        long lineByteStart = firstByte + k * byteStride;
        long lineByteEndExclusive = lineByteStart + (lineBaseEnd - lineBaseStart + 1);
        return new LineEntry(lineBaseStart, lineBaseEnd, lineByteStart, lineByteEndExclusive);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof LineSegment)) return false;
        LineSegment that = (LineSegment) o;
        return firstBase == that.firstBase
                && firstByte == that.firstByte
                && lineWidthBases == that.lineWidthBases
                && byteStride == that.byteStride
                && baseEnd == that.baseEnd;
    }

    @Override
    public int hashCode() {
        return Objects.hash(firstBase, firstByte, lineWidthBases, byteStride, baseEnd);
    }
}
