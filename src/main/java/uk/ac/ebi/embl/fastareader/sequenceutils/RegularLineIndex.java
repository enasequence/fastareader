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
 * Constant-memory line index for a record whose lines all share one base width and one byte stride,
 * except possibly a shorter final line. Any line's byte offset is computed arithmetically, so the
 * whole record costs four longs regardless of how many gigabases it contains.
 */
@JsonIgnoreProperties(ignoreUnknown = true)
public final class RegularLineIndex implements LineIndex {

    public long firstByteStart; // absolute byte offset of the first base
    public long lineWidthBases; // bases per full line (W)
    public long byteStride; // bytes between consecutive line starts (W + separator width)
    public long totalBaseCount; // total bases in the record

    /** Jackson needs this for JSON->OBJECT conversion. */
    public RegularLineIndex() {}

    public RegularLineIndex(long firstByteStart, long lineWidthBases, long byteStride, long totalBaseCount) {
        if (lineWidthBases <= 0) {
            throw new IllegalArgumentException("lineWidthBases must be positive: " + lineWidthBases);
        }
        this.firstByteStart = firstByteStart;
        this.lineWidthBases = lineWidthBases;
        this.byteStride = byteStride;
        this.totalBaseCount = totalBaseCount;
    }

    @Override
    @JsonIgnore
    public boolean isEmpty() {
        return totalBaseCount <= 0;
    }

    @Override
    public long totalBases() {
        return totalBaseCount;
    }

    @Override
    public LineEntry lineContainingBase(long base) {
        long k = (base - 1) / lineWidthBases; // 0-based line index
        long baseStart = 1 + k * lineWidthBases;
        long baseEnd = Math.min(baseStart + lineWidthBases - 1, totalBaseCount);
        long byteStart = firstByteStart + k * byteStride;
        long byteEndExclusive = byteStart + (baseEnd - baseStart + 1);
        return new LineEntry(baseStart, baseEnd, byteStart, byteEndExclusive);
    }

    @Override
    public List<LineEntry> linesView() {
        List<LineEntry> lines = new ArrayList<>();
        for (long baseStart = 1; baseStart <= totalBaseCount; baseStart += lineWidthBases) {
            lines.add(lineContainingBase(baseStart));
        }
        return List.copyOf(lines);
    }

    @Override
    public LineIndex copy() {
        return new RegularLineIndex(firstByteStart, lineWidthBases, byteStride, totalBaseCount);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof RegularLineIndex)) return false;
        RegularLineIndex that = (RegularLineIndex) o;
        return firstByteStart == that.firstByteStart
                && lineWidthBases == that.lineWidthBases
                && byteStride == that.byteStride
                && totalBaseCount == that.totalBaseCount;
    }

    @Override
    public int hashCode() {
        return Objects.hash(firstByteStart, lineWidthBases, byteStride, totalBaseCount);
    }
}
