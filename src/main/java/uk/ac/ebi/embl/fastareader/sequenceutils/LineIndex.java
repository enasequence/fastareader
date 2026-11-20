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

import com.fasterxml.jackson.annotation.JsonSubTypes;
import com.fasterxml.jackson.annotation.JsonTypeInfo;
import java.util.List;

/**
 * Maps 1-based sequence base coordinates to absolute byte offsets within a single sequence record.
 *
 * <p>Two representations exist: {@link RegularLineIndex} stores a uniformly wrapped record in O(1)
 * memory (four longs, regardless of record length), while {@link SegmentedLineIndex} stores one
 * {@link LineSegment} per regular run for records whose line widths or separators change, costing
 * memory proportional to the number of irregularities rather than the number of lines.
 */
@JsonTypeInfo(use = JsonTypeInfo.Id.NAME, include = JsonTypeInfo.As.PROPERTY, property = "kind")
@JsonSubTypes({
    @JsonSubTypes.Type(value = RegularLineIndex.class, name = "regular"),
    @JsonSubTypes.Type(value = SegmentedLineIndex.class, name = "segmented")
})
public interface LineIndex {

    boolean isEmpty();

    long totalBases();

    /** Return the (possibly synthesized) line entry containing the given 1-based base. */
    LineEntry lineContainingBase(long base);

    /**
     * Materialize all lines as explicit entries. O(lines) regardless of representation — intended
     * for inspection/tests, not the hot read path.
     */
    List<LineEntry> linesView();

    LineIndex copy();
}
