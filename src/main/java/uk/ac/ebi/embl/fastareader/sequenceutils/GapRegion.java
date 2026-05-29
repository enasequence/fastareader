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

import java.util.Objects;

/** A contiguous run of N/n bases in a sequence, using 1-based inclusive base coordinates. */
public final class GapRegion {
    public long startBase;
    public long endBase;

    /** Jackson needs this for JSON->OBJECT conversion. */
    public GapRegion() {}

    public GapRegion(long startBase, long endBase) {
        if (startBase < 1 || endBase < startBase) {
            throw new IllegalArgumentException("bad gap region: " + startBase + ".." + endBase);
        }
        this.startBase = startBase;
        this.endBase = endBase;
    }

    public long lengthBases() {
        return endBase - startBase + 1;
    }

    public boolean overlaps(long fromBase, long toBase) {
        if (fromBase < 1 || toBase < fromBase) {
            throw new IllegalArgumentException("bad base range: " + fromBase + ".." + toBase);
        }
        return startBase <= toBase && endBase >= fromBase;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GapRegion)) return false;
        GapRegion that = (GapRegion) o;
        return startBase == that.startBase && endBase == that.endBase;
    }

    @Override
    public int hashCode() {
        return Objects.hash(startBase, endBase);
    }
}
