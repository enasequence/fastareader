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

import java.util.*;
import uk.ac.ebi.embl.fastareader.SequenceRangeOption;

public final class SequenceIndex {

    public long firstBaseByte; // -1 if empty
    public long lastBaseByte; // -1 if empty
    public List<LineEntry> lines;
    public long nextHeaderByte; // byte offset of next '>' at line start, or fileSize (EOF)
    // Bases related counts
    public long startNBasesCount;
    public long endNBasesCount;
    public List<GapRegion> gapRegions;
    public Map<Character, Long> caseInsensitiveBaseCount;

    /** Jackson needs this for JSON->OBJECT conversion **/
    public SequenceIndex() {}

    /*** A constructor that is only used by the SequenceIndexBuilder before it counts the actual characters and builds the full index. ***/
    SequenceIndex(
            long firstBaseByte,
            long startNBasesCount,
            long lastBaseByte,
            long endNBasesCount,
            List<LineEntry> lines,
            long nextHeader,
            long[] baseCounts,
            List<Character> baseChars) {
        this(
                firstBaseByte,
                startNBasesCount,
                lastBaseByte,
                endNBasesCount,
                lines,
                nextHeader,
                baseCounts,
                baseChars,
                Collections.emptyList());
    }

    /*** A constructor used by the SequenceIndexBuilder once it has counted bases and gap regions. ***/
    SequenceIndex(
            long firstBaseByte,
            long startNBasesCount,
            long lastBaseByte,
            long endNBasesCount,
            List<LineEntry> lines,
            long nextHeader,
            long[] baseCounts,
            List<Character> baseChars,
            List<GapRegion> gapRegions) {
        this.firstBaseByte = firstBaseByte;
        this.startNBasesCount = startNBasesCount;
        this.lastBaseByte = lastBaseByte;
        this.endNBasesCount = endNBasesCount;
        this.lines = new ArrayList<>(lines);
        this.nextHeaderByte = nextHeader;
        this.gapRegions = new ArrayList<>(gapRegions);
        this.caseInsensitiveBaseCount = new HashMap<>();
        mapAllowedCharacters(baseCounts, baseChars);
    }

    /** Copy constructor **/
    public SequenceIndex(SequenceIndex other) {
        this.firstBaseByte = other.firstBaseByte;
        this.lastBaseByte = other.lastBaseByte;
        this.nextHeaderByte = other.nextHeaderByte;
        this.startNBasesCount = other.startNBasesCount;
        this.endNBasesCount = other.endNBasesCount;
        this.caseInsensitiveBaseCount =
                other.caseInsensitiveBaseCount != null ? new HashMap<>(other.caseInsensitiveBaseCount) : null;
        this.lines = other.lines != null ? new ArrayList<>(other.lines) : null;
        this.gapRegions = other.gapRegions != null ? new ArrayList<>(other.gapRegions) : new ArrayList<>();
    }

    private void mapAllowedCharacters(long[] baseCounts, List<Character> allowedCanonicalChars) {
        for (char c : allowedCanonicalChars) {
            // for alphabet, sum up lowercase once
            if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')) {
                char lowercase = Character.toLowerCase(c);
                char uppercase = Character.toUpperCase(c);
                if (!caseInsensitiveBaseCount.containsKey(uppercase)) {
                    caseInsensitiveBaseCount.put(c, baseCounts[lowercase] + baseCounts[uppercase]);
                }
            }
            // everything else is counted once
            else {
                caseInsensitiveBaseCount.put(c, baseCounts[c]);
            }
        }
    }

    public List<LineEntry> linesView() {
        return Collections.unmodifiableList(lines);
    }

    public long totalBases() {
        if (lines.isEmpty()) return 0;
        return lines.get(lines.size() - 1).baseEnd;
    }

    public long totalBasesExcludingEdgeNBases() {
        long bases = totalBases() - endNBasesCount - startNBasesCount;
        return Math.max(0, bases);
    }

    public List<GapRegion> gapRegionsView() {
        if (gapRegions == null) return Collections.emptyList();
        return Collections.unmodifiableList(gapRegions);
    }

    public List<GapRegion> gapRegionsView(SequenceRangeOption option) {
        switch (option) {
            case WHOLE_SEQUENCE:
                return gapRegionsView();
            case WITHOUT_EDGE_N_BASES:
                long trimmedTotal = totalBasesExcludingEdgeNBases();
                if (trimmedTotal == 0) return Collections.emptyList();
                return gapRegionsView(1, trimmedTotal, option);
            default:
                throw new IllegalStateException("Unknown option " + option);
        }
    }

    public List<GapRegion> gapRegionsView(long fromBase, long toBase) {
        return gapRegionsView(fromBase, toBase, SequenceRangeOption.WHOLE_SEQUENCE);
    }

    public List<GapRegion> gapRegionsView(long fromBase, long toBase, SequenceRangeOption option) {
        switch (option) {
            case WHOLE_SEQUENCE:
                return gapRegionsViewIncludingEdgeNBases(fromBase, toBase);
            case WITHOUT_EDGE_N_BASES:
                long trimmedTotal = totalBasesExcludingEdgeNBases();
                if (fromBase < 1 || toBase < fromBase || toBase > trimmedTotal) {
                    throw new IllegalArgumentException("bad base range: " + fromBase + ".." + toBase);
                }
                long actualFromBase = startNBasesCount + fromBase;
                long actualToBase = startNBasesCount + toBase;
                return shiftGapRegions(
                        gapRegionsViewIncludingEdgeNBases(actualFromBase, actualToBase), -startNBasesCount);
            default:
                throw new IllegalStateException("Unknown option " + option);
        }
    }

    private List<GapRegion> gapRegionsViewIncludingEdgeNBases(long fromBase, long toBase) {
        long total = totalBases();
        if (fromBase < 1 || toBase < fromBase || toBase > total) {
            throw new IllegalArgumentException("bad base range: " + fromBase + ".." + toBase);
        }
        if (gapRegions == null || gapRegions.isEmpty()) return Collections.emptyList();
        List<GapRegion> matches = new ArrayList<>();
        for (GapRegion gapRegion : gapRegions) {
            if (gapRegion.overlaps(fromBase, toBase)) {
                matches.add(gapRegion);
            }
        }
        return Collections.unmodifiableList(matches);
    }

    private static List<GapRegion> shiftGapRegions(List<GapRegion> gapRegions, long offset) {
        if (gapRegions.isEmpty()) return Collections.emptyList();
        List<GapRegion> shifted = new ArrayList<>(gapRegions.size());
        for (GapRegion gapRegion : gapRegions) {
            shifted.add(new GapRegion(gapRegion.startBase + offset, gapRegion.endBase + offset));
        }
        return Collections.unmodifiableList(shifted);
    }

    public ByteSpan byteSpanForBaseRangeIncludingEdgeNBases(long fromBase, long toBase) {
        long total = totalBases();
        if (fromBase < 1 || toBase < fromBase || toBase > total) {
            throw new IllegalArgumentException("bad base range: " + fromBase + ".." + toBase);
        }
        int i = findLineByBase(fromBase);
        int j = findLineByBase(toBase);

        LineEntry from = lines.get(i);
        long offStart = fromBase - from.baseStart;

        LineEntry to = lines.get(j);
        long offEndIncl = toBase - to.baseStart;

        long byteStart = from.byteStart + offStart;
        long byteEndEx = to.byteStart + offEndIncl + 1; // half-open

        return new ByteSpan(byteStart, byteEndEx);
    }

    public ByteSpan byteSpanForBaseRange(long fromBase, long toBase) {
        long trimmedTotal = totalBasesExcludingEdgeNBases();
        if (fromBase < 1 || toBase < fromBase || toBase > trimmedTotal) {
            throw new IllegalArgumentException("bad base range: " + fromBase + ".." + toBase);
        }
        long actualFromBase = startNBasesCount + fromBase;
        long actualToBase = startNBasesCount + toBase;
        return byteSpanForBaseRangeIncludingEdgeNBases(actualFromBase, actualToBase);
    }

    private int findLineByBase(long base) {
        int lo = 0, hi = lines.size() - 1, ans = hi;
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            LineEntry L = lines.get(mid);
            if (base < L.baseStart) hi = mid - 1;
            else if (base > L.baseEnd) lo = mid + 1;
            else return mid;
            ans = lo;
        }
        return Math.max(0, Math.min(ans, lines.size() - 1));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SequenceIndex)) return false;
        SequenceIndex that = (SequenceIndex) o;
        return firstBaseByte == that.firstBaseByte
                && lastBaseByte == that.lastBaseByte
                && nextHeaderByte == that.nextHeaderByte
                && startNBasesCount == that.startNBasesCount
                && endNBasesCount == that.endNBasesCount
                && Objects.equals(lines, that.lines)
                && Objects.equals(gapRegionsView(), that.gapRegionsView())
                && Objects.equals(caseInsensitiveBaseCount, that.caseInsensitiveBaseCount);
    }

    @Override
    public int hashCode() {
        return Objects.hash(
                firstBaseByte,
                lastBaseByte,
                nextHeaderByte,
                startNBasesCount,
                endNBasesCount,
                lines,
                gapRegionsView(),
                caseInsensitiveBaseCount);
    }
}
