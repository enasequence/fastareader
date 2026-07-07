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
import java.util.Collections;
import java.util.List;

public final class SequenceAlphabet {
    private final boolean[] allowed = new boolean[128];
    private final boolean[] specialChars = new boolean[128];
    private final boolean[] gapChars = new boolean[128];
    private final SequenceAlphabetSettings settings;

    /** Gap characters default to N/n (nucleotide behaviour) when not otherwise specified. */
    public SequenceAlphabet(String nucleotideString, String specialCharsString) {
        this.settings = new SequenceAlphabetSettings(nucleotideString, specialCharsString);
        setupAsciiArrays();
    }

    /**
     * @param baseString the allowed sequence characters (bases / residues)
     * @param specialCharsString characters tolerated but ignored within a sequence (e.g. {@code "\n\r"})
     * @param gapCharsString characters treated as gap / unknown for gap-region detection and edge trimming
     *     (e.g. {@code "N"} for nucleotides, {@code "X"} for protein); may be empty to disable gap tracking
     */
    public SequenceAlphabet(String baseString, String specialCharsString, String gapCharsString) {
        this.settings = new SequenceAlphabetSettings(baseString, specialCharsString, gapCharsString);
        setupAsciiArrays();
    }

    public SequenceAlphabet(SequenceAlphabetSettings settings) {
        this.settings = settings;
        setupAsciiArrays();
    }

    /** Fast ASCII check for is it an allowed char. */
    public boolean isAllowedBase(byte b) {
        int i = b & 0xFF;
        return i < 128 && allowed[i];
    }

    public boolean isNonSequenceAllowedChar(byte b) {
        int i = b & 0xFF;
        return i < 128 && specialChars[i];
    }

    /**
     * Fast ASCII check for whether {@code b} is a gap / unknown character for this alphabet (e.g. {@code 'N'}
     * for nucleotides, {@code 'X'} for protein). Drives gap-region detection and edge trimming.
     */
    public boolean isGapBase(byte b) {
        int i = b & 0xFF;
        return i < 128 && gapChars[i];
    }

    /**
     * @deprecated superseded by {@link #isGapBase(byte)}, which reflects the alphabet's configured gap
     *     characters rather than a hardcoded {@code 'N'}. Retained for source compatibility; for the
     *     default nucleotide alphabet it is behaviourally identical.
     */
    @Deprecated
    public boolean isNBase(byte b) {
        return isGapBase(b);
    }

    public static SequenceAlphabet defaultNucleotideAlphabet() {
        return new SequenceAlphabet(new SequenceAlphabetSettings("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "\n\r", "Nn"));
    }

    /**
     * Amino-acid alphabet for protein / translated sequences: the 20 standard residues plus ambiguity codes
     * {@code X B Z}, non-standard residues {@code U O}, and the stop character {@code *} (upper and lower
     * case). The gap / unknown character is {@code X} — note {@code N} is a valid residue (Asparagine) here,
     * not a gap.
     */
    public static SequenceAlphabet defaultProteinAlphabet() {
        return new SequenceAlphabet(
                new SequenceAlphabetSettings("ACDEFGHIKLMNPQRSTVWYXBZUO*acdefghiklmnpqrstvwyxbzuo", "\n\r", "Xx"));
    }

    /** Returns allowed bases as uppercase, de-duplicated (e.g., 'A' not both 'A' and 'a'). */
    public List<Character> getAllowedBaseCharList() {
        List<Character> out = new ArrayList<>();

        for (int i = 0; i < 128; i++) {
            if (!allowed[i]) continue;
            char c = (char) i;
            out.add(c);
        }

        return Collections.unmodifiableList(out);
    }

    public String describeAllowed() {
        StringBuilder sb = new StringBuilder();
        sb.append("[");

        boolean first = true;
        for (int i = 0; i < allowed.length; i++) {
            if (allowed[i]) {
                char c = (char) i;
                String display;
                display = Character.toString(c);
                if (!first) sb.append(", ");
                sb.append(display);
                first = false;
            }
        }

        sb.append("]");
        return sb.toString();
    }

    public SequenceAlphabetSettings exportAlphabetSettings() {
        return new SequenceAlphabetSettings(settings.chars(), settings.specialChars(), settings.gapChars());
    }

    private void setupAsciiArrays() {
        for (char c : this.settings.chars().toCharArray()) if (c < 128) allowed[c] = true;
        allowed['>'] = false;

        for (char c : this.settings.specialChars().toCharArray()) if (c < 128) this.specialChars[c] = true;

        String gaps = this.settings.gapChars();
        if (gaps != null) for (char c : gaps.toCharArray()) if (c < 128) this.gapChars[c] = true;
    }
}
