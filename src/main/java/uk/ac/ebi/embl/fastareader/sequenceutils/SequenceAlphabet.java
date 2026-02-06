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
    private static final byte N_UPPER = (byte) 'N'; // IUPAC defined any character
    private final boolean[] allowed = new boolean[128];
    private final boolean[] specialChars = new boolean[128];

    public SequenceAlphabet(String chars, String specialChars) {
        for (char c : chars.toCharArray()) if (c < 128) allowed[c] = true;
        allowed['>'] = false;

        for (char c : specialChars.toCharArray()) if (c < 128) this.specialChars[c] = true;
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

    /** Fast ASCII check for 'N' or 'n' without decoding. */
    public boolean isNBase(byte b) {
        return (byte) (b & ~0x20) == N_UPPER;
    }

    /** The canonical N base character (uppercase). */
    public char nBase() {
        return (char) N_UPPER; // returns 'N'
    }

    public static SequenceAlphabet defaultNucleotideAlphabet() {
        return new SequenceAlphabet("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "\n\r");
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
}
