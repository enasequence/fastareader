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

import org.junit.jupiter.api.Test;

class SequenceAlphabetTest {

    private static byte b(char c) {
        return (byte) c;
    }

    @Test
    void nucleotideDefaultTreatsNAsGap() {
        SequenceAlphabet a = SequenceAlphabet.defaultNucleotideAlphabet();
        assertTrue(a.isGapBase(b('N')));
        assertTrue(a.isGapBase(b('n')));
        assertFalse(a.isGapBase(b('A')));
        assertTrue(a.isAllowedBase(b('N')), "N is still an allowed base for nucleotides");
    }

    @Test
    void proteinDefaultTreatsXAsGapButNAsResidue() {
        SequenceAlphabet a = SequenceAlphabet.defaultProteinAlphabet();

        // X (upper/lower) is the protein gap / unknown character
        assertTrue(a.isGapBase(b('X')));
        assertTrue(a.isGapBase(b('x')));

        // N is Asparagine here: a real residue, NOT a gap
        assertFalse(a.isGapBase(b('N')), "N (Asn) must not be treated as a gap in protein");
        assertTrue(a.isAllowedBase(b('N')));

        // sample of residues + stop are allowed
        for (char c : "ACDEFGHIKLMNPQRSTVWYXBZUO*".toCharArray()) {
            assertTrue(a.isAllowedBase(b(c)), "expected allowed residue: " + c);
        }
        // a nucleotide-only ambiguity code that is not a protein residue is rejected
        assertFalse(a.isAllowedBase(b('J')));
    }

    @Test
    void backwardCompatibleTwoArgConstructorDefaultsGapToNn() {
        SequenceAlphabet a = new SequenceAlphabet("ACGTNacgtn", "\n\r");
        assertTrue(a.isGapBase(b('N')));
        assertTrue(a.isGapBase(b('n')));
        assertFalse(a.isGapBase(b('A')));
    }

    @Test
    void customGapCharsAreHonoured() {
        // aligned protein: both X (unknown) and '-' (alignment gap) count as gap
        SequenceAlphabet a = new SequenceAlphabet("ACDEFGHIKLMNPQRSTVWYX-", "\n\r", "X-");
        assertTrue(a.isGapBase(b('X')));
        assertTrue(a.isGapBase(b('-')));
        assertFalse(a.isGapBase(b('N')));
        assertFalse(a.isGapBase(b('M')));
    }

    @Test
    void emptyGapCharsDisablesGapTracking() {
        SequenceAlphabet a = new SequenceAlphabet("ACGTN", "\n\r", "");
        assertFalse(a.isGapBase(b('N')));
        assertFalse(a.isGapBase(b('X')));
    }

    @Test
    void exportSettingsRoundTripsGapChars() {
        SequenceAlphabet protein = SequenceAlphabet.defaultProteinAlphabet();
        SequenceAlphabetSettings exported = protein.exportAlphabetSettings();
        assertEquals("Xx", exported.gapChars());

        // reconstructing from exported settings preserves gap semantics
        SequenceAlphabet reloaded = new SequenceAlphabet(exported);
        assertTrue(reloaded.isGapBase(b('X')));
        assertFalse(reloaded.isGapBase(b('N')));
    }

    @Test
    void settingsNullGapCharsDefaultsToNn() {
        // simulates settings deserialized from JSON produced before gapChars existed
        SequenceAlphabetSettings legacy = new SequenceAlphabetSettings("ACGTNacgtn", "\n\r", null);
        assertEquals(SequenceAlphabetSettings.DEFAULT_GAP_CHARS, legacy.gapChars());

        SequenceAlphabet a = new SequenceAlphabet(legacy);
        assertTrue(a.isGapBase(b('N')));
    }
}
