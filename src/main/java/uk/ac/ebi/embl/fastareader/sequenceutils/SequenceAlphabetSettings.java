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

import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * Immutable description of a sequence alphabet.
 *
 * @param chars the allowed sequence characters (bases / residues), e.g. {@code "ACGT..."}
 * @param specialChars characters tolerated but ignored within a sequence (typically end-of-line: {@code "\n\r"})
 * @param gapChars characters treated as "gap / unknown" for gap-region detection and edge trimming
 *     (e.g. {@code "N"} for nucleotides, {@code "X"} for protein). May be empty to disable gap tracking.
 *     A {@code null} value defaults to {@code "Nn"} for backward compatibility with settings serialized
 *     before this field existed.
 */
public record SequenceAlphabetSettings(
        @JsonProperty("chars") String chars,
        @JsonProperty("specialChars") String specialChars,
        @JsonProperty("gapChars") String gapChars) {

    /** Default gap characters when none are supplied (matches historical N/n behaviour). */
    public static final String DEFAULT_GAP_CHARS = "Nn";

    public SequenceAlphabetSettings {
        if (gapChars == null) gapChars = DEFAULT_GAP_CHARS;
    }

    /** Backward-compatible constructor; gap characters default to {@value #DEFAULT_GAP_CHARS}. */
    public SequenceAlphabetSettings(String chars, String specialChars) {
        this(chars, specialChars, DEFAULT_GAP_CHARS);
    }
}
