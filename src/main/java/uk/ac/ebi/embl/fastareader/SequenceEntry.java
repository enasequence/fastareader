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
package uk.ac.ebi.embl.fastareader;

import java.util.Map;
import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class SequenceEntry {
    public long totalBases;
    public long totalBasesWithoutNBases;
    public long leadingNsCount;
    public long trailingNsCount;
    /**
     * Count of all appearances of permitted sequence characters in this entry.
     * <p>
     * Letter keys are canonicalized to uppercase (e.g., 'a' and 'A' are counted under 'A').
     * Non-letter allowed characters (e.g., '-', '.', '*') appear as-is.
     * Allowed characters not present have a count of 0.
     */
    public Map<Character, Long> baseCount;
}
