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
package uk.ac.ebi.embl.fastareader.api.rereading;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import java.nio.file.Path;
import java.util.HashMap;
import uk.ac.ebi.embl.fastareader.SequenceFileFormat;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabetSettings;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceIndex;

@JsonIgnoreProperties(ignoreUnknown = true)
public class SequenceInfoDTO {

    public Path filePath;
    public SequenceFileFormat sequenceFileFormat;
    public SequenceAlphabetSettings sequenceAlphabetSettings;
    public HashMap<Long, SequenceIndex> sequenceIndexesMap;
    public HashMap<Long, String> headerLines;
}
