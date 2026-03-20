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

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.nio.file.Path;
import java.util.HashMap;
import uk.ac.ebi.embl.fastareader.SequenceFileFormat;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceIndex;

@JsonIgnoreProperties(ignoreUnknown = true)
public class SequenceFormatReaderDTO {

    private final Path filePath;
    private final SequenceFileFormat sequenceFileFormat;
    private final SequenceAlphabet sequenceAlphabet;
    private final HashMap<Long, SequenceIndex>
            sequenceIndexesMap; // fastaReaderId -> sequenceIndex, should have one value  in case of SequenceFileFormat - Plain sequence
    private final HashMap<Long, String>
            headerLines; // fastaReaderId -> headerLine, should be null in case of SequenceFileFormat = Plain sequence

    @JsonCreator
    public SequenceFormatReaderDTO(
            @JsonProperty("filePath") Path filePath,
            @JsonProperty("sequenceFileFormat") SequenceFileFormat sequenceFileFormat,
            @JsonProperty("alphabet") SequenceAlphabet sequenceAlphabet,
            @JsonProperty("sequenceIndexesMap") HashMap<Long, SequenceIndex> sequenceIndexesMap,
            @JsonProperty("headerLines") HashMap<Long, String> headerLines) {
        this.filePath = filePath;
        this.sequenceFileFormat = sequenceFileFormat;
        this.sequenceAlphabet = sequenceAlphabet;
        this.sequenceIndexesMap = sequenceIndexesMap;
        this.headerLines = headerLines;
    }

    @JsonProperty("filePath")
    public Path getFilePath() {
        return filePath;
    }

    @JsonProperty("sequenceFileFormat")
    public SequenceFileFormat getSequenceFileFormat() {
        return sequenceFileFormat;
    }

    @JsonProperty("alphabet")
    public SequenceAlphabet getSequenceAlphabet() {
        return sequenceAlphabet;
    }

    @JsonProperty("sequenceIndexesMap")
    public HashMap<Long, SequenceIndex> getSequenceIndexesMap() {
        return sequenceIndexesMap;
    }

    @JsonProperty("headerLines")
    public HashMap<Long, String> getHeaderLines() {
        return headerLines;
    }
}
