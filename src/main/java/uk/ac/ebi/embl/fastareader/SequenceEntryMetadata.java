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

import lombok.Getter;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceIndex;

class SequenceEntryMetadata {

    /**
     * Assumed file format based on the input
     */
    @Getter
    FileFormat fileFormat;
    /**
     * information needed for accessing the file, if doesn't apply due to lack of header it's null
     */
    @Getter
    private String headerLine;
    /**
     * position of '>' in the file, if doesn't apply due to file not being in fasta file format it's going to be -1L
     */
    @Getter
    long fastaStartByte;
    /**
     * a smart index for querying ranges in the file
     */
    @Getter
    SequenceIndex sequenceIndex;

    /**
     * Constructor for sequence metadata about a plain sequence
     */
    SequenceEntryMetadata(SequenceIndex sequenceIndex) {
        this.fileFormat = FileFormat.PLAIN_SINGLE_SEQUENCE;
        this.sequenceIndex = sequenceIndex;
        this.fastaStartByte = -1;
        this.headerLine = null;
    }

    /**
     * Constructor for sequence metadata about a FASTA header+sequence
     */
    public SequenceEntryMetadata(String headerLine, long fastaStartByte, SequenceIndex sequenceIndex) {
        this.fileFormat = FileFormat.FASTA;
        this.headerLine = headerLine;
        this.fastaStartByte = fastaStartByte;
        this.sequenceIndex = sequenceIndex;
    }
}
