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
package uk.ac.ebi.embl.fastareader.api;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import uk.ac.ebi.embl.fastareader.SequenceFileFormat;
import uk.ac.ebi.embl.fastareader.SequenceRangeOption;
import uk.ac.ebi.embl.fastareader.SequenceReader;
import uk.ac.ebi.embl.fastareader.SequenceStats;
import uk.ac.ebi.embl.fastareader.api.rereading.SequenceInfoDTO;
import uk.ac.ebi.embl.fastareader.exception.SequenceFileException;
import uk.ac.ebi.embl.fastareader.sequenceutils.GapRegion;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceIndex;

public class SequenceReaderWrapper implements AutoCloseable, SequenceFormatReader {

    private final SequenceReader sequenceReader;
    private final long precodedId = 0L;

    public SequenceReaderWrapper(File sequenceFile) throws SequenceFileException, IOException {
        this.sequenceReader = new SequenceReader(sequenceFile);
    }

    public SequenceReaderWrapper(File sequenceFile, SequenceAlphabet alphabet)
            throws SequenceFileException, IOException {
        this.sequenceReader = new SequenceReader(sequenceFile, alphabet);
    }

    public SequenceReaderWrapper(
            File file, SequenceAlphabet sequenceAlphabet, HashMap<Long, SequenceIndex> sequenceIndexesMap)
            throws SequenceFileException, IOException {
        if (sequenceIndexesMap.size() != 1 && !sequenceIndexesMap.containsKey(precodedId)) {
            throw new IllegalArgumentException(
                    "Sequence files must have exactly one entry which should be indexed with " + precodedId + " .");
        }
        this.sequenceReader = new SequenceReader(file, sequenceAlphabet, sequenceIndexesMap.get(precodedId));
    }

    @Override
    public SequenceFileFormat getSequenceFileFormat() {
        return sequenceReader.getSequenceFileFormat();
    }

    @Override
    public File getFile() {
        return sequenceReader.getFile();
    }

    @Override
    public SequenceAlphabet getSequenceAlphabet() {
        return sequenceReader.getSequenceAlphabet();
    }

    @Override
    public List<Long> getOrderedIds() {
        return List.of(precodedId);
    }

    @Override
    public Optional<String> getHeaderline(long id) {
        return Optional.empty();
    }

    @Override
    public SequenceStats getStats(long id) {
        validateId(id);
        return sequenceReader.getStats();
    }

    @Override
    public List<GapRegion> getGapRegions(long id) {
        validateId(id);
        return sequenceReader.getGapRegions();
    }

    @Override
    public List<GapRegion> getGapRegions(long id, long fromBase, long toBase) {
        validateId(id);
        return sequenceReader.getGapRegions(fromBase, toBase);
    }

    @Override
    public String getSequenceSlice(long id, long fromBase, long toBase) throws Exception {
        return getSequenceSlice(id, fromBase, toBase, SequenceRangeOption.WHOLE_SEQUENCE);
    }

    @Override
    public String getSequenceSlice(long id, long fromBase, long toBase, SequenceRangeOption option) throws Exception {
        validateId(id);
        return sequenceReader.getSequenceSliceString(fromBase, toBase, option);
    }

    @Override
    public Reader getSequenceSliceReader(long id, long fromBase, long toBase) throws Exception {
        return getSequenceSliceReader(id, fromBase, toBase, SequenceRangeOption.WHOLE_SEQUENCE);
    }

    @Override
    public Reader getSequenceSliceReader(long id, long fromBase, long toBase, SequenceRangeOption option)
            throws Exception {
        validateId(id);
        return sequenceReader.getSequenceSliceReader(fromBase, toBase, option);
    }

    @Override
    public SequenceIndex getSequenceIndex(long id) {
        validateId(id);
        return sequenceReader.getSequenceIndex();
    }

    @Override
    public SequenceInfoDTO exportReaderSettings() {
        SequenceIndex indexCopy = new SequenceIndex(sequenceReader.getSequenceIndex()); // deep copy
        // requires copy constructor on SequenceIndex

        SequenceInfoDTO dto = new SequenceInfoDTO();
        dto.filePath = sequenceReader.getFile().toPath();
        dto.sequenceFileFormat = sequenceReader.getSequenceFileFormat();
        dto.sequenceAlphabetSettings = sequenceReader.getSequenceAlphabet().exportAlphabetSettings();
        dto.sequenceIndexesMap = new HashMap<>();
        dto.sequenceIndexesMap.put(precodedId, indexCopy);
        dto.headerLines = null; // null for PLAIN_SEQUENCE format
        return dto;
    }

    @Override
    public void close() throws Exception {
        sequenceReader.close();
    }

    private void validateId(long id) {
        if (id != precodedId) {
            throw new IllegalArgumentException("No record for id :" + id);
        }
    }
}
