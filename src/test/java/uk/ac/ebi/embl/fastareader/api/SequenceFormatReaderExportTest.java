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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.File;
import java.util.Optional;
import org.junit.jupiter.api.Test;
import uk.ac.ebi.embl.fastareader.*;
import uk.ac.ebi.embl.fastareader.api.rereading.SequenceInfoDTO;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;

public class SequenceFormatReaderExportTest {

    @Test
    void checkIfFastaReImportWorks() throws Exception {
        File sequenceFile = TestResources.file("fasta", "differing_headerline_fasta.txt");
        try (SequenceFormatReader service =
                SequenceFormatReaderFactory.readFasta(sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {

            SequenceStats stats = service.getStats(0L);
            String sequence = service.getSequenceSlice(0L, 1, stats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);

            // try reinitialising sequence reader from the exported settings and compare the internal information
            SequenceInfoDTO info = service.exportReaderSettings();
            service.close();

            try (SequenceFormatReader newService = SequenceFormatReaderFactory.readBySequenceInfo(info)) {
                assertEquals(newService.getFile().toPath(), info.filePath);
                assertEquals(newService.getSequenceFileFormat(), info.sequenceFileFormat);
                assertEquals(
                        newService.getSequenceAlphabet().describeAllowed(),
                        new SequenceAlphabet(info.sequenceAlphabetSettings).describeAllowed());
                // check headerlines
                assertEquals(newService.getOrderedIds().size(), info.headerLines.size());
                assertEquals(newService.getHeaderline(0L).get(), info.headerLines.get(0L));
                // check indexes
                assertEquals(newService.getOrderedIds().size(), info.sequenceIndexesMap.size());
                assertEquals(newService.getSequenceIndex(0L), info.sequenceIndexesMap.get(0L));
                // check sequence is read correctly
                SequenceStats newstats = newService.getStats(0L);
                String newSequence =
                        newService.getSequenceSlice(0L, 1, newstats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(newSequence, sequence);
            }
        }
    }

    @Test
    void checkIfSequenceReImportWorks() throws Exception {
        File sequenceFile = TestResources.file("sequence", "only_ns_example.txt");
        try (SequenceFormatReader service = SequenceFormatReaderFactory.readPlainSequence(
                sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {

            SequenceStats stats = service.getStats(0L);
            String sequence = service.getSequenceSlice(0L, 1, stats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals(43, stats.leadingNsCount(), "leading Ns");
            assertEquals(43, stats.trailingNsCount(), "trailing Ns");
            assertEquals(43, stats.baseCount().get('N'), "total Ns");

            // try reinitialising sequence reader from the exported settings and compare the internal information
            SequenceInfoDTO info = service.exportReaderSettings();
            service.close();
            try (SequenceFormatReader newService = SequenceFormatReaderFactory.readBySequenceInfo(info)) {
                assertEquals(newService.getFile().toPath(), info.filePath);
                assertEquals(newService.getSequenceFileFormat(), info.sequenceFileFormat);
                assertEquals(
                        newService.getSequenceAlphabet().describeAllowed(),
                        new SequenceAlphabet(info.sequenceAlphabetSettings).describeAllowed());
                // check headerlines
                assertNull(info.headerLines);
                assertEquals(newService.getHeaderline(0L), Optional.empty());
                // check indexes
                assertEquals(1, info.sequenceIndexesMap.size());
                assertEquals(newService.getSequenceIndex(0L), info.sequenceIndexesMap.get(0L));
                // check sequence reading works
                SequenceStats newstats = newService.getStats(0L);
                String newSequence =
                        newService.getSequenceSlice(0L, 1, newstats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(newSequence, sequence);
            }
        }
    }

    @Test
    void fastaReaderSurvivesJsonRoundTrip() throws Exception {
        File sequenceFile = TestResources.file("fasta", "differing_headerline_fasta.txt");
        try (SequenceFormatReader service =
                SequenceFormatReaderFactory.readFasta(sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {

            SequenceStats stats = service.getStats(0L);
            String sequence = service.getSequenceSlice(0L, 1, stats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);

            SequenceInfoDTO info = service.exportReaderSettings();
            service.close();

            // write to json and read back
            ObjectMapper mapper = new ObjectMapper();
            String json = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(info);
            // System.out.println(json);
            SequenceInfoDTO reImportedInfo = mapper.readValue(json, SequenceInfoDTO.class);

            try (SequenceFormatReader newService = SequenceFormatReaderFactory.readBySequenceInfo(reImportedInfo)) {
                assertEquals(newService.getFile().toPath(), reImportedInfo.filePath);
                assertEquals(newService.getSequenceFileFormat(), reImportedInfo.sequenceFileFormat);
                assertEquals(
                        newService.getSequenceAlphabet().describeAllowed(),
                        new SequenceAlphabet(reImportedInfo.sequenceAlphabetSettings).describeAllowed());
                // check headerlines
                assertEquals(newService.getOrderedIds().size(), reImportedInfo.headerLines.size());
                assertEquals(newService.getHeaderline(0L).get(), reImportedInfo.headerLines.get(0L));
                // check indexes
                assertEquals(newService.getOrderedIds().size(), reImportedInfo.sequenceIndexesMap.size());
                assertEquals(newService.getSequenceIndex(0L), reImportedInfo.sequenceIndexesMap.get(0L));
                // check sequence is read correctly
                SequenceStats newStats = newService.getStats(0L);
                String newSequence =
                        newService.getSequenceSlice(0L, 1, newStats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(newSequence, sequence);
            }
        }
    }

    @Test
    void plainSequenceReaderSurvivesJsonRoundTrip() throws Exception {
        File sequenceFile = TestResources.file("sequence", "only_ns_example.txt");
        try (SequenceFormatReader service = SequenceFormatReaderFactory.readPlainSequence(
                sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {

            SequenceStats stats = service.getStats(0L);
            String sequence = service.getSequenceSlice(0L, 1, stats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals(43, stats.leadingNsCount(), "leading Ns");
            assertEquals(43, stats.trailingNsCount(), "trailing Ns");
            assertEquals(43, stats.baseCount().get('N'), "total Ns");

            SequenceInfoDTO info = service.exportReaderSettings();
            service.close();

            // write to json and read back
            ObjectMapper mapper = new ObjectMapper();
            String json = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(info);
            // System.out.println(json);
            SequenceInfoDTO reImportedInfo = mapper.readValue(json, SequenceInfoDTO.class);

            try (SequenceFormatReader newService = SequenceFormatReaderFactory.readBySequenceInfo(reImportedInfo)) {
                assertEquals(newService.getFile().toPath(), reImportedInfo.filePath);
                assertEquals(newService.getSequenceFileFormat(), reImportedInfo.sequenceFileFormat);
                assertEquals(
                        newService.getSequenceAlphabet().describeAllowed(),
                        new SequenceAlphabet(reImportedInfo.sequenceAlphabetSettings).describeAllowed());
                // check headerlines
                assertNull(reImportedInfo.headerLines);
                assertEquals(newService.getHeaderline(0L), Optional.empty());
                // check indexes
                assertEquals(1, reImportedInfo.sequenceIndexesMap.size());
                assertEquals(newService.getSequenceIndex(0L), reImportedInfo.sequenceIndexesMap.get(0L));
                // check sequence reading works
                SequenceStats newStats = newService.getStats(0L);
                String newSequence =
                        newService.getSequenceSlice(0L, 1, newStats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(newSequence, sequence);
            }
        }
    }
}
