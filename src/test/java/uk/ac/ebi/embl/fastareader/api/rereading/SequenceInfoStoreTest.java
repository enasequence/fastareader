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

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import uk.ac.ebi.embl.fastareader.SequenceFileFormat;
import uk.ac.ebi.embl.fastareader.TestResources;
import uk.ac.ebi.embl.fastareader.api.SequenceFormatReader;
import uk.ac.ebi.embl.fastareader.api.SequenceFormatReaderFactory;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;

public class SequenceInfoStoreTest {

    @TempDir
    Path tempDir;

    private SequenceInfoStore store;

    @BeforeEach
    void setUp() {
        store = new SequenceInfoStore();
    }

    @Test
    void fastaSequenceInfoSurvivesJsonRoundTrip() throws Exception {
        File sequenceFile = TestResources.file("fasta", "differing_headerline_fasta.txt");
        List<SequenceInfoDTO> infos = new ArrayList<>();

        try (SequenceFormatReader service =
                SequenceFormatReaderFactory.readFasta(sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {
            infos.add(service.exportReaderSettings());
        }

        Path outputPath = tempDir.resolve("sequenceInfo.json");
        store.write(infos, outputPath);

        List<SequenceInfoDTO> reImportedInfos = store.read(outputPath);

        assertEquals(infos.size(), reImportedInfos.size());
        SequenceInfoDTO original = infos.get(0);
        SequenceInfoDTO reImported = reImportedInfos.get(0);

        assertEquals(original.filePath, reImported.filePath);
        assertEquals(original.sequenceFileFormat, reImported.sequenceFileFormat);
        assertEquals(original.sequenceAlphabetSettings, reImported.sequenceAlphabetSettings);
        assertEquals(original.sequenceIndexesMap, reImported.sequenceIndexesMap);
        assertEquals(original.headerLines, reImported.headerLines);
    }

    @Test
    void plainSequenceInfoSurvivesJsonRoundTrip() throws Exception {
        File sequenceFile = TestResources.file("sequence", "only_ns_example.txt");
        List<SequenceInfoDTO> infos = new ArrayList<>();

        try (SequenceFormatReader service = SequenceFormatReaderFactory.readPlainSequence(
                sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {
            infos.add(service.exportReaderSettings());
        }

        Path outputPath = tempDir.resolve("sequenceInfo.json");
        store.write(infos, outputPath);

        List<SequenceInfoDTO> reImportedInfos = store.read(outputPath);

        assertEquals(infos.size(), reImportedInfos.size());
        SequenceInfoDTO original = infos.get(0);
        SequenceInfoDTO reImported = reImportedInfos.get(0);

        assertEquals(original.filePath, reImported.filePath);
        assertEquals(original.sequenceFileFormat, reImported.sequenceFileFormat);
        assertEquals(original.sequenceAlphabetSettings, reImported.sequenceAlphabetSettings);
        assertEquals(original.sequenceIndexesMap, reImported.sequenceIndexesMap);
        assertNull(reImported.headerLines);
    }

    @Test
    void multipleSequenceInfosSurviveJsonRoundTrip() throws Exception {
        File fastaFile = TestResources.file("fasta", "differing_headerline_fasta.txt");
        File sequenceFile = TestResources.file("sequence", "only_ns_example.txt");
        List<SequenceInfoDTO> infos = new ArrayList<>();

        try (SequenceFormatReader service =
                SequenceFormatReaderFactory.readFasta(fastaFile, SequenceAlphabet.defaultNucleotideAlphabet())) {
            infos.add(service.exportReaderSettings());
        }
        try (SequenceFormatReader service = SequenceFormatReaderFactory.readPlainSequence(
                sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {
            infos.add(service.exportReaderSettings());
        }

        Path outputPath = tempDir.resolve("sequenceInfo.json");
        store.write(infos, outputPath);

        List<SequenceInfoDTO> reImportedInfos = store.read(outputPath);

        assertEquals(2, reImportedInfos.size());
        assertEquals(SequenceFileFormat.FASTA, reImportedInfos.get(0).sequenceFileFormat);
        assertEquals(SequenceFileFormat.PLAIN_SEQUENCE, reImportedInfos.get(1).sequenceFileFormat);
    }

    @Test
    void readThrowsOnInvalidPath() {
        Path invalidPath = tempDir.resolve("nonexistent.json");
        assertThrows(IllegalArgumentException.class, () -> store.read(invalidPath));
    }

    @Test
    void writeCreatesParentDirectoriesIfMissing() throws Exception {
        File sequenceFile = TestResources.file("fasta", "differing_headerline_fasta.txt");
        List<SequenceInfoDTO> infos = new ArrayList<>();

        try (SequenceFormatReader service =
                SequenceFormatReaderFactory.readFasta(sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {
            infos.add(service.exportReaderSettings());
        }

        Path nestedPath = tempDir.resolve("nested/dirs/sequenceInfo.json");
        store.write(infos, nestedPath);

        assertTrue(Files.exists(nestedPath));
    }

    @Test
    void writeToPathWithNullParentDoesNotNpe() throws Exception {
        // Paths.get("bare.json").getParent() == null; previously caused NullPointerException
        Path bare = Paths.get("seq_info_test_bare.json");
        assertNull(bare.getParent());
        try {
            SequenceInfoStore.write(List.of(), bare);
            assertTrue(Files.exists(bare));
        } finally {
            Files.deleteIfExists(bare);
        }
    }
}
