package uk.ac.ebi.embl.fastareader.api;

import org.junit.jupiter.api.Test;
import uk.ac.ebi.embl.fastareader.*;
import uk.ac.ebi.embl.fastareader.api.rereading.SequenceInfoDTO;
import uk.ac.ebi.embl.fastareader.exception.SequenceFileException;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;

import java.io.File;
import java.io.IOException;
import java.util.Optional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

public class SequenceFormatReaderExportTest {

    @Test
    void checkIfFastaReImportWorks() throws Exception {
        File sequenceFile = TestResources.file("fasta", "differing_headerline_fasta.txt");
        try (SequenceFormatReader service = SequenceFormatReaderFactory.readFasta(sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {

            SequenceStats stats = service.getStats(0L);
            String sequence = service.getSequenceSlice(0L, 1, stats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);

            //try reinitialising sequence reader from the exported settings and compare the internal information
            SequenceInfoDTO info = service.exportReaderSettings();
            service.close();
            try (SequenceFormatReader newService = SequenceFormatReaderFactory.readBySequenceInfo(info)) {
                assertEquals(newService.getFile().toPath(), info.filePath);
                assertEquals(newService.getSequenceFileFormat(), info.sequenceFileFormat);
                assertEquals(newService.getSequenceAlphabet().describeAllowed(), info.sequenceAlphabet.describeAllowed());
                //check headerlines
                assertEquals(newService.getOrderedIds().size(), info.headerLines.size());
                assertEquals(newService.getHeaderline(0L).get(), info.headerLines.get(0L));
                //check indexes
                assertEquals(newService.getOrderedIds().size(), info.sequenceIndexesMap.size());
                assertEquals(newService.getSequenceIndex(0L), info.sequenceIndexesMap.get(0L));
                //check sequence is read correctly
                SequenceStats newstats = newService.getStats(0L);
                String newSequence = newService.getSequenceSlice(0L, 1, newstats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(newSequence, sequence);
            }

        }
    }
    @Test
    void checkIfSequenceReImportWorks() throws Exception {
        File sequenceFile = TestResources.file("sequence", "only_ns_example.txt");
        try (SequenceFormatReader service = SequenceFormatReaderFactory.readPlainSequence(sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet())) {

            SequenceStats stats = service.getStats(0L);
            String sequence = service.getSequenceSlice(0L, 1, stats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
            assertEquals(43, stats.leadingNsCount(), "leading Ns");
            assertEquals(43, stats.trailingNsCount(), "trailing Ns");
            assertEquals(43, stats.baseCount().get('N'), "total Ns");

            //try reinitialising sequence reader from the exported settings and compare the internal information
            SequenceInfoDTO info = service.exportReaderSettings();
            service.close();
            try (SequenceFormatReader newService = SequenceFormatReaderFactory.readBySequenceInfo(info)) {
                assertEquals(newService.getFile().toPath(), info.filePath);
                assertEquals(newService.getSequenceFileFormat(), info.sequenceFileFormat);
                assertEquals(newService.getSequenceAlphabet().describeAllowed(), info.sequenceAlphabet.describeAllowed());
                //check headerlines
                assertNull(info.headerLines);
                assertEquals(newService.getHeaderline(0L), Optional.empty());
                //check indexes
                assertEquals(1, info.sequenceIndexesMap.size());
                assertEquals(newService.getSequenceIndex(0L), info.sequenceIndexesMap.get(0L));
                //check sequence reading works
                SequenceStats newstats = newService.getStats(0L);
                String newSequence = newService.getSequenceSlice(0L, 1, newstats.totalBases(), SequenceRangeOption.WHOLE_SEQUENCE);
                assertEquals(newSequence, sequence);
            }

        }
    }

}
