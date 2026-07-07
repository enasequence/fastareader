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

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.util.List;
import org.junit.jupiter.api.Test;
import uk.ac.ebi.embl.fastareader.sequenceutils.GapRegion;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;

/** End-to-end FASTA reading with the protein alphabet: X is the gap character, N is a normal residue. */
class FastaReaderProteinIntegrationTest {

    @Test
    void readsProteinFastaWithXGapsAndNResidues() throws Exception {
        File fasta = TestResources.file("fasta", "protein_example.txt");

        try (FastaReader service = new FastaReader(fasta, SequenceAlphabet.defaultProteinAlphabet())) {
            assertEquals(List.of(0L, 1L), service.getOrderedIds());

            // prot1: X X M N N K Y X X
            SequenceStats stats = service.getStats(0L);
            assertEquals(9, stats.totalBases());
            assertEquals(2, stats.leadingNsCount(), "leading X gap bases");
            assertEquals(2, stats.trailingNsCount(), "trailing X gap bases");
            assertEquals(5, stats.totalBasesWithoutNBases(), "MNNKY once edge X trimmed");

            // gap regions cover only the X runs; the internal NN (Asn) is NOT a gap
            assertEquals(List.of(new GapRegion(1, 2), new GapRegion(8, 9)), service.getGapRegions(0L));

            assertEquals("XXMNNKYXX", service.getSequenceSlice(0L, 1, 9, SequenceRangeOption.WHOLE_SEQUENCE));
            // edge-trimmed view: bases 1..5 map to the residues between the X runs
            assertEquals("MNNKY", service.getSequenceSlice(0L, 1, 5, SequenceRangeOption.WITHOUT_EDGE_N_BASES));

            // prot2 has a stop character and no gaps
            SequenceStats stats2 = service.getStats(1L);
            assertEquals(11, stats2.totalBases());
            assertEquals(0, stats2.leadingNsCount());
            assertEquals(0, stats2.trailingNsCount());
            assertTrue(service.getGapRegions(1L).isEmpty(), "no gaps in prot2");
            assertEquals("MKTAYIAKQR*", service.getSequenceSlice(1L, 1, 11, SequenceRangeOption.WHOLE_SEQUENCE));
        }
    }
}
