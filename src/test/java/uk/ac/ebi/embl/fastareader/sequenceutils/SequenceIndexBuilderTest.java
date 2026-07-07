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

import static org.junit.jupiter.api.Assertions.*;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import uk.ac.ebi.embl.fastareader.exception.SequenceReadingException;
import uk.ac.ebi.embl.fastareader.io.FileChannelByteReader;
import uk.ac.ebi.embl.fastareader.io.SeekableByteReader;

public class SequenceIndexBuilderTest {

    private static final byte GT = (byte) '>';

    @TempDir
    Path tempDir;

    private static SeekableByteReader openRead(Path p) throws IOException {
        return new FileChannelByteReader(p.toFile());
    }

    private static Path writeAscii(Path dir, String filename, String content) throws IOException {
        Path p = dir.resolve(filename);
        Files.write(p, content.getBytes(StandardCharsets.US_ASCII));
        return p;
    }

    @Test
    void buildsIndexCorrectly() throws Exception {
        // Layout (US-ASCII):
        // >ID1 | {"d":"x"}\n
        // NNAC\n
        // acgt\n
        // ttnN\n
        // \n
        // \t\n
        // \n
        // >NEXT\n
        String header = ">ID1 | {\"d\":\"x\"}\n";
        String l1 = "NNAC\n"; // leading N=2
        String l2 = "acgt\n";
        String l3 = "ttnN\n"; // trailing N=2
        String empties = "\n\n\n";
        String nextHead = ">NEXT\n";

        String fasta = header + l1 + l2 + l3 + empties + nextHead;
        Path p = writeAscii(tempDir, "idx1.fa", fasta);

        try (SeekableByteReader ch = openRead(p)) {
            long seqStartPos = header.getBytes(StandardCharsets.US_ASCII).length; // first byte after header line

            SequenceAlphabet alpha = SequenceAlphabet.defaultNucleotideAlphabet();
            SequenceIndexBuilder sib = new SequenceIndexBuilder(ch, alpha, Optional.of(GT));

            long beforePos = ch.position(); // should remain unchanged
            SequenceIndex idx = sib.buildFrom(seqStartPos);
            long afterPos = ch.position();

            // builder must not touch channel.position()
            assertEquals(beforePos, afterPos, "builder must not change channel.position()");

            // Lines: only 3 sequence lines; empties ignored
            List<LineEntry> lines = idx.linesView();
            assertEquals(3, lines.size(), "only non-empty sequence lines must be indexed");

            // Base numbering should be contiguous across lines (4 bases per line)
            assertEquals(1, lines.get(0).baseStart);
            assertEquals(4, lines.get(0).baseEnd);
            assertEquals(5, lines.get(1).baseStart);
            assertEquals(8, lines.get(1).baseEnd);
            assertEquals(9, lines.get(2).baseStart);
            assertEquals(12, lines.get(2).baseEnd);

            // Byte math: each line has 4 letters; byteEndExclusive = lastBaseByte + 1
            long l1Start = seqStartPos; // begins right after header line
            long l1EndEx = l1Start + 4;
            long l2Start = l1EndEx + 1; // + LF between lines
            long l2EndEx = l2Start + 4;
            long l3Start = l2EndEx + 1;
            long l3EndEx = l3Start + 4;

            assertEquals(l1Start, lines.get(0).byteStart);
            assertEquals(l1EndEx, lines.get(0).byteEndExclusive);

            assertEquals(l2Start, lines.get(1).byteStart);
            assertEquals(l2EndEx, lines.get(1).byteEndExclusive);

            assertEquals(l3Start, lines.get(2).byteStart);
            assertEquals(l3EndEx, lines.get(2).byteEndExclusive);

            // first/last base bytes
            assertEquals(l1Start, idx.firstBaseByte);
            assertEquals(l3EndEx - 1, idx.lastBaseByte);

            // Edge N counting: only first and last lines are inspected
            assertEquals(2, idx.startNBasesCount, "leading Ns only from first sequence line");
            assertEquals(2, idx.endNBasesCount, "trailing Ns only from last sequence line");
            // assert total count of N bases is correctly done
            assertEquals(4, idx.caseInsensitiveBaseCount.get('N'));
            assertEquals(2, idx.caseInsensitiveBaseCount.get('A'));
            assertEquals(2, idx.caseInsensitiveBaseCount.get('C'));
            assertEquals(3, idx.caseInsensitiveBaseCount.get('T'));
            assertEquals(1, idx.caseInsensitiveBaseCount.get('G'));

            // nextHeaderByte should point to '>' of NEXT header
            long expectedNextHeader = header.length() + l1.length() + l2.length() + l3.length() + empties.length();
            assertEquals(expectedNextHeader, idx.nextHeaderByte);
        }
    }

    @Test
    void buildsIndexCorrectlyTest2() throws Exception {
        String header = ">ID2\n";
        String l1 = "NNxx".replace('x', 'A') + "\n";
        String l2 = "gggg\n";
        String next = ">H2\n";

        String fasta = header + l1 + l2 + next;
        Path p = writeAscii(tempDir, "idx2.fa", fasta);

        try (SeekableByteReader ch = openRead(p)) {
            long seqStart = header.getBytes(StandardCharsets.US_ASCII).length;
            SequenceIndexBuilder sib =
                    new SequenceIndexBuilder(ch, SequenceAlphabet.defaultNucleotideAlphabet(), Optional.of(GT));

            SequenceIndex idx = sib.buildFrom(seqStart);

            // Two non-empty lines only
            assertEquals(2, idx.linesView().size());

            // leading Ns counted only on first line (here: 2)
            assertEquals(2, idx.startNBasesCount);
            // no trailing Ns on last line (all 'g')
            assertEquals(0, idx.endNBasesCount);
            // assert total count of N bases is correctly done
            assertEquals(2, idx.caseInsensitiveBaseCount.get('N'));
            assertEquals(4, idx.caseInsensitiveBaseCount.get('G'));
            assertEquals(2, idx.caseInsensitiveBaseCount.get('A'));

            // nextHeader should be at the '>' byte of H2
            long expectedNext = fasta.lastIndexOf(">H2\n"); // ascii index
            assertEquals(expectedNext, idx.nextHeaderByte);
        }
    }

    @Test
    void readsOnlyNsCorrectly() throws Exception {
        String header = ">ID2\n";
        String l1 = "NNnnnnnNNnnnnNNN\n";
        String l2 = "nnnnNNNNnnnNNNNnn\n";
        String next = ">H2\n";

        String fasta = header + l1 + l2 + next;
        Path p = writeAscii(tempDir, "idx3.fa", fasta);

        try (SeekableByteReader ch = openRead(p)) {
            long seqStart = header.getBytes(StandardCharsets.US_ASCII).length;
            SequenceIndexBuilder sib =
                    new SequenceIndexBuilder(ch, SequenceAlphabet.defaultNucleotideAlphabet(), Optional.of(GT));

            SequenceIndex idx = sib.buildFrom(seqStart);

            // Two non-empty lines only
            assertEquals(2, idx.linesView().size());

            // leading Ns until the end
            assertEquals(33, idx.startNBasesCount);
            // trailing Ns until the start
            assertEquals(33, idx.endNBasesCount);
            // assert total count of N bases is correctly done
            assertEquals(idx.totalBases(), idx.caseInsensitiveBaseCount.get('N'));

            // nextHeader should be at the '>' byte of H2
            long expectedNext = fasta.lastIndexOf(">H2\n"); // ascii index
            assertEquals(expectedNext, idx.nextHeaderByte);
        }
    }

    @Test
    void ignoresEmptyLinesCorrectly() throws Exception {
        String header = ">ID3\n";
        String l1 = "NACG\n"; // leading N = 1
        String l2 = "NNNN\n"; // middle line of Ns — must NOT affect start/end N counts
        String blanks = "\n\n";
        String l3 = "GGGn\n"; // trailing n = 1
        String next = ">K\n";

        String fasta = header + l1 + l2 + blanks + l3 + next;
        Path p = writeAscii(tempDir, "idx4.fa", fasta);

        try (SeekableByteReader ch = openRead(p)) {
            long seqStart = header.getBytes(StandardCharsets.US_ASCII).length;
            SequenceIndexBuilder sib =
                    new SequenceIndexBuilder(ch, SequenceAlphabet.defaultNucleotideAlphabet(), Optional.of(GT));

            long before = ch.position();
            SequenceIndex idx = sib.buildFrom(seqStart);
            long after = ch.position();
            assertEquals(before, after, "builder must not move channel position");

            // three non-empty sequence lines: l1, l2, l3
            assertEquals(3, idx.linesView().size());

            // Edge N counts: only first and last lines considered
            assertEquals(1, idx.startNBasesCount, "only first line leading Ns");
            assertEquals(1, idx.endNBasesCount, "only last line trailing Ns");
            // assert total count of N bases is correctly done
            assertEquals(6, idx.caseInsensitiveBaseCount.get('N'));
            assertEquals(4, idx.caseInsensitiveBaseCount.get('G'));
            assertEquals(1, idx.caseInsensitiveBaseCount.get('A'));
            assertEquals(1, idx.caseInsensitiveBaseCount.get('C'));

            // Middle line of Ns shouldn't change edge counts
            assertEquals(idx.linesView().get(1).lengthBases(), 4);

            // Total base numbering should be contiguous: 4 + 4 + 4 = 12
            assertEquals(12, idx.totalBases());
        }
    }

    @Test
    void readsJustSequenceCorrectly() throws Exception {
        String l1 = "NACG\n"; // leading N = 1
        String l2 = "NNNN\n"; // middle line of Ns — must NOT affect start/end N counts
        String blanks = "\n\n";
        String l3 = "GGGn\n"; // trailing n = 1

        String sequence = l1 + l2 + blanks + l3;
        Path p = writeAscii(tempDir, "idx5.fa", sequence);

        try (SeekableByteReader ch = openRead(p)) {
            SequenceIndexBuilder sib =
                    new SequenceIndexBuilder(ch, SequenceAlphabet.defaultNucleotideAlphabet(), Optional.empty());

            long before = ch.position();
            SequenceIndex idx = sib.buildFrom(0);
            long after = ch.position();
            assertEquals(before, after, "builder must not move channel position");

            // three non-empty sequence lines: l1, l2, l3
            assertEquals(3, idx.linesView().size());

            // Edge N counts: only first and last lines considered
            assertEquals(1, idx.startNBasesCount, "only first line leading Ns");
            assertEquals(1, idx.endNBasesCount, "only last line trailing Ns");
            // assert total count of N bases is correctly done
            assertEquals(6, idx.caseInsensitiveBaseCount.get('N'));
            assertEquals(4, idx.caseInsensitiveBaseCount.get('G'));
            assertEquals(1, idx.caseInsensitiveBaseCount.get('A'));
            assertEquals(1, idx.caseInsensitiveBaseCount.get('C'));

            // Middle line of Ns shouldn't change edge counts
            assertEquals(idx.linesView().get(1).lengthBases(), 4);

            // Total base numbering should be contiguous: 4 + 4 + 4 = 12
            assertEquals(12, idx.totalBases());
        }
    }

    @Test
    void recordsGapRegionsAcrossLinesAndEmptyLines() throws Exception {
        String header = ">ID4\n";
        String l1 = "ACNN\n";
        String blanks = "\n\n";
        String l2 = "nnGTN\n";
        String l3 = "NAC\n";
        String next = ">K\n";

        String fasta = header + l1 + blanks + l2 + l3 + next;
        Path p = writeAscii(tempDir, "idx_gap.fa", fasta);

        try (SeekableByteReader ch = openRead(p)) {
            long seqStart = header.getBytes(StandardCharsets.US_ASCII).length;
            SequenceIndexBuilder sib =
                    new SequenceIndexBuilder(ch, SequenceAlphabet.defaultNucleotideAlphabet(), Optional.of(GT));

            SequenceIndex idx = sib.buildFrom(seqStart);

            assertEquals(List.of(new GapRegion(3, 6), new GapRegion(9, 10)), idx.gapRegionsView());
            assertEquals(List.of(new GapRegion(3, 6)), idx.gapRegionsView(1, 8));
            assertEquals(List.of(new GapRegion(9, 10)), idx.gapRegionsView(7, 11));
        }
    }

    @Test
    void recordsSingleGapWhenSequenceIsOnlyNs() throws Exception {
        String sequence = "NNnn\nNN";
        Path p = writeAscii(tempDir, "idx_only_gap.txt", sequence);

        try (SeekableByteReader ch = openRead(p)) {
            SequenceIndexBuilder sib =
                    new SequenceIndexBuilder(ch, SequenceAlphabet.defaultNucleotideAlphabet(), Optional.empty());

            SequenceIndex idx = sib.buildFrom(0);

            assertEquals(List.of(new GapRegion(1, 6)), idx.gapRegionsView());
        }
    }

    @Test
    void proteinAlphabetTreatsXAsGapAndNAsResidue() throws Exception {
        // X X M N N K Y | X X   (bases 1..9)
        // Leading gap = XX (1-2), trailing gap = XX (8-9). The internal NN are Asparagine residues,
        // NOT gaps, so they must not appear as a gap region.
        String sequence = "XXMNNKY\nXX";
        Path p = writeAscii(tempDir, "idx_protein_gap.txt", sequence);

        try (SeekableByteReader ch = openRead(p)) {
            SequenceIndexBuilder sib =
                    new SequenceIndexBuilder(ch, SequenceAlphabet.defaultProteinAlphabet(), Optional.empty());

            SequenceIndex idx = sib.buildFrom(0);

            assertEquals(List.of(new GapRegion(1, 2), new GapRegion(8, 9)), idx.gapRegionsView());
            assertEquals(2, idx.startNBasesCount, "leading X gap bases");
            assertEquals(2, idx.endNBasesCount, "trailing X gap bases");
            assertEquals(9, idx.totalBases());
        }
    }

    @Test
    void recordsGapRegionsAcrossLineBreaksAndTabsWithoutOffByOneErrors() throws Exception {
        String sequence = "ACGTNN\t\n" + "NNNCGTN\n" + "NNACTGNN";
        Path p = writeAscii(tempDir, "idx_gap_off_by_one.txt", sequence);

        try (SeekableByteReader ch = openRead(p)) {
            SequenceAlphabet alphabet = new SequenceAlphabet("ACGTNacgtn", "\n\r\t");
            SequenceIndexBuilder sib = new SequenceIndexBuilder(ch, alphabet, Optional.empty());

            SequenceIndex idx = sib.buildFrom(0);

            assertEquals(
                    List.of(new GapRegion(5, 9), new GapRegion(13, 15), new GapRegion(20, 21)), idx.gapRegionsView());
        }
    }

    @Test
    void doesNotTolerateBreaksWhenReadingContinuousSequence() throws Exception {
        String l1 = "NACG\n"; // leading N = 1
        String l2 = "N>NNN\n"; // middle line of Ns — must NOT affect start/end N counts
        String blanks = "\n\n";
        String l3 = "GGGn\n"; // trailing n = 1

        String sequence = l1 + l2 + blanks + l3;
        Path p = writeAscii(tempDir, "idx6.fa", sequence);

        try (SeekableByteReader ch = openRead(p)) {
            SequenceIndexBuilder sib =
                    new SequenceIndexBuilder(ch, SequenceAlphabet.defaultNucleotideAlphabet(), Optional.empty());

            assertThrows(SequenceReadingException.class, () -> sib.buildFrom(0));
        }
    }
}
