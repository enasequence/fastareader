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
package uk.ac.ebi.embl.fastareader.headerutils;

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.StandardCharsets;
import java.nio.file.StandardOpenOption;
import java.util.HashSet;
import java.util.Set;
import org.junit.jupiter.api.Test;
import uk.ac.ebi.embl.fastareader.FastaTestResources;

class HeaderLineDecoderTest {

    private static final byte LF = (byte) '\n';
    private static final byte CR = (byte) '\r';

    private static Set<Byte> endChars() {
        Set<Byte> s = new HashSet<>();
        s.add(LF);
        s.add(CR);
        return s;
    }

    private static FileChannel openRead(File f) throws IOException {
        return FileChannel.open(f.toPath(), StandardOpenOption.READ);
    }

    @Test
    void readsAsciiLineEndingWithLf_andAdvancesPositionPastLf() throws Exception {
        // File content: ">abc\nNEXT"
        File fasta = FastaTestResources.file("headerutils", "ascii_lf.txt");

        HeaderLineDecoder d = new HeaderLineDecoder(StandardCharsets.UTF_8, endChars(), 8);

        try (FileChannel ch = openRead(fasta)) {
            long from = 0;
            String line = d.readHeaderLine(ch, from);

            assertEquals(">abc", line, "Should read header without LF");
            assertEquals(5, ch.position(), "Position should advance past LF (byte index 4 + 1)");
        }
    }

    @Test
    void readsUtf8LineWhereMultibyteCharIsSplitAcrossBuffers() throws Exception {
        // Ensure BUFFER_SIZE is small so the multibyte sequence is split.
        // File content: ">A¢B\n"
        // UTF-8 bytes for € are E2 82 AC (3 bytes). With BUFFER_SIZE=4 we can force splits.
        File fasta = FastaTestResources.file("headerutils", "utf8_split.txt");

        HeaderLineDecoder d = new HeaderLineDecoder(StandardCharsets.UTF_8, endChars(), 4);

        try (FileChannel ch = openRead(fasta)) {
            String line = d.readHeaderLine(ch, 0);
            assertEquals(">A¢B", line);
        }
    }

    @Test
    void stopsAtCr_andAdvancesPositionPastCrOnly() throws Exception {
        // File content: ">abc\rNEXT"
        File fasta = FastaTestResources.file("headerutils", "ascii_cr.txt");

        HeaderLineDecoder d = new HeaderLineDecoder(StandardCharsets.UTF_8, endChars(), 8);

        try (FileChannel ch = openRead(fasta)) {
            String line = d.readHeaderLine(ch, 0);

            assertEquals(">abc", line);
            assertEquals(5, ch.position(), "Position should advance past CR (byte index 4 + 1)");
        }
    }

    @Test
    void crlfConsumesOnlyCr_andLeavesLfAsNextByte_sameAsCurrentBehavior() throws Exception {
        // File content: ">abc\r\nNEXT"
        File fasta = FastaTestResources.file("headerutils", "crlf.txt");

        HeaderLineDecoder d = new HeaderLineDecoder(StandardCharsets.UTF_8, endChars(), 8);

        try (FileChannel ch = openRead(fasta)) {
            String line = d.readHeaderLine(ch, 0);
            assertEquals(">abc", line);
            assertEquals(5, ch.position());
            String second = d.readHeaderLine(ch, ch.position());
            assertEquals("", second, "Because LF is treated as end char at start position");
        }
    }

    @Test
    void readsLineAtEofWithoutNewline_andPositionsAtFileEnd() throws Exception {
        // File content: ">no_newline_at_eof"
        File fasta = FastaTestResources.file("headerutils", "no_newline.txt");

        HeaderLineDecoder d = new HeaderLineDecoder(StandardCharsets.UTF_8, endChars(), 8);

        try (FileChannel ch = openRead(fasta)) {
            long size = ch.size();
            String line = d.readHeaderLine(ch, 0);

            assertEquals(">no_newline_at_eof", line);
            assertEquals(size, ch.position(), "Should move to EOF when no delimiter is found");
        }
    }

    @Test
    void returnsNullIfFromAtOrPastEof() throws Exception {
        File fasta = FastaTestResources.file("headerutils", "ascii_lf.txt");

        HeaderLineDecoder d = new HeaderLineDecoder(StandardCharsets.UTF_8, endChars(), 8);

        try (FileChannel ch = openRead(fasta)) {
            long size = ch.size();
            assertNull(d.readHeaderLine(ch, size), "from == fileSize should return null");
            assertNull(d.readHeaderLine(ch, size + 10), "from > fileSize should return null");
        }
    }

    @Test
    void canReadSecondHeaderFromOffset() throws Exception {
        // File content:
        // >first
        // >second
        File fasta = FastaTestResources.file("headerutils", "two_headers.txt");

        HeaderLineDecoder d = new HeaderLineDecoder(StandardCharsets.UTF_8, endChars(), 4);

        try (FileChannel ch = openRead(fasta)) {
            String first = d.readHeaderLine(ch, 0);
            assertEquals(">first", first);
            long secondFrom = ch.position();
            String second = d.readHeaderLine(ch, secondFrom);
            assertEquals(">second", second);
        }
    }

    @Test
    void throwsOnMalformedUtf8_whenReportIsEnabled() throws Exception {
        // File content: ">bad \xC3\x28\n" (bytes that are not valid UTF-8 in the header before newline)
        // In UTF-8, 0xC3 must be followed by 0x80..0xBF; 0x28 '(' is invalid.
        File fasta = FastaTestResources.file("headerutils", "malformed_utf8.txt");

        HeaderLineDecoder d = new HeaderLineDecoder(StandardCharsets.UTF_8, endChars(), 8);

        try (FileChannel ch = openRead(fasta)) {
            assertThrows(CharacterCodingException.class, () -> d.readHeaderLine(ch, 0));
        }
    }
}
