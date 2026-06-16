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
package uk.ac.ebi.embl.fastareader.encoding;

import static org.junit.jupiter.api.Assertions.*;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import uk.ac.ebi.embl.fastareader.TestResources;

class Utf8DetectorTest {

    @Test
    void asciiIsUtf8() throws IOException {
        Path p = TestResources.path("encoding", "ascii.txt");
        assertTrue(Utf8Detector.isProbablyUtf8(p));
    }

    @Test
    void validUtf8MultibyteIsUtf8() throws IOException {
        Path p = TestResources.path("encoding", "utf8.txt");
        assertTrue(Utf8Detector.isProbablyUtf8(p));
    }

    @Test
    void utf8WithBomIsUtf8() throws IOException {
        Path p = TestResources.path("encoding", "utf8_bom.txt");
        assertTrue(Utf8Detector.isProbablyUtf8(p));
    }

    @Test
    void detectsNonUtf8EncodingCorrectly() throws IOException {
        Path p = TestResources.path(
                "encoding", "non_utf8.txt"); // contains 0x80 by itself, which is an invalid UTF-8 start byte.
        Path p1 = TestResources.path(
                "encoding",
                "non_utf8_1.txt"); // contains 0xE9 by itself - in ISO-8859-1 this is "é", but as a lone byte it's
        // invalid UTF-8.
        Path p2 = TestResources.path("encoding", "non_utf8_2.txt");

        assertFalse(Utf8Detector.isProbablyUtf8(p));
        assertFalse(Utf8Detector.isProbablyUtf8(p1));
        assertFalse(Utf8Detector.isProbablyUtf8(p2));
    }

    // ---- InputStream overload (used by the BGZF-aware UTF-8 gate) ----

    @Test
    void streamOverloadValidAsciiPasses() throws IOException {
        byte[] content = "ACGTACGT\nACGT\n".getBytes(StandardCharsets.US_ASCII);
        assertTrue(Utf8Detector.isProbablyUtf8(new ByteArrayInputStream(content), content.length));
    }

    @Test
    void streamOverloadValidMultibyteUtf8Passes() throws IOException {
        // U+00E9 encoded as 2-byte UTF-8 sequence (0xC3 0xA9)
        byte[] content = {0x41, (byte) 0xC3, (byte) 0xA9, 0x0A};
        assertTrue(Utf8Detector.isProbablyUtf8(new ByteArrayInputStream(content), content.length));
    }

    @Test
    void streamOverloadNonUtf8ByteFails() throws IOException {
        // 0x80 alone is an invalid UTF-8 start byte
        byte[] content = {0x41, 0x43, (byte) 0x80, 0x54};
        assertFalse(Utf8Detector.isProbablyUtf8(new ByteArrayInputStream(content), content.length));
    }

    @Test
    void streamOverloadRespectsMaxBytesLimit() throws IOException {
        // First maxBytes bytes are valid ASCII, last byte is 0x80 (invalid UTF-8)
        int limit = 8;
        byte[] content = new byte[limit + 1];
        for (int i = 0; i < limit; i++) content[i] = 0x41; // 'A'
        content[limit] = (byte) 0x80;
        assertTrue(
                Utf8Detector.isProbablyUtf8(new ByteArrayInputStream(content), limit),
                "should pass when maxBytes covers only the valid prefix");
        assertFalse(
                Utf8Detector.isProbablyUtf8(new ByteArrayInputStream(content), limit + 1),
                "should fail when maxBytes covers the invalid byte");
    }

    // ---- Path-based overload (existing tests) ----

    @Test
    void respectsMaxBytesLimit() throws IOException {
        Path p = TestResources.path(
                "encoding",
                "large_non_utf8.txt"); // 1Mb file that is utf8 for the first 1Mb, and non utf8 as the last char
        int limit = 1024 * 1024; // 1 MiB limit to test the file

        assertTrue(Utf8Detector.isProbablyUtf8(p, limit)); // first 1Mb should be detected as Utf8
        assertFalse(Utf8Detector.isProbablyUtf8(p, limit + 1)); // if detecting later characters, should fail
    }
}
