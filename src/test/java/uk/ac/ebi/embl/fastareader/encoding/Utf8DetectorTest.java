package uk.ac.ebi.embl.fastareader.encoding;
import org.junit.jupiter.api.Test;
import uk.ac.ebi.embl.fastareader.FastaTestResources;

import java.io.IOException;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;

class Utf8DetectorTest {

    @Test
    void asciiIsUtf8() throws IOException {
        Path p = FastaTestResources.path("encoding", "ascii.txt");
        assertTrue(Utf8Detector.isProbablyUtf8(p));
    }

    @Test
    void validUtf8MultibyteIsUtf8() throws IOException {
        Path p = FastaTestResources.path("encoding", "utf8.txt");
        assertTrue(Utf8Detector.isProbablyUtf8(p));
    }

    @Test
    void utf8WithBomIsUtf8() throws IOException {
        Path p = FastaTestResources.path("encoding", "utf8_bom.txt");
        assertTrue(Utf8Detector.isProbablyUtf8(p));
    }

    @Test
    void detectsNonUtf8EncodingCorrectly() throws IOException {
        Path p = FastaTestResources.path("encoding", "non_utf8.txt"); // contains 0x80 by itself, which is an invalid UTF-8 start byte.
        Path p1 = FastaTestResources.path("encoding", "non_utf8_1.txt"); // contains 0xE9 by itself - in ISO-8859-1 this is "é", but as a lone byte it's invalid UTF-8.
        Path p2 = FastaTestResources.path("encoding", "non_utf8_2.txt");

        assertFalse(Utf8Detector.isProbablyUtf8(p));
        assertFalse(Utf8Detector.isProbablyUtf8(p1));
        assertFalse(Utf8Detector.isProbablyUtf8(p2));
    }

    @Test
    void respectsMaxBytesLimit() throws IOException {
        Path p = FastaTestResources.path("encoding", "large_non_utf8.txt"); //1Mb file that is utf8 for the first 1Mb, and non utf8 as the last char
        int limit = 1024 * 1024; // 1 MiB limit to test the file

        assertTrue(Utf8Detector.isProbablyUtf8(p, limit)); //first 1Mb should be detected as Utf8
        assertFalse(Utf8Detector.isProbablyUtf8(p, limit + 1)); //if detecting later characters, should fail
    }
}

