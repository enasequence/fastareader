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

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.charset.*;
import java.util.Set;
import uk.ac.ebi.embl.fastareader.io.SeekableByteReader;

public class HeaderLineDecoder {
    private Set<Byte> endChars;
    private int BUFFER_SIZE;
    private CharsetDecoder decoder;

    public HeaderLineDecoder(Charset charset, Set<Byte> endChars, int BUFFER_SIZE) {
        this.endChars = endChars;
        this.BUFFER_SIZE = BUFFER_SIZE;
        this.decoder = charset.newDecoder()
                .onMalformedInput(CodingErrorAction.REPORT)
                .onUnmappableCharacter(CodingErrorAction.REPORT);
    }

    /** Reads one Unicode line from input position, assuming the position handed to it contains '>', advances past LF/CR or to EOF. */
    public String readHeaderLine(SeekableByteReader reader, long from) throws IOException {
        long fileSize = reader.size();
        if (from >= fileSize) return null;

        reader.position(from);
        long scanPos = reader.position();

        decoder.reset();
        StringBuilder sb = new StringBuilder(512);
        ByteBuffer buffer = ByteBuffer.allocateDirect(BUFFER_SIZE);
        CharBuffer charBuf = CharBuffer.allocate(BUFFER_SIZE);

        while (scanPos < fileSize) {
            buffer.clear();
            int want = (int) Math.min(buffer.capacity(), fileSize - scanPos);
            buffer.limit(want);

            int n = reader.read(buffer, scanPos);
            if (n <= 0) break; // no characters to read, eg. hit end of file

            buffer.flip();

            int lineEndIndex = indexOfHeaderLineEnd(buffer);
            if (lineEndIndex >= 0) {
                // found headerLine end, decode remaining bytes up to (but not including) LF/CR. & reposition reader
                decodeBytes(sb, charBuf, buffer, lineEndIndex, true);
                long nextLineStart = scanPos + lineEndIndex + 1;
                reader.position(nextLineStart);
                return sb.toString();

            } else {
                // Decode entire buffer
                decodeBytes(sb, charBuf, buffer, buffer.remaining(), false);
                scanPos += n;
            }
        }

        // reached EOF: finalize decoder to flush any buffered bytes -> chars.
        finishDecoding(sb, charBuf);

        reader.position(fileSize);
        return sb.toString();
    }

    private int indexOfHeaderLineEnd(ByteBuffer buf) {
        for (int i = 0; i < buf.remaining(); i++) {
            byte b = buf.get(buf.position() + i);
            if (endChars.contains(b)) return i;
        }
        return -1;
    }

    /**
     * Decodes exactly 'len' bytes from 'buf' (starting at current position) into sb.
     * If endOfInput is true, we finalize decoder state for this line (but do not flush decoder globally).
     */
    private void decodeBytes(StringBuilder sb, CharBuffer charBuf, ByteBuffer buf, int len, boolean endOfInput)
            throws CharacterCodingException {
        ByteBuffer slice = buf.slice();
        slice.limit(len);
        buf.position(buf.position() + len);

        while (true) {
            CoderResult r = decoder.decode(slice, charBuf, endOfInput);
            charBuf.flip();
            sb.append(charBuf);
            charBuf.clear();

            if (r.isUnderflow()) break; // need more input or we're done
            if (r.isOverflow()) continue; // charBuf filled; loop to drain
            r.throwException(); // malformed/unmappable
        }

        if (endOfInput) {
            flushDecoder(sb, charBuf);
        }
    }

    /** Flushes decoder at EOF **/
    private void finishDecoding(StringBuilder sb, CharBuffer charBuf) throws CharacterCodingException {
        // Signal end of input with an empty buffer to emit any buffered state.
        decoder.decode(ByteBuffer.allocate(0), charBuf, true);
        charBuf.flip();
        sb.append(charBuf);
        charBuf.clear();

        flushDecoder(sb, charBuf);
    }

    /** Flushes decoder **/
    private void flushDecoder(StringBuilder sb, CharBuffer charBuf) throws CharacterCodingException {
        while (true) {
            CoderResult r = decoder.flush(charBuf);
            charBuf.flip();
            sb.append(charBuf);
            charBuf.clear();

            if (r.isUnderflow()) break;
            if (r.isOverflow()) continue;
            r.throwException();
        }
    }
}
