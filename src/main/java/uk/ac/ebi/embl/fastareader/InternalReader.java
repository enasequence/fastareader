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

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.util.*;
import uk.ac.ebi.embl.fastareader.exception.SequenceReadingException;
import uk.ac.ebi.embl.fastareader.headerutils.HeaderLineDecoder;
import uk.ac.ebi.embl.fastareader.io.SeekableByteReader;
import uk.ac.ebi.embl.fastareader.io.SeekableByteReaderFactory;
import uk.ac.ebi.embl.fastareader.sequenceutils.*;

class InternalReader implements AutoCloseable {

    private static final int BUFFER_SIZE = 4 * 1024 * 1024; // 4 MB
    private static final int CHAR_BUF_SIZE = 512 * 1024; // 512 KB
    private static final byte GT = (byte) '>';
    private static final byte LF = (byte) '\n';
    private static final byte CR = (byte) '\r';

    private final SeekableByteReader reader;
    private final long fileSize;
    private final SequenceAlphabet alphabet;
    private final SequenceFileFormat fileFormat;
    HeaderLineDecoder headerLineDecoder;

    InternalReader(File file, SequenceFileFormat fileFormat) throws IOException {
        this(file, SequenceAlphabet.defaultNucleotideAlphabet(), fileFormat);
    }

    InternalReader(File file, SequenceAlphabet alphabet, SequenceFileFormat fileFormat) throws IOException {
        this(file, alphabet, fileFormat, SeekableByteReaderFactory.open(file));
    }

    InternalReader(File file, SequenceAlphabet alphabet, SequenceFileFormat fileFormat, SeekableByteReader reader)
            throws IOException {
        Objects.requireNonNull(file, "Input file is null");
        if (!file.exists()) throw new FileNotFoundException(file.getAbsolutePath());
        if (file.isDirectory()) throw new FileNotFoundException("Directory: " + file.getAbsolutePath());
        if (!file.canRead()) throw new IllegalArgumentException("No read permission: " + file.getAbsolutePath());
        this.alphabet = Objects.requireNonNull(alphabet, "alphabet");
        this.fileFormat = Objects.requireNonNull(fileFormat, "fileFormat");
        this.reader = Objects.requireNonNull(reader, "reader");
        this.fileSize = reader.size();
        this.headerLineDecoder =
                new HeaderLineDecoder(StandardCharsets.UTF_8, new HashSet<>(Arrays.asList(LF, CR)), BUFFER_SIZE);
    }

    @Override
    public void close() throws IOException {
        reader.close();
    }

    boolean readingFile() {
        return reader.isOpen();
    }

    String getSequenceSliceString(ByteSpan span) throws IOException {
        long byteStart = span.start;
        long byteEndExclusive = span.endEx;

        if (byteStart < 0 || byteEndExclusive < byteStart || byteEndExclusive > fileSize) {
            throw new IllegalArgumentException("Bad byte window: " + byteStart + ".." + byteEndExclusive);
        }
        long remain = byteEndExclusive - byteStart;
        long off = byteStart;

        // pre-size builder with a sane cap (skip newlines, so content <= remain)
        int expect = (int) Math.min(remain, 1_000_000L);
        StringBuilder sb = new StringBuilder(expect);

        ByteBuffer buf = ByteBuffer.allocateDirect(BUFFER_SIZE);
        while (remain > 0) {
            buf.clear();
            int want = (int) Math.min(buf.capacity(), remain);
            buf.limit(want);
            int n = reader.read(buf, off);
            if (n <= 0) break;
            buf.flip();
            while (buf.hasRemaining()) {
                byte b = buf.get();
                if (alphabet.isNonSequenceAllowedChar(b)) continue;
                sb.append((char) (b & 0xFF));
            }
            remain -= n;
            off += n;
        }
        return sb.toString().toUpperCase();
    }

    /** Char-stream view over [span.start, span.endEx): ASCII decode, skip LF/CR.
     *  Uses absolute reads; does NOT change reader.position(). */
    java.io.Reader getSequenceSliceReader(ByteSpan span) {
        final long start = span.start;
        final long endEx = span.endEx;

        return new java.io.Reader() {
            private long pos = start;

            private final java.nio.ByteBuffer buf = java.nio.ByteBuffer.allocateDirect(CHAR_BUF_SIZE);

            {
                // allocate buffer and mark it EMPTY so the very first read() refills it from the reader.
                // Without this, hasRemaining() is true and we'll read uninitialized bytes (→ '\0').
                buf.limit(0);
            }

            @Override
            public int read(
                    char[] characterBuffer, int startingWriteIndexInCharacterBuffer, int maximumNumberOfCharsToRead)
                    throws java.io.IOException {
                // --- Validate caller's target window [off .. off + len) ---
                validateTargetWindow(characterBuffer, startingWriteIndexInCharacterBuffer, maximumNumberOfCharsToRead);
                if (maximumNumberOfCharsToRead == 0) return 0;

                int out = 0;
                while (out < maximumNumberOfCharsToRead) {
                    // --- Prep the buffer for next read & fill it out ---
                    if (!buf.hasRemaining()) {
                        if (pos >= endEx) break; // if end of slice reached, stop reading

                        buf.clear();
                        int toRead = (int) Math.min(buf.capacity(), endEx - pos);
                        buf.limit(toRead);

                        int n = reader.read(buf, pos);
                        if (n <= 0) break; // if no bytes were read, break
                        pos += n;
                        buf.flip();
                    }
                    // Drain bytes + ASCII decode -> writees chars into caller's window [off .. off+len)
                    while (buf.hasRemaining() && out < maximumNumberOfCharsToRead) {
                        byte b = buf.get();
                        if (alphabet.isNonSequenceAllowedChar(b)) continue; // skip irrelevant bytes
                        char c = Character.toUpperCase((char) (b & 0xFF)); // decode byte and shift it to uppercase
                        characterBuffer[startingWriteIndexInCharacterBuffer + out] = c;

                        out++;
                    }
                }
                // If we produced nothing AND we're at EOF, signal -1
                return (out == 0) ? -1 : out;
            }

            private void validateTargetWindow(char[] buffer, int offset, int length) {

                Objects.requireNonNull(buffer, "buffer");

                if (offset < 0 || length < 0) {
                    throw new IndexOutOfBoundsException(
                            "offset and length must be non-negative (offset=" + offset + ", length=" + length + ')');
                }

                if (offset + length > buffer.length) {
                    throw new IndexOutOfBoundsException("Requested range exceeds buffer size "
                            + "(offset=" + offset
                            + ", length=" + length
                            + ", bufferLength=" + buffer.length + ')');
                }
            }

            @Override
            public int read() throws java.io.IOException {
                char[] one = new char[1];
                int n = read(one, 0, 1);
                return (n == -1) ? -1 : one[0];
            }

            @Override
            public boolean ready() {
                return buf.hasRemaining() || pos < endEx;
            }

            @Override
            public void close() {
                /* no-op, reader is kept alive */
            }
        };
    }

    List<SequenceEntryMetadata> readFile() throws SequenceReadingException, IOException {
        long position = 0;
        List<SequenceEntryMetadata> entries = new ArrayList<>();

        switch (fileFormat) {
            case FASTA:
                while (true) {
                    Optional<SequenceEntryMetadata> fastaEntryMetadata = readNextFastaEntry(position);
                    if (fastaEntryMetadata.isEmpty()) break;
                    entries.add(fastaEntryMetadata.get());
                    position = fastaEntryMetadata
                            .get()
                            .getSequenceIndex()
                            .lastBaseByte; // read from the end of last sequence
                }
                break;
            case PLAIN_SEQUENCE:
                Optional<SequenceEntryMetadata> sequenceEntryMetadata = readFileAsSequence();
                if (sequenceEntryMetadata.isEmpty()) break;
                entries.add(sequenceEntryMetadata.get());
                break;
            default:
                throw new SequenceReadingException("Unsupported file format: " + fileFormat);
        }
        return entries;
    }

    /** Reads the next FASTA entry starting at or after 'from'. */
    private Optional<SequenceEntryMetadata> readNextFastaEntry(long from) throws SequenceReadingException {
        try {
            OptionalLong headerPosOpt = seekToNextHeader(from);
            if (headerPosOpt.isEmpty()) return Optional.empty();

            long headerPos = headerPosOpt.getAsLong();
            String headerLine = headerLineDecoder.readHeaderLine(reader, headerPos);
            if (headerLine == null)
                throw new SequenceReadingException("FASTA header is malformed at byte " + headerPos);

            long sequenceStartPos = reader.position(); // first byte after header line is the sequence position
            SequenceIndexBuilder sib = new SequenceIndexBuilder(reader, alphabet, Optional.of(GT));
            SequenceIndex index = sib.buildFrom(sequenceStartPos);

            // Move reader cursor to the sequence start position
            reader.position(sequenceStartPos);

            return Optional.of(new SequenceEntryMetadata(headerLine, headerPos, index));
        } catch (IOException io) {
            long pos = safePos();
            throw new SequenceReadingException("I/O while reading FASTA at byte " + pos + ": " + io.getMessage(), io);
        }
    }

    /** Reads the single sequence from start to end */
    private Optional<SequenceEntryMetadata> readFileAsSequence() throws SequenceReadingException {
        try {
            long sequenceStartPos = reader.position(); // first byte should be sequence start
            SequenceIndexBuilder sib = new SequenceIndexBuilder(reader, alphabet, Optional.empty());
            SequenceIndex index = sib.buildFrom(sequenceStartPos);

            // Move reader cursor to the sequence start position
            reader.position(sequenceStartPos);

            SequenceEntryMetadata e = new SequenceEntryMetadata(index);

            if (e.sequenceIndex.totalBases() == 0) return Optional.empty();
            return Optional.of(e);
        } catch (IOException io) {
            long pos = safePos();
            throw new SequenceReadingException(
                    "I/O while reading sequence at byte " + pos + ": " + io.getMessage(), io);
        }
    }

    // ------------------ header seeking & line reading ------------------

    private OptionalLong seekToNextHeader(long from) throws IOException {
        if (from >= fileSize) return OptionalLong.empty();
        ByteBuffer buf = ByteBuffer.allocateDirect(BUFFER_SIZE);

        while (from < fileSize) {
            buf.clear();
            int want = (int) Math.min(buf.capacity(), fileSize - from);
            buf.limit(want);
            int n = reader.read(buf, from);
            if (n <= 0) break;
            buf.flip();
            while (buf.hasRemaining()) {
                long abs = from + buf.position();
                if (buf.get() == GT && isLineStart(abs)) {
                    return OptionalLong.of(abs);
                }
            }
            from += n;
        }
        return OptionalLong.empty();
    }

    private boolean isLineStart(long abs) throws IOException {
        if (abs == 0) return true;
        if (abs > fileSize) return false;
        return peek(abs - 1) == LF;
    }

    private byte peek(long abs) throws IOException {
        if (abs < 0 || abs >= fileSize) return 0;
        ByteBuffer one = ByteBuffer.allocate(1);
        int n = reader.read(one, abs);
        return (n == 1) ? one.get(0) : 0;
    }

    private long safePos() {
        try {
            return reader.position();
        } catch (IOException e) {
            return -1;
        }
    }
}
