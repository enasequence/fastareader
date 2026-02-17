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

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.*;
import lombok.Getter;
import uk.ac.ebi.embl.fastareader.encoding.Utf8Detector;
import uk.ac.ebi.embl.fastareader.exception.SequenceFileException;
import uk.ac.ebi.embl.fastareader.exception.SequenceReadingException;
import uk.ac.ebi.embl.fastareader.sequenceutils.ByteSpan;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceIndex;

/**
 * Reads plain-text biological sequence files efficiently.
 *
 * <p>This reader provides fast, buffered access to a file that contains only sequence characters
 * (e.g. A/C/G/T/N/etc), optionally skipping non-sequence characters such as newlines and carriage returns.
 * The alphabet of allowed bases and skippable allowed characters can be modified with input parameters.
 * Suitable for high-throughput parsing and random-access workflows. </p>
 *
 * <p>Intended for UTF-8 encoded "sequence-only" inputs (no FASTA headers, no quality scores).</p>
 */
public class SequenceReader implements AutoCloseable {
    private static final int UTF_8_CHECK_MAXIMUM_BYTES =
            1024 * 1024; // check just preliminary first 1Mb to confirm encoding is likely UTF8

    /** Returns Sequence entry data */
    @Getter
    private SequenceEntry sequenceEntry = null;

    private SequenceIndex sequenceIndex = null;
    private File file;
    private InternalReader reader;
    private SequenceAlphabet alphabet;

    /** Initializes Sequence reader, skimming through the whole file right away. */
    public SequenceReader(File sequenceFile) throws SequenceFileException, IOException {
        this(sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet());
    }

    /**
     * Initializes Sequence reader, skimming through the whole file right away.
     * Adds the option to define your own desired SequenceAlphabet and a list of tolerable characters in the sequence (usually eg. \n, \r)
     * */
    public SequenceReader(File sequenceFile, SequenceAlphabet alphabet) throws SequenceFileException, IOException {
        this.file = Objects.requireNonNull(sequenceFile, "sequenceFile");
        this.sequenceIndex = null;
        this.sequenceEntry = null;
        this.alphabet = alphabet;
        this.reader = new InternalReader(sequenceFile, this.alphabet, FileFormat.PLAIN_SINGLE_SEQUENCE);

        checkIfUtf8(sequenceFile);
        loadSequence();
    }

    // ---------------------------- queries ----------------------------

    /** Return a sequence slice as a String (no EOLs) for [fromBase..toBase] inclusive. */
    public String getSequenceSliceString(long fromBase, long toBase) throws SequenceFileException {
        return getSequenceSliceString(fromBase, toBase, SequenceRangeOption.WHOLE_SEQUENCE);
    }

    /**
     *  Return a sequence slice as a String (no EOLs) for [fromBase..toBase] inclusive.
     *
     * You can choose whether to read the whole sequence or the interval not including edge N's using the last parameter.
     * */
    public String getSequenceSliceString(long fromBase, long toBase, SequenceRangeOption option)
            throws SequenceFileException {
        ensureFileReaderOpen();

        final ByteSpan span;
        switch (option) {
            case WHOLE_SEQUENCE:
                span = sequenceIndex.byteSpanForBaseRangeIncludingEdgeNBases(fromBase, toBase);
                break;
            case WITHOUT_EDGE_N_BASES:
                span = sequenceIndex.byteSpanForBaseRange(fromBase, toBase);
                break;
            default:
                throw new IllegalStateException("Unknown option " + option);
        }

        try {
            return reader.getSequenceSliceString(span);
        } catch (IOException ioe) {
            throw new SequenceFileException(
                    "I/O while reading slice sequence bytes " + span.start + ".." + (span.endEx - 1), ioe);
        }
    }

    /**
     * Return a sequence slice for reader [fromBase..toBase] (1-based, inclusive).
     * Uses the cached index to translate bases -> bytes, then asks the reader to stream
     * ASCII bytes while skipping '\n' and '\r' on the fly.
     */
    public Reader getSequenceSliceReader(long fromBase, long toBase) {
        return getSequenceSliceReader(fromBase, toBase, SequenceRangeOption.WHOLE_SEQUENCE);
    }

    /**
     * Return a sequence slice for reader [fromBase..toBase] (1-based, inclusive).
     * Uses the cached index to translate bases -> bytes, then asks the reader to stream
     * ASCII bytes while skipping non-base tolerated characters (by default \n, \r) on the fly.
     *
     * You can choose whether to read the whole sequence or the interval not including edge N's using the last parameter.
     */
    public Reader getSequenceSliceReader(long fromBase, long toBase, SequenceRangeOption option) {
        ensureFileReaderOpen();

        ByteSpan span;
        switch (option) {
            case WHOLE_SEQUENCE:
                span = sequenceIndex.byteSpanForBaseRangeIncludingEdgeNBases(fromBase, toBase);
                break;
            case WITHOUT_EDGE_N_BASES:
                span = sequenceIndex.byteSpanForBaseRange(fromBase, toBase);
                break;
            default:
                throw new IllegalStateException("Unknown option " + option);
        }

        return reader.getSequenceSliceReader(span);
    }

    // ---------------------------- interactions with the reader ----------------------------

    public void openNewFile(File sequenceFile) throws SequenceFileException, IOException {
        close(); // if already open, close first
        this.file = Objects.requireNonNull(sequenceFile, "file");
        reader = new InternalReader(sequenceFile, this.alphabet, FileFormat.PLAIN_SINGLE_SEQUENCE);
        checkIfUtf8(sequenceFile);
        loadSequence();
    }

    /** Close the reader. Safe to call multiple times. */
    @Override
    public void close() throws IOException {
        this.sequenceEntry = null;
        this.sequenceIndex = null;
        if (reader != null) {
            reader.close();
            reader = null;
        }
    }

    // ----------------------------- helper methods for actually loading the plain sequence ------------------

    private void checkIfUtf8(File file) throws IOException, SequenceFileException {
        if (!Utf8Detector.isProbablyUtf8(file.toPath(), UTF_8_CHECK_MAXIMUM_BYTES)) {
            throw new SequenceFileException("File is not a UTF-8 compliant file, and as such cannot be processed");
        }
    }

    /**
     * Performs a one-time scan of the sequence file to build in-memory sequence index.
     * The cached index is later used to translate base ranges into byte spans for efficient random-access reads.
     *
     * This method is called once during construction and requires exclusive
     * ownership of the underlying reader.
     */
    private void loadSequence() throws IOException, SequenceFileException {
        List<SequenceEntryMetadata> readEntries;
        try {
            readEntries = reader.readFile();
        } catch (SequenceReadingException e) {
            throw new SequenceFileException(e);
        }

        if (readEntries.isEmpty()) throw new SequenceFileException("No sequence entry found in the file.");
        if (readEntries.size() > 1)
            throw new SequenceFileException("More than one sequence entry found in file, which shouldn't happen.");

        var entry = readEntries.get(0);
        long adjustedBases = entry.getSequenceIndex().totalBases()
                - entry.getSequenceIndex().startNBasesCount
                - entry.getSequenceIndex().endNBasesCount;

        sequenceEntry = new SequenceEntry(
                adjustedBases,
                entry.getSequenceIndex().totalBases(),
                entry.getSequenceIndex().startNBasesCount,
                entry.getSequenceIndex().endNBasesCount,
                entry.getSequenceIndex().caseInsensitiveBaseCount);

        sequenceIndex = entry.getSequenceIndex();
    }

    private void ensureFileReaderOpen() {
        if (reader == null || !reader.readingFile())
            throw new IllegalStateException("Service is not open. Call open() first.");
    }
}
