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
import uk.ac.ebi.embl.fastareader.encoding.Utf8Detector;
import uk.ac.ebi.embl.fastareader.exception.FastaFileException;
import uk.ac.ebi.embl.fastareader.exception.SequenceFileException;
import uk.ac.ebi.embl.fastareader.sequenceutils.ByteSpan;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceIndex;

public class SequenceReader implements AutoCloseable {
    private int UTF_8_CHECK_MAXIMUM_BYTES =
            1024 * 1024; // check just preliminary first 1Mb to confirm encoding is likely UTF8

    public SequenceEntry sequenceEntry = null;
    private SequenceIndex sequenceIndex = null;
    private File file;
    private SequentialFileReader reader;

    /** Initializes Sequence reader, skimming through the whole file right away. */
    public SequenceReader(File sequenceFile) throws FastaFileException, IOException {
        this(sequenceFile, SequenceAlphabet.defaultNucleotideAlphabet());
    }

    /**
     * Initializes Sequence reader, skimming through the whole file right away.
     * Adds the option to define your own desired SequenceAlphabet and a list of tolerable characters in the sequence (usually eg. \n, \r)
     * */
    public SequenceReader(File sequenceFile, SequenceAlphabet alphabet) throws SequenceFileException, IOException {
        this.file = Objects.requireNonNull(sequenceFile, "sequenceFile");
        this.reader = new SequentialFileReader(sequenceFile, alphabet, FileFormat.SINGLE_SEQUENCE);
        this.sequenceIndex = null;
        this.sequenceEntry = null;

        try {
            checkIfUtf8(sequenceFile);
            loadSequence();
        } catch (FastaFileException e) {
            throw new SequenceFileException(e); // small wrap to avoid fasta file exception confusion
        }
    }

    // ---------------------------- queries ----------------------------

    /** Returns Sequence entry data */
    public SequenceEntry getSequenceInfo() {
        return sequenceEntry;
    }

    /** Return a sequence slice as a String (no EOLs) for [fromBase..toBase] inclusive. */
    public String getSequenceSliceString(long fromBase, long toBase) throws FastaFileException {
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
     * Return a sequence slice for reader [fromBase..toBase] (1-based, inclusive) for the given ID.
     * Uses the cached index to translate bases -> bytes, then asks the reader to stream
     * ASCII bytes while skipping '\n' and '\r' on the fly.
     */
    public Reader getSequenceSliceReader(long fromBase, long toBase) throws FastaFileException {
        return getSequenceSliceReader(fromBase, toBase, SequenceRangeOption.WHOLE_SEQUENCE);
    }

    /**
     * Return a sequence slice for reader [fromBase..toBase] (1-based, inclusive) for the given ID.
     * Uses the cached index to translate bases -> bytes, then asks the reader to stream
     * ASCII bytes while skipping '\n' and '\r' on the fly.
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

    public void openNewFile(File fastaFile) throws FastaFileException, IOException {
        close(); // if already open, close first
        this.file = Objects.requireNonNull(fastaFile, "file");
        reader = new SequentialFileReader(
                fastaFile, SequenceAlphabet.defaultNucleotideAlphabet(), FileFormat.SINGLE_SEQUENCE);
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

    // ----------------------------- helper methods for actually loading the fastaEntries ------------------

    private void checkIfUtf8(File file) throws IOException, FastaFileException {
        if (!Utf8Detector.isProbablyUtf8(file.toPath(), UTF_8_CHECK_MAXIMUM_BYTES)) {
            throw new FastaFileException("File is not a UTF-8 compliant file, and as such cannot be processed");
        }
    }

    /**
     * Performs a one-time scan of the sequence file to build in-memory sequence indexes.
     * The cached indexes are later used to translate base ranges into byte spans for efficient random-access reads.
     *
     * This method is called once during construction and requires exclusive
     * ownership of the underlying reader.
     */
    private void loadSequence() throws IOException, FastaFileException {
        List<SequenceEntryMetadata> readEntries = reader.readFile();

        if (readEntries.isEmpty()) throw new SequenceFileException("No sequence entry found in the file.");
        if (readEntries.size() > 1)
            throw new SequenceFileException("More than one sequence entry found in file, which shouldn't happen.");

        var entry = readEntries.get(0);
        long adjustedBases = entry.sequenceIndex.totalBases()
                - entry.sequenceIndex.startNBasesCount
                - entry.sequenceIndex.endNBasesCount;

        sequenceEntry = new SequenceEntry();
        sequenceEntry.setTotalBases(entry.sequenceIndex.totalBases());
        sequenceEntry.setLeadingNsCount(entry.sequenceIndex.startNBasesCount);
        sequenceEntry.setTrailingNsCount(entry.sequenceIndex.endNBasesCount);
        sequenceEntry.setBaseCount(entry.sequenceIndex.caseInsensitiveBaseCount);
        sequenceEntry.setBaseCount(entry.sequenceIndex.caseInsensitiveBaseCount);
        sequenceEntry.setTotalBasesWithoutNBases(adjustedBases);

        sequenceIndex = entry.sequenceIndex;
    }

    private void ensureFileReaderOpen() {
        if (reader == null || !reader.readingFile())
            throw new IllegalStateException("Service is not open. Call open() first.");
    }
}
