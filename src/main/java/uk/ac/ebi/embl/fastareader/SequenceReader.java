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

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.*;
import lombok.Getter;
import uk.ac.ebi.embl.fastareader.encoding.Utf8Detector;
import uk.ac.ebi.embl.fastareader.exception.SequenceFileException;
import uk.ac.ebi.embl.fastareader.exception.SequenceReadingException;
import uk.ac.ebi.embl.fastareader.io.BgzfByteReader;
import uk.ac.ebi.embl.fastareader.io.BgzfDetector;
import uk.ac.ebi.embl.fastareader.io.FileChannelByteReader;
import uk.ac.ebi.embl.fastareader.io.SeekableByteReader;
import uk.ac.ebi.embl.fastareader.sequenceutils.ByteSpan;
import uk.ac.ebi.embl.fastareader.sequenceutils.GapRegion;
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

    private static final SequenceFileFormat FILE_FORMAT = SequenceFileFormat.PLAIN_SEQUENCE;

    private static final int UTF_8_CHECK_MAXIMUM_BYTES =
            1024 * 1024; // check just preliminary first 1Mb to confirm encoding is likely UTF8

    private InternalReader reader;

    /** File & the accompanying alphabet it's read with */
    @Getter
    private File file;

    @Getter
    private SequenceAlphabet sequenceAlphabet;
    /** Returns Sequence entry data */
    @Getter
    private SequenceStats stats;

    @Getter
    private SequenceIndex sequenceIndex;

    /** Initializes Sequence reader, skimming through the whole file right away. */
    public SequenceReader(File sequenceFile) throws SequenceFileException, IOException {
        this.file = Objects.requireNonNull(sequenceFile, "sequenceFile");
        resetData();
        this.sequenceAlphabet = SequenceAlphabet.defaultNucleotideAlphabet();
        initAndLoad(sequenceFile);
    }

    /**
     * Initializes Sequence reader, skimming through the whole file right away.
     * Adds the option to define your own desired SequenceAlphabet and a list of tolerable characters in the sequence (usually eg. \n, \r)
     * */
    public SequenceReader(File sequenceFile, SequenceAlphabet sequenceAlphabet)
            throws SequenceFileException, IOException {
        this.file = Objects.requireNonNull(sequenceFile, "sequenceFile");
        resetData();
        this.sequenceAlphabet = sequenceAlphabet;
        initAndLoad(sequenceFile);
    }

    private void initAndLoad(File file) throws SequenceFileException, IOException {
        this.reader = new InternalReader(file, this.sequenceAlphabet, FILE_FORMAT, openAndValidateUtf8(file));
        try {
            loadSequence();
        } catch (SequenceFileException | IOException e) {
            try { reader.close(); } catch (IOException ignored) {}
            reader = null;
            throw e;
        }
    }

    /**
     * Initializes Sequence reader, loading the sequenceIndex that was previously read with the provided alphabet from the since unchanged file.
     * */
    public SequenceReader(File file, SequenceAlphabet sequenceAlphabet, SequenceIndex sequenceIndex)
            throws SequenceFileException, IOException {

        this.file = Objects.requireNonNull(file, "sequenceFile");

        resetData();
        this.sequenceAlphabet = sequenceAlphabet;
        this.reader = new InternalReader(file, this.sequenceAlphabet, FILE_FORMAT, openAndValidateUtf8(file));

        this.sequenceIndex = sequenceIndex;
        setUpSequenceStatsFromIndex();
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
        final ByteSpan span = getByteSpan(fromBase, toBase, option);

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
        ByteSpan span = getByteSpan(fromBase, toBase, option);
        return reader.getSequenceSliceReader(span);
    }

    public List<GapRegion> getGapRegions() {
        return getGapRegions(SequenceRangeOption.WHOLE_SEQUENCE);
    }

    public List<GapRegion> getGapRegions(SequenceRangeOption option) {
        return sequenceIndex.gapRegionsView(option);
    }

    public List<GapRegion> getGapRegions(long fromBase, long toBase) {
        return getGapRegions(fromBase, toBase, SequenceRangeOption.WHOLE_SEQUENCE);
    }

    public List<GapRegion> getGapRegions(long fromBase, long toBase, SequenceRangeOption option) {
        return sequenceIndex.gapRegionsView(fromBase, toBase, option);
    }

    // ---------------------------- interactions with the reader ----------------------------

    public SequenceFileFormat getSequenceFileFormat() {
        return FILE_FORMAT;
    }

    public void openNewFile(File sequenceFile) throws SequenceFileException, IOException {
        close(); // if already open, close first
        this.file = Objects.requireNonNull(sequenceFile, "file");
        this.reader =
                new InternalReader(sequenceFile, this.sequenceAlphabet, FILE_FORMAT, openAndValidateUtf8(sequenceFile));
        loadSequence();
    }

    /** Close the reader. Safe to call multiple times. */
    @Override
    public void close() throws IOException {
        resetData();
        if (reader != null) {
            reader.close();
            reader = null;
        }
    }

    // ----------------------------- helper methods for actually loading the plain sequence ------------------

    /**
     * Detects compression, rejects plain gzip, validates UTF-8, and returns an open
     * {@link SeekableByteReader} over the logical byte stream. For BGZF files this reuses the
     * same reader instance for validation and parsing, avoiding a second block-directory walk.
     */
    private static SeekableByteReader openAndValidateUtf8(File file) throws IOException, SequenceFileException {
        BgzfDetector.Compression comp = BgzfDetector.detect(file);
        if (comp == BgzfDetector.Compression.PLAIN_GZIP) {
            throw new SequenceFileException(
                    "Plain gzip is not supported; recompress with bgzip: " + file.getAbsolutePath());
        }
        if (comp == BgzfDetector.Compression.BGZF) {
            BgzfByteReader reader = new BgzfByteReader(file);
            boolean ok;
            try {
                ok = Utf8Detector.isProbablyUtf8(
                        new ByteArrayInputStream(BgzfReaderUtils.decompressedPrefix(reader, UTF_8_CHECK_MAXIMUM_BYTES)),
                        UTF_8_CHECK_MAXIMUM_BYTES);
            } catch (IOException e) {
                reader.close();
                throw e;
            }
            if (!ok) {
                reader.close();
                throw new SequenceFileException("File is not a UTF-8 compliant file, and as such cannot be processed");
            }
            return reader;
        }
        boolean ok = Utf8Detector.isProbablyUtf8(file.toPath(), UTF_8_CHECK_MAXIMUM_BYTES);
        if (!ok) {
            throw new SequenceFileException("File is not a UTF-8 compliant file, and as such cannot be processed");
        }
        return new FileChannelByteReader(file);
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

        sequenceIndex = entry.getSequenceIndex();
        setUpSequenceStatsFromIndex();
    }

    private void setUpSequenceStatsFromIndex() {
        long adjustedBases = sequenceIndex.totalBases() - sequenceIndex.startNBasesCount - sequenceIndex.endNBasesCount;
        stats = new SequenceStats(
                sequenceIndex.totalBases(),
                adjustedBases,
                sequenceIndex.startNBasesCount,
                sequenceIndex.endNBasesCount,
                sequenceIndex.caseInsensitiveBaseCount);
    }

    private void resetData() {
        this.sequenceIndex = null;
        this.stats = null;
    }

    private void ensureFileReaderOpen() {
        if (reader == null || !reader.readingFile())
            throw new IllegalStateException("Service is not open. Call open() first.");
    }

    private ByteSpan getByteSpan(long fromBase, long toBase, SequenceRangeOption option) {
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
        return span;
    }
}
