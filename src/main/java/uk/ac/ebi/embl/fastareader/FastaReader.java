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
import lombok.Setter;
import uk.ac.ebi.embl.fastareader.api.SequenceFormatReader;
import uk.ac.ebi.embl.fastareader.encoding.Utf8Detector;
import uk.ac.ebi.embl.fastareader.exception.FastaFileException;
import uk.ac.ebi.embl.fastareader.exception.SequenceReadingException;
import uk.ac.ebi.embl.fastareader.sequenceutils.ByteSpan;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceIndex;

/**
 * Reads generic FASTA formats with sequences efficiently.
 *
 * <p>This reader provides fast, buffered access to a file that contains header lines,
 * and keeps minimum metadata in memory (headers, base counts), while allowing quick random access to slices
 * or streaming of the sequence. Especially suitable for large sequences, or grabbing headers for further processing. </p>
 *
 * <p>One FASTA entry is considered to be one header with an accompanying sequence, in generic format with the only
 * prerequisite being the greater than '>' character at the start of the header and the next-line '\n' character at
 * the end of the header line. The base alphabet of the sequence is modifiable, as well as the tolerated non-base characters (commonly \n, \r) in the sequence
 * which will be ignored when enumerating bases or streaming the sequence. </p>
 *
 * <p> One file can contain as many FASTA entries as desired while keeping reading and retrieval fast, as
 * long as they are UTF-8 encoded and the bases are within the allowed pre-set alphabet. The entries are enumerated from 0 onwards in order of appearance in the file. </p>
 */
@Getter
@Setter
public final class FastaReader implements AutoCloseable, SequenceFormatReader {

    private static final SequenceFileFormat FILE_FORMAT = SequenceFileFormat.FASTA;

    private static final int UTF_8_CHECK_MAXIMUM_BYTES =
            1024 * 1024; // check just preliminary first 1Mb to confirm encoding is likely UTF8
    private File file;
    private InternalReader reader;
    private SequenceAlphabet alphabet;

    /** FASTA entries, each one with a unique fastaReaderId, header line string and basic sequence information */
    @Getter
    private List<Long> orderedIds = new ArrayList<>();

    private HashMap<Long, String> headerLinesMap = new HashMap<>();
    private HashMap<Long, SequenceIndex> sequenceIndexesMap = new HashMap<>();
    private HashMap<Long, SequenceStats> sequenceStatsMap = new HashMap<>();

    /** Initializes FASTA reader, skimming through the whole file right away. */
    public FastaReader(File fastaFile) throws FastaFileException, IOException {
        this(fastaFile, SequenceAlphabet.defaultNucleotideAlphabet());
    }

    /**
     * Initializes FASTA reader, skimming through the whole file right away.
     * Adds the option to define your own desired SequenceAlphabet and a list of tolerable characters in the sequence (usually eg. \n, \r)
     * */
    public FastaReader(File fastaFile, SequenceAlphabet alphabet) throws FastaFileException, IOException {
        this.file = Objects.requireNonNull(fastaFile, "sequenceFile");
        checkIfUtf8(file);

        resetData();
        this.alphabet = alphabet;
        this.reader = new InternalReader(fastaFile, this.alphabet, FILE_FORMAT);

        loadEntries();
    }

    /**
     * Initializes FASTA reader, loading the sequenceIndexes.
     * Adds the option to define your own desired SequenceAlphabet and a list of tolerable characters in the sequence (usually eg. \n, \r)
     * */
    public FastaReader(File fastaFile, HashMap<Long, SequenceIndex> indexes, HashMap<Long, String> headerLines)
            throws FastaFileException, IOException {
        this.file = Objects.requireNonNull(fastaFile, "sequenceFile");
        checkIfUtf8(file);

        resetData();
        this.sequenceIndexesMap = indexes;
        setUpSequenceStatsFromIndexes();
        this.headerLinesMap = headerLines==null?new HashMap<>():headerLines;


        this.reader = new InternalReader(fastaFile, this.alphabet, FILE_FORMAT);

        loadEntries();
    }

    @Override
    public SequenceFileFormat getSequenceFileFormat() {
        return FILE_FORMAT;
    }

    @Override
    public File getFile() {
        return file;
    }

    @Override
    public SequenceAlphabet getSequenceAlphabet() {
        return alphabet;
    }

    // ---------------------------- queries ----------------------------

    @Override
    public Optional<String> getHeaderline(long id) {
        return Optional.ofNullable(headerLinesMap.get(id));
    }

    @Override
    public SequenceStats getStats(long id) {
        return sequenceStatsMap.get(id);
    }

    /**
     *  Return a sequence slice as a String (no EOLs) for [fromBase..toBase] inclusive.
     *
     * You can choose whether to read the whole sequence or the interval not including edge N's using the last parameter.
     * */
    @Override
    public String getSequenceSlice(long id, long fromBase, long toBase, SequenceRangeOption option) throws Exception {
        ensureFileReaderOpen();
        SequenceIndex index = sequenceIndexesMap.get(id);
        if (index == null) {
            throw new FastaFileException("No sequence index found for fasta reader Id " + id);
        }

        final ByteSpan span = getByteSpan(fromBase, toBase, option, index);
        try {
            return reader.getSequenceSliceString(span);
        } catch (IOException ioe) {
            throw new FastaFileException(
                    "I/O while reading slice for " + id + " bytes " + span.start + ".." + (span.endEx - 1), ioe);
        }
    }

    /** Return a sequence slice as a String (no EOLs) for [fromBase..toBase] inclusive. */
    @Override
    public String getSequenceSlice(long id, long fromBase, long toBase) throws Exception {
        return getSequenceSlice(id, fromBase, toBase, SequenceRangeOption.WHOLE_SEQUENCE);
    }

    /**
     * Return a sequence slice for reader [fromBase..toBase] (1-based, inclusive) for the given ID.
     * Uses the cached index to translate bases -> bytes, then asks the reader to stream
     * ASCII bytes while skipping '\n' and '\r' on the fly.
     *
     * You can choose whether to read the whole sequence or the interval not including edge N's using the last parameter.
     */
    @Override
    public Reader getSequenceSliceReader(long id, long fromBase, long toBase, SequenceRangeOption option)
            throws Exception {
        ensureFileReaderOpen();
        var index = sequenceIndexesMap.get(id);
        if (index == null) {
            throw new FastaFileException("No sequence index found for fasta reader Id " + id);
        }

        ByteSpan span = getByteSpan(fromBase, toBase, option, index);
        return reader.getSequenceSliceReader(span);
    }

    /**
     * Return a sequence slice for reader [fromBase..toBase] (1-based, inclusive) for the given ID.
     * Uses the cached index to translate bases -> bytes, then asks the reader to stream
     * ASCII bytes while skipping '\n' and '\r' on the fly.
     */
    @Override
    public Reader getSequenceSliceReader(long id, long fromBase, long toBase) throws Exception {
        return getSequenceSliceReader(id, fromBase, toBase, SequenceRangeOption.WHOLE_SEQUENCE);
    }

    @Override
    public SequenceIndex getSequenceIndex(long id) {
        return sequenceIndexesMap.get(id);
    }

    // ---------------------------- interactions with the reader ----------------------------

    public void openNewFile(File fastaFile) throws FastaFileException, IOException {
        close(); // if already open, close first
        this.file = Objects.requireNonNull(fastaFile, "file");
        reader = new InternalReader(fastaFile, this.alphabet, SequenceFileFormat.FASTA);
        checkIfUtf8(fastaFile);
        loadEntries();
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

    // ----------------------------- helper methods for actually loading the fastaEntries ------------------

    private void checkIfUtf8(File file) throws IOException, FastaFileException {
        if (!Utf8Detector.isProbablyUtf8(file.toPath(), UTF_8_CHECK_MAXIMUM_BYTES)) {
            throw new FastaFileException("File is not a UTF-8 compliant file, and as such cannot be processed");
        }
    }

    private void resetData() {
        this.orderedIds = new ArrayList<>();
        this.headerLinesMap = new HashMap<>();
        this.sequenceStatsMap = new HashMap<>();
        this.sequenceIndexesMap = new HashMap<>();
    }

    /**
     * Performs a one-time scan of the FASTA file to build in-memory sequence indexes.
     * The cached indexes are later used to translate base ranges into byte spans for efficient random-access reads.
     *
     * This method is called once during construction and requires exclusive
     * ownership of the underlying reader.
     */
    private void loadEntries() throws IOException, FastaFileException {
        List<SequenceEntryMetadata> readEntries;
        try {
            readEntries = reader.readFile();
        } catch (SequenceReadingException e) {
            throw new FastaFileException(e);
        }

        long currentId = 0;
        for (var entry : readEntries) {
            orderedIds.add(currentId);
            headerLinesMap.put(currentId, entry.getHeaderLine());
            sequenceIndexesMap.put(currentId, entry.getSequenceIndex());
            currentId++;
        }

        setUpSequenceStatsFromIndexes();
    }

    private void setUpSequenceStatsFromIndexes() {
        for (var fastaReaderId : sequenceIndexesMap.keySet()) {
            var sequenceIndex = sequenceIndexesMap.get(fastaReaderId);
            long adjustedBases =
                    sequenceIndex.totalBases() - sequenceIndex.startNBasesCount - sequenceIndex.endNBasesCount;
            var sequenceStats = new SequenceStats(
                    sequenceIndex.totalBases(),
                    adjustedBases,
                    sequenceIndex.startNBasesCount,
                    sequenceIndex.endNBasesCount,
                    sequenceIndex.caseInsensitiveBaseCount);
            sequenceStatsMap.put(fastaReaderId, sequenceStats);
        }
    }

    private void ensureFileReaderOpen() {
        if (reader == null || !reader.readingFile())
            throw new IllegalStateException("Service is not open. Call open() first.");
    }

    private static ByteSpan getByteSpan(long fromBase, long toBase, SequenceRangeOption option, SequenceIndex index) {
        ByteSpan span;
        switch (option) {
            case WHOLE_SEQUENCE:
                span = index.byteSpanForBaseRangeIncludingEdgeNBases(fromBase, toBase);
                break;
            case WITHOUT_EDGE_N_BASES:
                span = index.byteSpanForBaseRange(fromBase, toBase);
                break;
            default:
                throw new IllegalStateException("Unknown option " + option);
        }
        return span;
    }
}
