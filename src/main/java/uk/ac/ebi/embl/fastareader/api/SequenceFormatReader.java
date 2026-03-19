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
package uk.ac.ebi.embl.fastareader.api;

import java.io.File;
import java.io.Reader;
import java.util.List;
import java.util.Optional;
import uk.ac.ebi.embl.fastareader.SequenceFileFormat;
import uk.ac.ebi.embl.fastareader.SequenceRangeOption;
import uk.ac.ebi.embl.fastareader.SequenceStats;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceIndex;

/**
 * Unified read-only abstraction over submission sequence inputs.
 *
 * <p>Supports two submission types:
 * <ul>
 *   <li>{@link SequenceFileFormat#FASTA}: multi-entry FASTA file </li>
 *   <li>{@link SequenceFileFormat#PLAIN_SEQUENCE}: single-sequence file </li>
 * </ul>
 *
 */
public interface SequenceFormatReader extends AutoCloseable {

    /**
     * Returns the type of input backing this reader.
     */
    SequenceFileFormat getSequenceFileFormat();

    /**
     * Returns the file the reader is reading
     */
    File getFile();

    /**
     * Returns the alphabet the reader is using to read the file
     */
    SequenceAlphabet getSequenceAlphabet();

    /**
     * Returns record ids in a stable, deterministic order.
     *
     * <p>This order must match the order in the underlying file and must remain stable for the lifetime of the reader.
     *
     * <p>For plain sequence:
     * <ul>
     *   <li>returns a single-element list containing the accession id for both id types (or throws for SUBMITTERID,
     *       depending on your chosen semantics).</li>
     * </ul>
     */
    List<Long> getOrderedIds();

    /**
     * Returns optional FASTA header metadata for the record.
     *
     * <p>FASTA:
     * <ul>
     *   <li>Typically present and parsed from the header line.</li>
     *   <li>Missing id should return {@link Optional#empty()} or throw depending on implementation contract.</li>
     * </ul>
     *
     * <p>Plain sequence:
     * <ul>
     *   <li>May be absent because header metadata is not embedded in the sequence file.</li>
     *   <li>If provided externally, this returns it; otherwise {@link Optional#empty()}.</li>
     * </ul>
     *
     * <p><b>Important:</b> The returned Optional itself should never be null.
     */
    Optional<String> getHeaderline(long id);

    /**
     * Returns normalized sequence statistics for the record (lengths, edge Ns, base counts).
     *
     * @throws IllegalArgumentException if the id does not exist or cannot be resolved
     */
    SequenceStats getStats(long id);

    /**
     * Returns a substring slice of the sequence as a String.
     *
     * <p>{@link SequenceRangeOption} controls whether edge N bases are included or excluded when interpreting the range.
     * @param fromBase start base (starts from 1)
     * @param toBase end base (end base should be smaller than total length of the sequence)
     *
     * @throws Exception format-specific exceptions on I/O, invalid ranges, or id resolution failures
     */
    String getSequenceSlice(long id, long fromBase, long toBase, SequenceRangeOption option) throws Exception;

    /**
     * Returns a substring slice of the WHOLE sequence as a String.
     *
     * @param fromBase start base (starts from 1)
     * @param toBase end base (end base should be smaller than total length of the sequence)
     *
     * @throws Exception format-specific exceptions on I/O, invalid ranges, or id resolution failures
     */
    String getSequenceSlice(long id, long fromBase, long toBase) throws Exception;

    /**
     * Returns a streaming {@link Reader} over a slice of the sequence.
     *
     * <p>Useful for large slices to avoid allocating large Strings.
     * The returned reader must be closed by the caller.
     *
     * @throws Exception format-specific exceptions on I/O, invalid ranges, or id resolution failures
     */
    Reader getSequenceSliceReader(long id, long fromBase, long toBase, SequenceRangeOption option) throws Exception;

    /**
     * Returns a streaming {@link Reader} over a slice of the WHOLW sequence.
     *
     * <p>Useful for large slices to avoid allocating large Strings.
     * The returned reader must be closed by the caller.
     *
     * @throws Exception format-specific exceptions on I/O, invalid ranges, or id resolution failures
     */
    Reader getSequenceSliceReader(long id, long fromBase, long toBase) throws Exception;

    /**
     * Gets the sequence index of the sequence entry with the associated Id. Useful for re-reading the file.
     *
     * @throws Exception if the specific index read is not present
     */
    SequenceIndex getSequenceIndex(long id);
}
