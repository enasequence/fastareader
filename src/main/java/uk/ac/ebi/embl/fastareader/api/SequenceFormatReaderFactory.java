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
import java.nio.file.Files;
import java.nio.file.Path;
import uk.ac.ebi.embl.fastareader.FastaReader;
import uk.ac.ebi.embl.fastareader.SequenceFileFormat;
import uk.ac.ebi.embl.fastareader.SequenceReader;
import uk.ac.ebi.embl.fastareader.api.rereading.SequenceInfoDTO;
import uk.ac.ebi.embl.fastareader.sequenceutils.SequenceAlphabet;

/**
 * Factory methods for opening {@link SequenceReader} instances.
 *
 * <p>This is the single entry point for reading sequence submissions in a unified way,
 * regardless of whether the submission is a multi-record FASTA file or a single-record
 * plain sequence file.</p>
 *
 * <h2>Supported submission types</h2>
 * <ul>
 *   <li><b>FASTA</b>: a FASTA file containing one or more entries. Each entry is addressable by
 *       its fastaReader ID (first entry in the file has the id 0, second one has id 1, etc.).
 *   <li><b>PLAIN_SEQUENCE</b>: a single sequence file containing exactly one sequence record.
 *       The record is addressed by the id 0.
 * </ul>
 *
 * <h2>Resource management</h2>
 * <p>Readers returned by this factory hold file handles. Always close them when done.</p>
 *
 * <pre>{@code
 * try (SequenceReader reader = SequenceReaderFactory.openFasta(fastaFile)) {
 *     // use reader...
 * }
 * }</pre>
 */
public class SequenceFormatReaderFactory {

    private SequenceFormatReaderFactory() {
        // Utility class; do not instantiate.
    }

    /**
     * Opens a reader for a FASTA submission.
     *
     * <p>The FASTA file may contain multiple sequence records. Each record is expected to include
     * a headerline from which the metadata can be parsed. Uses default alphabet which is defined in {@link SequenceAlphabet}.</p>
     *
     * @param fastaFile FASTA file to open (must not be {@code null})
     * @return a {@link SequenceFormatReader} backed by the FASTA file
     * @throws Exception if the file cannot be opened, parsed, or validated (e.g. invalid FASTA,
     *                   invalid headers, duplicate IDs, I/O errors)
     */
    public static SequenceFormatReader readFasta(File fastaFile) throws Exception {
        return new FastaReader(fastaFile, SequenceAlphabet.defaultNucleotideAlphabet());
    }

    /**
     * Opens a reader for a FASTA submission.
     *
     * <p>The FASTA file may contain multiple sequence records. Each record is expected to include
     * a headerline from which the metadata can be parsed, and can be fetched using ids which are assigned according to the order of the fasta entries itself.</p>
     *
     * @param fastaFile FASTA file to open (must not be {@code null})
     * @param alphabet FASTA file to open (must not be {@code null})
     * @return a {@link SequenceReader} backed by the FASTA file
     * @throws Exception if the file cannot be opened, parsed, or validated (e.g. invalid FASTA,
     *                   invalid headers, duplicate IDs, I/O errors)
     */
    public static SequenceFormatReader readFasta(File fastaFile, SequenceAlphabet alphabet) throws Exception {
        return new FastaReader(fastaFile, alphabet);
    }

    /**
     * Opens a reader for a FASTA submission whose content begins at {@code startByteOffset} rather than at
     * the start of the file.
     *
     * <p>This supports partial files that embed FASTA after a leading non-FASTA prefix — for example a GFF3
     * file whose annotation lines are followed by an embedded FASTA section. Bytes before
     * {@code startByteOffset} are ignored while scanning. Uses the default alphabet defined in
     * {@link SequenceAlphabet}.</p>
     *
     * <p>The offset is measured against the reader's <em>decompressed</em> byte stream, so for a
     * BGZF-compressed file it is not a raw compressed-file offset; see
     * {@link uk.ac.ebi.embl.fastareader.FastaReader#FastaReader(File, SequenceAlphabet, long)} for the full
     * coordinate-space and header-boundary contract.
     *
     * @param fastaFile file containing the embedded FASTA content (must not be {@code null})
     * @param startByteOffset absolute offset into the decompressed byte stream to begin scanning from; must be
     *     within {@code [0, decompressedSize]}
     * @return a {@link SequenceFormatReader} backed by the FASTA content beginning at {@code startByteOffset}
     * @throws Exception if the file cannot be opened, parsed, or validated
     */
    public static SequenceFormatReader readFasta(File fastaFile, long startByteOffset) throws Exception {
        return new FastaReader(fastaFile, SequenceAlphabet.defaultNucleotideAlphabet(), startByteOffset);
    }

    /**
     * Opens a reader for a FASTA submission whose content begins at {@code startByteOffset} rather than at
     * the start of the file, using a caller-supplied alphabet.
     *
     * <p>The offset is measured against the reader's <em>decompressed</em> byte stream; see
     * {@link uk.ac.ebi.embl.fastareader.FastaReader#FastaReader(File, SequenceAlphabet, long)} for the full
     * coordinate-space and header-boundary contract.
     *
     * @param fastaFile file containing the embedded FASTA content (must not be {@code null})
     * @param alphabet the {@link SequenceAlphabet} used to interpret the sequence bytes
     * @param startByteOffset absolute offset into the decompressed byte stream to begin scanning from; must be
     *     within {@code [0, decompressedSize]}
     * @return a {@link SequenceFormatReader} backed by the FASTA content beginning at {@code startByteOffset}
     * @throws Exception if the file cannot be opened, parsed, or validated
     */
    public static SequenceFormatReader readFasta(File fastaFile, SequenceAlphabet alphabet, long startByteOffset)
            throws Exception {
        return new FastaReader(fastaFile, alphabet, startByteOffset);
    }

    /**
     * Opens a reader for a plain (single-record) sequence submission.
     *
     * <p>This uses a wrapper {@link SequenceReaderWrapper} which enables interoperability in case both FASTA, plain sequences and/or other formats are ingested simoultaneously.
     * As such, the sequence record is addressed with the id 0 (long). The headerline string is assumed to be non-existent here, so
     * metadata access via {@link SequenceFormatReader#getHeaderline(long)} will return an empty optional value,
     * which can allready be predicted by using the {@link SequenceFormatReader#getSequenceFileFormat()} and judgement on whether the entry contains any metadata.</p>
     *
     * @param sequenceFile plain sequence file to open (must not be {@code null})
     * @return a {@link SequenceReader} backed by the plain sequence file
     * @throws Exception if the file cannot be opened, parsed, or validated (e.g. not UTF-8, wrong format,
     *                   empty file, more than one record, I/O errors)
     */
    public static SequenceFormatReader readPlainSequence(File sequenceFile) throws Exception {
        return new SequenceReaderWrapper(sequenceFile);
    }

    /**
     * Opens a reader for a plain (single-record) sequence submission.
     *
     * <p>This uses a wrapper {@link SequenceReaderWrapper} which enables interoperability in case both FASTA, plain sequences and/or other formats are ingested simoultaneously.
     * As such, the sequence record is addressed with the id 0 (long). The headerline string is assumed to be non-existent here, so
     * metadata access via {@link SequenceFormatReader#getHeaderline(long)} will return an empty optional value,
     * which can allready be predicted by using the {@link SequenceFormatReader#getSequenceFileFormat()} and judgement on whether the entry contains any metadata.</p>
     *
     * @param sequenceFile plain sequence file to open (must not be {@code null})
     * @return a {@link SequenceReader} backed by the plain sequence file
     * @throws Exception if the file cannot be opened, parsed, or validated (e.g. not UTF-8, wrong format,
     *                   empty file, more than one record, I/O errors)
     */
    public static SequenceFormatReader readPlainSequence(File sequenceFile, SequenceAlphabet alphabet)
            throws Exception {
        return new SequenceReaderWrapper(sequenceFile, alphabet);
    }

    /**
     * Opens an appropriate reader from a {@link SequenceInfoDTO}.
     */
    public static SequenceFormatReader readBySequenceInfo(SequenceInfoDTO dto) throws Exception {
        Path filePath = dto.filePath;

        if (!Files.exists(filePath) || !Files.isRegularFile(filePath)) {
            throw new IllegalArgumentException("Not a valid file: " + filePath);
        }

        switch (dto.sequenceFileFormat) {
            case FASTA:
                return new FastaReader(
                        filePath.toFile(),
                        new SequenceAlphabet(dto.sequenceAlphabetSettings),
                        dto.sequenceIndexesMap,
                        dto.headerLines);
            case PLAIN_SEQUENCE:
                if (dto.headerLines != null && !dto.headerLines.isEmpty()) {
                    throw new IllegalArgumentException("Header lines are not supported");
                }
                return new SequenceReaderWrapper(
                        filePath.toFile(), new SequenceAlphabet(dto.sequenceAlphabetSettings), dto.sequenceIndexesMap);
            default:
                throw new IllegalArgumentException("Unrecognized sequence format: " + dto.sequenceFileFormat
                        + ". Allowed values: " + SequenceFileFormat.describe());
        }
    }
}
