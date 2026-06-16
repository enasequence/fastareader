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
package uk.ac.ebi.embl.fastareader.io;

import java.io.File;
import java.io.IOException;

/**
 * Chooses the appropriate {@link SeekableByteReader} implementation for a file based on its
 * detected compression: plain files use {@link FileChannelByteReader}; BGZF files use
 * {@link BgzfByteReader}; plain (non-BGZF) gzip is rejected because it is not seekable.
 */
public final class SeekableByteReaderFactory {

    private SeekableByteReaderFactory() {}

    /** Opens a {@link SeekableByteReader} over the logical (uncompressed) stream of {@code file}. */
    public static SeekableByteReader open(File file) throws IOException {
        switch (BgzfDetector.detect(file)) {
            case BGZF:
                return new BgzfByteReader(file);
            case PLAIN_GZIP:
                throw new IOException("Plain gzip is not seekable; recompress with bgzip: " + file.getAbsolutePath());
            case UNCOMPRESSED:
            default:
                return new FileChannelByteReader(file);
        }
    }
}
