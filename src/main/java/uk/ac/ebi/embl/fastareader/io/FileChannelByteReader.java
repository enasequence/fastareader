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
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.StandardOpenOption;

/**
 * {@link SeekableByteReader} adapter over a plain (uncompressed) {@link FileChannel}.
 *
 * <p>All operations delegate 1-to-1 to the underlying channel. The logical stream is
 * identical to the raw file bytes, so logical offsets equal physical file offsets.</p>
 */
public final class FileChannelByteReader implements SeekableByteReader {

    private final FileChannel channel;

    /**
     * Opens the given file for reading and wraps it in a {@code FileChannelByteReader}.
     * The caller is responsible for closing this reader when done.
     */
    public FileChannelByteReader(File file) throws IOException {
        this.channel = FileChannel.open(file.toPath(), StandardOpenOption.READ);
    }

    /**
     * Wraps an already-open {@link FileChannel}. Ownership of the channel is transferred
     * to this reader; closing this reader closes the channel.
     */
    public FileChannelByteReader(FileChannel channel) {
        this.channel = channel;
    }

    @Override
    public int read(ByteBuffer dst, long position) throws IOException {
        return channel.read(dst, position);
    }

    @Override
    public long size() throws IOException {
        return channel.size();
    }

    @Override
    public long position() throws IOException {
        return channel.position();
    }

    @Override
    public SeekableByteReader position(long newPosition) throws IOException {
        channel.position(newPosition);
        return this;
    }

    @Override
    public boolean isOpen() {
        return channel.isOpen();
    }

    @Override
    public void close() throws IOException {
        channel.close();
    }
}
