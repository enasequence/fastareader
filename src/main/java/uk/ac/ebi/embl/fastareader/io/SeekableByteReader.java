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

import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * Minimal abstraction over the subset of {@code FileChannel} used by this library, expressed
 * against the <em>logical (uncompressed)</em> byte stream.
 *
 * <p>For plain files the logical stream is the file itself; for BGZF-compressed files the
 * logical stream is the concatenation of all decompressed block payloads. All positional
 * arithmetic in the library operates on logical offsets, which are flat and additive.</p>
 *
 * <p>Mirrors the positional-read contract of {@link java.nio.channels.FileChannel}:
 * {@link #read(ByteBuffer, long)} does <em>not</em> affect the stateful cursor.</p>
 */
public interface SeekableByteReader extends AutoCloseable {

    /**
     * Absolute positional read against the logical (uncompressed) stream.
     * Mirrors {@code FileChannel.read(ByteBuffer, long)}: returns the number of bytes read,
     * or {@code -1} if {@code position} is at or past EOF. Does <em>not</em> affect
     * {@link #position()}.
     */
    int read(ByteBuffer dst, long position) throws IOException;

    /** Total length of the logical (uncompressed) stream in bytes. */
    long size() throws IOException;

    /** Returns the current stateful cursor position. */
    long position() throws IOException;

    /**
     * Sets the stateful cursor to {@code newPosition} and returns {@code this} for chaining.
     * Does not perform any I/O.
     */
    SeekableByteReader position(long newPosition) throws IOException;

    /** Returns {@code true} if this reader is open and ready to serve reads. */
    boolean isOpen();

    @Override
    void close() throws IOException;
}
