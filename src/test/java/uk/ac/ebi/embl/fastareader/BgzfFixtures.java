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
import java.nio.file.Files;
import java.nio.file.Path;
import uk.ac.ebi.embl.fastareader.io.BgzfTestWriter;

/** Test helper that produces BGZF copies of uncompressed test resources. */
final class BgzfFixtures {

    private BgzfFixtures() {}

    /** Creates a BGZF copy of the given test resource in {@code dir}, using {@code blockSize}. */
    static File bgzfCopyOf(String resourceDir, String fileName, Path dir, int blockSize) throws IOException {
        File source = TestResources.file(resourceDir, fileName);
        byte[] raw = Files.readAllBytes(source.toPath());
        Path target = dir.resolve(fileName + ".bgzf");
        return BgzfTestWriter.writeBgzfFile(target, raw, blockSize);
    }
}
