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
package uk.ac.ebi.embl.fastareader.api.rereading;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class SequenceInfoStore {

    private static final ObjectMapper MAPPER = new ObjectMapper();

    public void write(List<SequenceInfoDTO> infos, Path outputPath) throws IOException {
        Files.createDirectories(outputPath.getParent());
        MAPPER.writerWithDefaultPrettyPrinter().writeValue(outputPath.toFile(), infos);
    }

    public List<SequenceInfoDTO> read(Path inputPath) throws IOException {
        if (!Files.exists(inputPath) || !Files.isRegularFile(inputPath)) {
            throw new IllegalArgumentException("Not a valid file: " + inputPath);
        }
        return MAPPER.readValue(
                inputPath.toFile(), MAPPER.getTypeFactory().constructCollectionType(List.class, SequenceInfoDTO.class));
    }
}
