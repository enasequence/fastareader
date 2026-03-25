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
        return MAPPER.readValue(inputPath.toFile(),
                MAPPER.getTypeFactory().constructCollectionType(List.class, SequenceInfoDTO.class));
    }
}
