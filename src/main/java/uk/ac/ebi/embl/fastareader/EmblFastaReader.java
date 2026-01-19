package uk.ac.ebi.embl.fastareader;

import uk.ac.ebi.embl.fastareader.exception.FastaFileException;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

public class EmblFastaReader {

    public List<FastaEntry> fastaEntries = new ArrayList<>();

    /*
    public Optional<FastaEntry> setAccessionId(String submissionId, String accessionId) throws FastaFileException {
        Optional<FastaEntry> target = fastaEntries.stream()
                .filter(entry -> entry.getSubmissionId().equals(submissionId))
                .findFirst();
        target.ifPresent(entry -> entry.setAccessionId(accessionId));
        return target;
    }

    public Optional<FastaEntry> getFastaWithSubmissionId(String submissionId) throws FastaFileException {
        return fastaEntries.stream()
                .filter(entry -> entry.getSubmissionId().equals(submissionId))
                .findFirst();
    }
*/
}
