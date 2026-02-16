package uk.ac.ebi.embl.fastareader.exception;

public class SequenceFileException extends Exception {
    public SequenceFileException() {}

    public SequenceFileException(String message) { super(message);}

    public SequenceFileException(Throwable cause) {super(cause);}

    public SequenceFileException(String message, Throwable cause) {super(message, cause);}
}
