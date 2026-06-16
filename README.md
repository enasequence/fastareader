# fastareader

A Java library for indexed read of FASTA files, providing efficient base-range sequence retrieval with explicit resource management.

## Features

- Random-access sequence slice retrieval (string or streaming `Reader`) with 1-based, inclusive base coordinates
- Multi-entry FASTA support via `FastaReader`; single-sequence plain-text files via `SequenceReader`
- In-memory `SequenceIndex` built once at open time for fast repeated range queries
- Gap region tracking (leading/trailing N-bases, internal gaps)
- Serialisation round-trip via `SequenceInfoDTO` / `SequenceInfoStore`
- **Transparent `bgzip`-compressed input** — BGZF files are decompressed on demand; no code changes required

## Supported input formats

| Format | Extension (typical) | Supported |
|---|---|---|
| Plain FASTA (UTF-8) | `.fa`, `.fasta`, `.txt` | ✅ |
| Plain sequence (UTF-8) | `.txt`, `.seq` | ✅ |
| `bgzip`-compressed FASTA | `.fa.gz`, `.fasta.gz` | ✅ |
| `bgzip`-compressed sequence | `.txt.gz` | ✅ |
| Plain gzip | `.fa.gz`, `.gz` | ❌ (see below) |

## BGZF / bgzip support

`FastaReader` and `SequenceReader` automatically detect and transparently decompress
[BGZF](https://samtools.github.io/hts-specs/SAMv1.pdf)-compressed files (the block-gzip
format produced by `bgzip` from samtools/htslib). No changes to calling code are needed:

```java
// Uncompressed — works as before
try (FastaReader reader = new FastaReader(new File("genome.fa"))) {
    String slice = reader.getSequenceSlice(0, 1, 100);
}

// bgzip-compressed — identical API, identical results
try (FastaReader reader = new FastaReader(new File("genome.fa.gz"))) {
    String slice = reader.getSequenceSlice(0, 1, 100);
}
```

Detection is automatic: the library reads the first few bytes and checks for the BGZF
`BC` extra subfield in the gzip header. No file-extension convention is assumed.

### How it works

A BGZF file is a sequence of independent, fixed-size gzip blocks (≤ 64 KiB of
decompressed payload each). The library:

1. Walks the block headers at open time to build an in-memory **block directory**
   mapping logical (uncompressed) byte offsets to compressed block locations.
2. Satisfies all reads against the **logical uncompressed stream** — the same
   flat, additive address space used for plain files — so every internal
   offset calculation is unchanged.
3. Decompresses individual blocks on demand using a small **LRU block cache**
   (16 blocks, ≤ ~1 MiB), reusing cached payloads for adjacent or repeated slice queries.

As a result the `SequenceIndex` produced from a BGZF file is **identical** to the one
produced from its decompressed equivalent, and serialised indexes (`SequenceInfoDTO`)
round-trip correctly.

### Plain gzip is not supported

Plain gzip files (`.gz` files *not* produced by `bgzip`) are **not seekable** and are
therefore rejected at open time with a clear error:

```
Plain gzip is not supported; recompress with bgzip: /path/to/file.fa.gz
```

To convert a plain gzip file:

```bash
gunzip file.fa.gz
bgzip  file.fa          # produces file.fa.gz in BGZF format
```

### No external dependencies

BGZF decompression is implemented entirely using the JDK standard library
(`java.util.zip.Inflater` in nowrap mode). No HTSJDK, no Hadoop, no additional
runtime dependencies are introduced.

## Testing

```bash
./gradlew spotlessApply test
```
