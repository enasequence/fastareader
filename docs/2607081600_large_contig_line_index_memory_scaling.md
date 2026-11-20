# Issue: per-line/per-gap in-memory index does not scale to gigabase-scale single contigs

- Document Date: 2026-07-08
- Status: reported, not yet fixed or scoped

## Summary

`SequenceIndexBuilder` builds, and `SequenceIndex` retains for the lifetime of the reader, one
`LineEntry` object per physical FASTA sequence line and one `GapRegion` object per contiguous
gap-base (`N`/`n`) run. Both are plain `ArrayList`s held fully in memory with no eviction or
paging. For a single FASTA record whose length approaches or exceeds roughly 1–2 Gbp, the
`LineEntry` list alone can consume multiple gigabytes of JVM heap, causing consumers to hit
`OutOfMemoryError` **before any actual sequence data is read or written**, regardless of how
memory-efficient the consumer's own read/write logic is.

This is orthogonal to (and, in practice, ended up masking) unrelated memory-efficiency work in a
downstream consumer: see "How this was found" below.

## Root cause

`SequenceIndexBuilder.buildFrom` (`src/main/java/uk/ac/ebi/embl/fastareader/sequenceutils/SequenceIndexBuilder.java:43`)
scans the whole FASTA record once and, for every physical line, calls `commitOpenLineIfAny`
(`SequenceIndexBuilder.java:193`) which appends a `LineEntry` to `ScanState.lines`
(`SequenceIndexBuilder.java:90`, an `ArrayList<LineEntry>`). Likewise, every contiguous run of gap
bases (`N`/`n`, or whatever the configured alphabet treats as a gap) is committed as a `GapRegion`
via `commitOpenGapIfAny` (`SequenceIndexBuilder.java:186`) into `ScanState.gapRegions`. Both lists
are handed unmodified into the resulting `SequenceIndex` (`SequenceIndexBuilder.java:63-72`) and
stored as `SequenceIndex.lines` / `SequenceIndex.gapRegions`
(`src/main/java/uk/ac/ebi/embl/fastareader/sequenceutils/SequenceIndex.java:20,25`), where they
live for the lifetime of the reader and are used for O(log n) byte-offset lookups
(`SequenceIndex.findLineByBase`, `SequenceIndex.java:213`, used by
`byteSpanForBaseRangeIncludingEdgeNBases`, `SequenceIndex.java:183`).

Per-instance cost (JVM object header + field alignment, no compressed-oops assumptions beyond
typical 64-bit defaults):

- `LineEntry` (`LineEntry.java`) — 4 `long` fields → ≈ 48 bytes/instance.
- `GapRegion` (`GapRegion.java`) — 2 `long` fields → ≈ 32 bytes/instance.
- Plus the backing `ArrayList<E>` reference array: ~8 bytes/entry (amortized over doubling growth).

For a FASTA file wrapped at the common 60-characters-per-line convention, the number of
`LineEntry` instances is `totalBases / 60`. The line index alone therefore scales as
approximately:

```
(48 + 8) bytes/line × (totalBases / 60 lines) ≈ 0.93 bytes of JVM heap per base of sequence
```

For a single 2.5 Gbp contig that is **≈ 2.3 GB for the line index alone**, before accounting for
anything else the consuming process needs (other JVM state, GC headroom, `gapRegions`, etc.).

## How this was found

While validating a new sequence-streaming write path in `gff3tools` (PR
[enasequence/gff3tools#148](https://github.com/enasequence/gff3tools/pull/148), design doc
`gff3tools/docs/2607081330_stream_entry_sequence_write.md`, currently consuming
`uk.ac.ebi.ena:fastareader:1.2.0`), a synthetic 2.5 Gbp single-contig FASTA + matching GFF3 was
generated to smoke-test the CLI conversion tool end-to-end at "big file" scale.

1. `gff3tools conversion --sequence big.fasta big.gff3 out.embl` was run at `-Xmx512m`, `768m`,
   `1200m`, and `2G`. **All four runs hit `OutOfMemoryError`** (`CLIExitCode.OUT_OF_MEMORY`), with
   peak RSS scaling roughly linearly with `-Xmx` (i.e. the process consumes whatever heap it is
   given and then dies, rather than plateauing at some fixed working set).
2. Live process monitoring showed the crash happens **before any output is written** — the
   destination `.embl` file remained at 0 bytes for the entire run in every case.
3. A `jcmd <pid> GC.class_histogram` snapshot taken mid-run (before the JVM subsequently OOM'd)
   showed the two dominant live object types:

   ```
   num     #instances         #bytes  class name
   1:      11,249,137     539,958,576  uk.ac.ebi.embl.fastareader.sequenceutils.LineEntry
   2:      15,451,294     494,441,408  uk.ac.ebi.embl.fastareader.sequenceutils.GapRegion
   ```

   (measured mid-load — the full record has ~41.6M lines, so both counts were still climbing when
   the process ran out of heap.)

   Note on the `GapRegion` count: the first version of the synthetic-data generator scattered
   single-base `N`s uniformly at random (~2.3% of bases) rather than in realistic contiguous
   assembly-gap blocks, which pathologically inflated the number of (very short) gap regions
   relative to real genomic data. After fixing the generator to produce ~1 gap block per 20 Mbp
   (125 gap blocks total for 2.5 Gbp — a much more realistic assembly-gap density), the OOM still
   reproduced: `LineEntry` alone (≈ 48 bytes × 41.6M ≈ 2.3 GB) is sufficient to explain it at up to
   2 GB heap, independent of gap-region count.
4. As a control, the same tool converting a more realistic 300 Mbp single contig (comparable to
   the largest real human chromosome, chr1 ≈ 248 Mbp) succeeded consistently from ~384 MB heap
   upward. The per-line index cost there (~4.1M lines × 56 bytes ≈ 230 MB) is proportionally
   significant but tolerable — this scale is where the downstream gff3tools streaming-write
   optimisation this investigation was validating actually showed its intended ~4x reduction in
   minimum required heap versus the prior non-streaming code path.

## Impact

- Any consumer reading a single FASTA record whose length approaches or exceeds roughly 1–2 Gbp
  needs proportionally large heap (~1 byte of JVM heap per base, just for the line index) no
  matter how memory-conscious the consumer's own logic is. This is uncommon for vertebrate
  chromosomes but real for some plant/amphibian genomes (e.g. individual wheat, axolotl, and
  lungfish chromosomes are known to exceed 1 Gbp).
- It masks and complicates testing of memory-efficiency work in downstream consumers: a
  read-side memory ceiling that scales with input size defeats the purpose of a bounded-memory
  write path, and made it impossible to fully validate the gff3tools streaming write path at true
  gigabase scale within this investigation.
- The current design (one `LineEntry` per physical line, enabling O(log n) byte-offset lookup via
  `SequenceIndex.findLineByBase`) is reasonable for typical chromosome-scale inputs but does not
  scale to arbitrarily large single contigs.

## Possible directions (not evaluated in depth — for discussion, not a committed plan)

1. **Regular-line-width fast path.** The overwhelming majority of well-formed FASTA files use a
   single fixed line width for every line in a record except (sometimes) the last line. If
   `SequenceIndexBuilder` detected that a record's lines all share one constant width, it could
   store `(lineWidth, totalBases, byteStart)` instead of one `LineEntry` per line, computing any
   line's byte offset arithmetically instead of storing it, and falling back to explicit per-line
   entries only for the irregular tail line (or degrading gracefully to a full per-line index if
   irregular widths are detected mid-record). This would take the dominant cost from O(lines) to
   O(1) for the common case.
2. **Chunked/sparse indexing.** Store an explicit `LineEntry` only every N lines (e.g. every
   1024th), plus the regular line width between checkpoints, trading a small amount of extra
   linear scanning per lookup for a large (~1000x) reduction in index memory.
3. **Gap-region memory** is a secondary concern for realistic biological data (assembly gaps are
   naturally few and long), but pathologically fragmented inputs (e.g. densely scattered
   ambiguity bases) could still produce a large `gapRegions` list. The same
   chunking/summarization idea could apply if this proves to be a real-world problem rather than
   a synthetic-data artifact.

## References

- `src/main/java/uk/ac/ebi/embl/fastareader/sequenceutils/SequenceIndexBuilder.java` —
  `buildFrom:43`, `commitOpenLineIfAny:193`, `commitOpenGapIfAny:186`
- `src/main/java/uk/ac/ebi/embl/fastareader/sequenceutils/SequenceIndex.java` — `lines:20`,
  `gapRegions:25`, `findLineByBase:213`, `byteSpanForBaseRangeIncludingEdgeNBases:183`
- `src/main/java/uk/ac/ebi/embl/fastareader/sequenceutils/LineEntry.java`,
  `src/main/java/uk/ac/ebi/embl/fastareader/sequenceutils/GapRegion.java`
- Discovered while validating
  [enasequence/gff3tools#148](https://github.com/enasequence/gff3tools/pull/148)
  (`gff3tools/docs/2607081330_stream_entry_sequence_write.md`)
