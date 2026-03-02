# Command Reference

All commands support the `--verbose` global flag for detailed logging.

```bash
kmtools [--version] [--verbose] <command> [args]
```

---

## `kmtools chunk`

Run `km find_mutation` in parallel across chunked target sequences.

```bash
kmtools chunk \
    --threads <N> \
    --km-find-mutation-options "<options>" \
    --km-target-directory <dir> \
    --km-jellyfish-file <file.jf> \
    [--output-dir <dir>] \
    [--prefix <prefix>] \
    [--merge] \
    [--merge-output <file>] \
    [--merge-keep]
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--threads` | Yes | - | Number of parallel `km` processes |
| `--km-find-mutation-options` | Yes | - | Options passed directly to `km find_mutation` (quote the string) |
| `--km-target-directory` | Yes | - | Directory containing pre-split target subdirectories. Must have exactly `--threads` subdirectories with roughly equal file counts. |
| `--km-jellyfish-file` | Yes | - | Path to the jellyfish k-mer counts file (`.jf`) |
| `--output-dir` | No | `.` | Directory to save chunk output files |
| `--prefix` | No | `km_find_mutation_output` | Prefix for output filenames (produces `<prefix>_0.txt`, `<prefix>_1.txt`, etc.) |
| `--merge` | No | `false` | Merge chunk outputs into a single file after processing |
| `--merge-output` | No | `km_find_mutation_merged_output.txt` | Filename for merged output when `--merge` is used |
| `--merge-keep` | No | `false` | Keep intermediate chunk files after merging |

### How it works

1. Validates that `km` is installed and the target directory has the correct structure
2. Splits work across threads (one subdirectory per thread)
3. Runs `km find_mutation` on each chunk in parallel via `ThreadPoolExecutor`
4. Optionally merges all results when complete

---

## `kmtools merge`

Combine multiple `km` output files into a single result.

```bash
kmtools merge <input_files...> \
    --output <file> \
    [--keep] \
    [--sort-by <column>] \
    [--drop-duplicates]
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `inputs` | Yes | - | Input files (positional arguments, supports glob patterns) |
| `--output` | Yes | - | Output file for merged results |
| `--keep` | No | `false` | Keep input files after merging (default: delete them) |
| `--sort-by` | No | `None` | Sort merged output by this column name |
| `--drop-duplicates` | No | `false` | Remove exact duplicate rows from merged output |

### Notes

- Input files are read as tab-separated values
- Column consistency is validated across files (mismatches produce warnings)
- Without `--keep`, input files are deleted after successful merge

---

## `kmtools filter`

Filter `km find_mutation` output against a reference variant list.

```bash
kmtools filter \
    --reference <file> \
    --km-output <file> \
    --output <file> \
    [--output-type tsv|csv|xlsx] \
    [--count-threshold <N>] \
    [--use-alt]
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--reference` | Yes | - | Reference file with known variants (`.csv` or `.tsv`) |
| `--km-output` | Yes | - | `km find_mutation` output file to filter |
| `--output` | Yes | - | Output file for filtered results |
| `--output-type` | No | `tsv` | Output format: `tsv`, `csv`, or `xlsx` |
| `--count-threshold` | No | `2` | Minimum k-mer count (`Min_coverage`) for a variant to pass |
| `--use-alt` | No | `false` | Enable alternate sequence matching mode |

### Filter conditions (standard mode)

A km row passes the filter when ALL of these are true:

1. **TYPE match**: km `Type` matches reference `TYPE`
2. **COUNT threshold**: km `Min_coverage` >= `--count-threshold`
3. **ALT match**: km `Sequence` matches the calculated expected variant sequence
4. **Not reference type**: reference `TYPE` is not `"Reference"`
5. **vs_ref info**: km `Info` field is `"vs_ref"`

### Filter conditions (use-alt mode)

1. **TYPE match**: km `Type` matches alt reference `TYPE`
2. **COUNT threshold**: km `Min_coverage` >= `--count-threshold`
3. **Sequence match**: km `Sequence` matches `ALT_SEQUENCE` from the reference
4. **vs_ref info**: km `Info` field is `"vs_ref"`

### Output columns

**Standard mode:** `SAMPLE`, `CHROM`, `POS`, `REF`, `ALT`, `TYPE`, `FOUND`, `FILTER_NOTES`, `KMER_VAF`, `KMER_MIN_COVERAGE`, `KMER_EXPRESSION`, `REF_SEQUENCE`, `VARIANT_SEQUENCE`

**Use-alt mode:** `SAMPLE`, `CHROM`, `ALT_SEQUENCE`, `TYPE`, `FOUND`, `FILTER_NOTES`, `KMER_VAF`, `KMER_MIN_COVERAGE`, `KMER_EXPRESSION`, `REF_SEQUENCE`, `VARIANT_SEQUENCE`

---

## `kmtools plot`

Generate charts from filtered results.

```bash
kmtools plot <file> \
    [--output-dir <dir>] \
    [--charts <types>]
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `file` | Yes | - | Filtered results file (positional) |
| `--output-dir` | No | `.` | Directory to save plot images |
| `--charts` | No | `all` | Comma-separated chart types |

### Chart types

| Type | Output file | Description |
|------|-------------|-------------|
| `vaf` | `vaf_distribution.png` | Histogram of variant allele frequencies (rVAF) for found variants |
| `type` | `type_distribution.png` | Bar chart of variant type counts |
| `sample` | `sample_summary.png` | Stacked bar: found vs not-found per sample |
| `overall` | `overall_summary.png` | Pie chart of overall found vs not-found ratio |
| `all` | All of the above | Generate all chart types |

---

## `kmtools runall`

Run the complete pipeline in one command: `chunk` -> `merge` -> `filter` -> `plot`.

```bash
kmtools runall \
    --threads <N> \
    --km-find-mutation-options "<options>" \
    --km-target-directory <dir> \
    --km-jellyfish-file <file.jf> \
    --merge-output <file> \
    --reference <file> \
    --filtered-output <file> \
    [--output-type tsv|csv|xlsx] \
    [--count-threshold <N>] \
    [--output-dir <dir>] \
    [--charts <types>]
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--threads` | Yes | - | Number of parallel threads for chunking |
| `--km-find-mutation-options` | Yes | - | Options passed to `km find_mutation` |
| `--km-target-directory` | Yes | - | Directory with pre-split target subdirectories |
| `--km-jellyfish-file` | Yes | - | Path to the jellyfish k-mer counts file |
| `--merge-output` | Yes | - | Output file for merged chunk results |
| `--reference` | Yes | - | Reference file with known variants |
| `--filtered-output` | Yes | - | Output file for filtered results |
| `--output-type` | No | `tsv` | Output format for filtered results |
| `--count-threshold` | No | `2` | Minimum k-mer count threshold |
| `--output-dir` | No | `.` | Directory to save plots |
| `--charts` | No | `all` | Chart types to generate |

### Pipeline steps

1. **Chunk**: Split targets and run `km find_mutation` in parallel
2. **Merge**: Combine chunk outputs into `--merge-output`
3. **Filter**: Filter merged results against `--reference` into `--filtered-output`
4. **Plot**: Generate charts from filtered results into `--output-dir`
