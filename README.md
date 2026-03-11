# kmtools

A toolkit for running, filtering, merging, and plotting [`km find_mutation`](https://github.com/iric-soft/km) results.

`kmtools` streamlines large-scale `km` workflows by allowing you to:

- Run `km find_mutation` across multiple threads via target sequence chunking
- Merge multiple output files into a single result
- Filter results against a reference variant list (with optional `--use-alt` mode)
- Generate summary plots for downstream analysis

---

## Installation

### From PyPI

```bash
pip install kmtools-bio
```

### Using pipx (recommended for CLI tools)

```bash
pipx install kmtools-bio
```

### From source

```bash
git clone https://github.com/trentzz/kmtools.git
cd kmtools
pip install .
```

### Prerequisites

- Python >= 3.10
- [`km`](https://github.com/iric-soft/km) must be installed and accessible in your `PATH` (required for `chunk` and `runall` commands).

To install `km` via `pipx`:

```bash
pipx install km-walk
```

---

## Quick Start

```bash
# Check version
kmtools --version

# Filter km results against a reference
kmtools filter \
    --reference variants.tsv \
    --km-output km_results.txt \
    --output filtered.tsv

# Merge multiple km output files
kmtools merge chunk_0.txt chunk_1.txt chunk_2.txt \
    --output merged.txt --keep

# Generate plots from filtered results
kmtools plot filtered.tsv --output-dir plots --charts all
```

---

## Commands

| Command          | Description                                                                         |
| ---------------- | ----------------------------------------------------------------------------------- |
| `kmtools chunk`  | Run `km find_mutation` in parallel across chunked target sequences                  |
| `kmtools merge`  | Merge multiple `km` output files into a single result                               |
| `kmtools filter` | Filter `km` results against a reference variant list                                |
| `kmtools plot`   | Generate visualisations from filtered results                                       |
| `kmtools runall` | Run the full pipeline (`chunk` -> `merge` -> `filter` -> `plot`) in a single command |

### Global Flags

| Flag        | Description            |
| ----------- | ---------------------- |
| `--version` | Show version and exit  |
| `--verbose` | Enable verbose logging |

---

### `kmtools chunk`

Run `km find_mutation` concurrently across multiple threads by splitting target sequences into pre-split subdirectories.

```bash
kmtools chunk \
    --threads 8 \
    --km-find-mutation-options "--ratio 0.0001" \
    --km-target-directory targets_split \
    --km-jellyfish-file database.jf \
    --output-dir output \
    --prefix km_output \
    --merge \
    --merge-output merged.txt \
    --merge-keep \
    --verbose
```

| Argument                     | Required | Description                                                             |
| ---------------------------- | -------- | ----------------------------------------------------------------------- |
| `--threads`                  | Yes      | Number of parallel `km` processes                                       |
| `--km-find-mutation-options` | Yes      | Options passed directly to `km find_mutation`                           |
| `--km-target-directory`      | Yes      | Directory containing pre-split target subdirectories (one per thread)   |
| `--km-jellyfish-file`        | Yes      | Path to the jellyfish k-mer counts file (`.jf`)                        |
| `--output-dir`               | No       | Directory to save output files (default: `.`)                           |
| `--prefix`                   | No       | Prefix for output files (default: `km_find_mutation_output`)            |
| `--merge`                    | No       | Automatically merge chunk outputs after processing                      |
| `--merge-output`             | No       | Filename for merged output (default: `km_find_mutation_merged_output.txt`) |
| `--merge-keep`               | No       | Keep intermediate chunk files after merging                             |

---

### `kmtools merge`

Combine multiple `km` output files into a single result file.

```bash
kmtools merge chunk_0.txt chunk_1.txt chunk_2.txt \
    --output merged.txt \
    --keep \
    --sort-by Min_coverage \
    --drop-duplicates \
    --verbose
```

| Argument            | Required | Description                                        |
| ------------------- | -------- | -------------------------------------------------- |
| `inputs`            | Yes      | Input files (positional, supports glob patterns)    |
| `--output`          | Yes      | Output file for merged results                     |
| `--keep`            | No       | Keep input files after merging (default: delete)    |
| `--sort-by`         | No       | Sort merged output by this column name             |
| `--drop-duplicates` | No       | Remove duplicate rows from merged output           |

---

### `kmtools filter`

Filter `km` output against a reference variant list to identify true positives.

```bash
kmtools filter \
    --reference variants.tsv \
    --km-output merged.txt \
    --output filtered.tsv \
    --output-type tsv \
    --count-threshold 2 \
    --verbose
```

With `--use-alt` mode (matches km sequences directly against provided ALT_SEQUENCE):

```bash
kmtools filter \
    --reference alt_variants.tsv \
    --km-output merged.txt \
    --output filtered.tsv \
    --use-alt
```

| Argument            | Required | Description                                                                 |
| ------------------- | -------- | --------------------------------------------------------------------------- |
| `--reference`       | Yes      | Reference file with known variants (`.csv` or `.tsv`)                       |
| `--km-output`       | Yes      | `km find_mutation` output file to filter                                    |
| `--output`          | Yes      | Output file for filtered results                                            |
| `--output-type`     | No       | Output format: `tsv`, `csv`, or `xlsx` (default: `tsv`)                    |
| `--count-threshold` | No       | Minimum k-mer count for a variant to pass (default: `2`)                    |
| `--use-alt`         | No       | Use alternate sequence mode (reference needs `CHROM`, `ALT_SEQUENCE`, `TYPE`) |

**Reference file format (standard mode):** TSV/CSV with columns `CHROM`, `POS`, `REF`, `ALT`, `TYPE`

**Reference file format (`--use-alt` mode):** TSV/CSV with columns `CHROM`, `ALT_SEQUENCE`, `TYPE`

See [docs/use-alt-mode.md](docs/use-alt-mode.md) for details on the `--use-alt` workflow.

---

### `kmtools plot`

Generate charts from filtered results.

```bash
kmtools plot filtered.tsv \
    --output-dir plots \
    --charts vaf,type,sample,overall \
    --verbose
```

| Argument       | Required | Description                                                                      |
| -------------- | -------- | -------------------------------------------------------------------------------- |
| `file`         | Yes      | Filtered results file (positional)                                               |
| `--output-dir` | No       | Directory to save plots (default: `.`)                                           |
| `--charts`     | No       | Comma-separated chart types: `vaf`, `type`, `sample`, `overall`, or `all` (default: `all`) |

**Available chart types:**

| Chart     | Output file             | Description                                       |
| --------- | ----------------------- | ------------------------------------------------- |
| `vaf`     | `vaf_distribution.png`  | Histogram of variant allele frequencies            |
| `type`    | `type_distribution.png` | Bar chart of variant type counts                   |
| `sample`  | `sample_summary.png`    | Stacked bar chart of found vs not-found per sample |
| `overall` | `overall_summary.png`   | Pie chart of overall detection summary             |

---

### `kmtools runall`

Run the full pipeline (`chunk` -> `merge` -> `filter` -> `plot`) in a single command.

```bash
kmtools runall \
    --threads 8 \
    --km-find-mutation-options "--ratio 0.0001" \
    --km-target-directory targets_split \
    --km-jellyfish-file database.jf \
    --merge-output merged.txt \
    --reference variants.tsv \
    --filtered-output filtered.tsv \
    --output-type tsv \
    --count-threshold 2 \
    --output-dir plots \
    --charts all \
    --verbose
```

| Argument                     | Required | Description                                       |
| ---------------------------- | -------- | ------------------------------------------------- |
| `--threads`                  | Yes      | Number of parallel threads for chunking            |
| `--km-find-mutation-options` | Yes      | Options passed to `km find_mutation`               |
| `--km-target-directory`      | Yes      | Directory containing pre-split target subdirectories |
| `--km-jellyfish-file`        | Yes      | Path to the jellyfish k-mer counts file (`.jf`)   |
| `--merge-output`             | Yes      | Output file for merged chunk results               |
| `--reference`                | Yes      | Reference file with known variants (`.csv`/`.tsv`) |
| `--filtered-output`          | Yes      | Output file for filtered results                   |
| `--output-type`              | No       | Output format for filtered results (default: `tsv`) |
| `--count-threshold`          | No       | Minimum k-mer count threshold (default: `2`)       |
| `--output-dir`               | No       | Directory to save plots (default: `.`)             |
| `--charts`                   | No       | Chart types to generate (default: `all`)           |

---

## Project Structure

```text
src/kmtools/
├── __init__.py
├── kmtools.py        # CLI entry point and argument parsing
├── chunk.py          # Parallel km find_mutation execution
├── merge.py          # Merge multiple km output files
├── filter.py         # Filter km results against reference variants
├── filter_types.py   # Dataclasses for filter types
├── plot.py           # Generate charts from filtered results
├── exceptions.py     # Custom exception classes
└── utils.py          # Shared utilities (logging, timing)
```

---

## Documentation

- [Getting Started](docs/getting-started.md) - Installation and quickstart
- [Command Reference](docs/commands.md) - Full argument tables for all commands
- [Use-Alt Mode](docs/use-alt-mode.md) - Guide for `--use-alt` filtering
- [File Formats](docs/file-formats.md) - All input/output file format specifications
- [Testing & Benchmarking](docs/testing-and-benchmarking.md) - Test suite and performance benchmarks

---

## Future Plans

- Automatic detection of chunk outputs in `runall`
- Support for alternate `km` subcommands beyond `find_mutation`
- Additional plot types and export formats
