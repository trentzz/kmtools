# Getting Started

## Prerequisites

- Python >= 3.10
- [`km`](https://github.com/iric-soft/km) must be installed and accessible in your `PATH` (required only for `chunk` and `runall` commands)

## Installation

### From PyPI (recommended)

```bash
pip install kmtools
```

### Using pipx (isolated environment)

```bash
pipx install kmtools
```

### From source

```bash
git clone https://github.com/trentzz/kmtools.git
cd kmtools
pip install .
```

### Development install

```bash
git clone https://github.com/trentzz/kmtools.git
cd kmtools
pip install poetry
poetry install
```

## Verify installation

```bash
kmtools --version
kmtools --help
```

## Quickstart

### 1. Filter km results against a reference

Given a reference TSV (`variants.tsv`) with known variants and km output (`km_results.txt`):

```bash
kmtools filter \
    --reference variants.tsv \
    --km-output km_results.txt \
    --output filtered.tsv
```

### 2. Merge multiple km output files

```bash
kmtools merge chunk_0.txt chunk_1.txt chunk_2.txt \
    --output merged.txt \
    --keep
```

### 3. Generate plots

```bash
kmtools plot filtered.tsv --output-dir plots --charts all
```

### 4. Full pipeline (requires km installed)

```bash
kmtools runall \
    --threads 8 \
    --km-find-mutation-options "--ratio 0.0001" \
    --km-target-directory targets_split \
    --km-jellyfish-file database.jf \
    --merge-output merged.txt \
    --reference variants.tsv \
    --filtered-output filtered.tsv \
    --output-dir plots \
    --charts all \
    --verbose
```

## Next steps

- [Command Reference](commands.md) - Full argument documentation
- [File Formats](file-formats.md) - Input/output file specifications
- [Use-Alt Mode](use-alt-mode.md) - Alternate sequence matching workflow
