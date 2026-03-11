# Testing & Benchmarking

## Test suite overview

kmtools uses [pytest](https://pytest.org/) for testing. The test suite covers unit tests, comprehensive component tests, and integration tests.

### Running tests

```bash
# Run all tests
poetry run pytest -v

# Run a specific test file
poetry run pytest tests/test_filter.py -v

# Run tests matching a pattern
poetry run pytest -k "test_merge" -v
```

### Test structure

```
tests/
├── conftest.py                  # Shared fixtures
├── data/                        # Realistic test data files
│   ├── reference_substitution.tsv
│   ├── reference_deletion.tsv
│   ├── reference_mixed.tsv
│   ├── alt_reference.tsv
│   ├── km_output_substitution.txt
│   ├── km_output_deletion.txt
│   └── km_output_mixed.txt
├── test_chunk.py                # Chunk class unit tests
├── test_cli.py                  # CLI argument parsing tests
├── test_data_integration.py     # Tests using data/ files
├── test_filter.py               # Filter calculation tests
├── test_filter_comprehensive.py # Comprehensive filter component tests
├── test_filter_types.py         # Dataclass tests
├── test_integration.py          # End-to-end integration tests
├── test_merge.py                # Merge class unit tests
├── test_use_alt.py              # Use-alt mode tests
└── test_utils.py                # Utility function tests
```

### Test categories

| Category | File(s) | Tests | What it covers |
|----------|---------|------:|----------------|
| Unit - Chunk | `test_chunk.py` | 11 | km validation, target splitting, subprocess calls |
| Unit - CLI | `test_cli.py` | 7 | Argument parsing for all subcommands |
| Unit - Filter | `test_filter.py`, `test_filter_comprehensive.py` | 25 | Filter logic, reference validation, output writing |
| Unit - Filter types | `test_filter_types.py` | 6 | Dataclass construction |
| Unit - Merge | `test_merge.py` | 6 | File merging, glob patterns, keep/delete |
| Unit - Use-alt | `test_use_alt.py` | 16 | Alt reference validation, alt filtering, CLI flag |
| Unit - Utils | `test_utils.py` | 6 | Logging, timing |
| Integration | `test_integration.py` | 34 | End-to-end pipelines, CLI invocation, error handling |
| Data integration | `test_data_integration.py` | 10 | Tests using realistic data files |
| **Total** | | **~136** | |

### Test data files

The `tests/data/` directory contains biologically plausible test data with internally consistent sequences:

- **reference_substitution.tsv** - 3 substitution variants on chr1 and chr2
- **reference_deletion.tsv** - 2 deletion variants on chr3 and chr4
- **reference_mixed.tsv** - 4 variants covering substitution, deletion, and insertion types
- **alt_reference.tsv** - 3 alt sequences for use-alt mode testing
- **km_output_substitution.txt** - 6 km rows including passing, low-coverage, and non-vs_ref rows
- **km_output_deletion.txt** - 4 km rows for deletion variants
- **km_output_mixed.txt** - 10 km rows covering all variant types

All sequences are generated with guaranteed internal consistency: the variant offset, reference sequence, and variant sequence are mathematically verified to be correct.

---

## Benchmark results

Benchmarks run on the development machine. Results will vary by hardware.

To run benchmarks yourself:

```bash
poetry run python benchmarks/run_benchmarks.py
```

### Merge performance

| Files | Rows/File | Total Rows | Time (s) | Throughput (rows/s) |
|------:|----------:|-----------:|---------:|--------------------:|
| 10 | 10 | 100 | 0.1305 | 766 |
| 50 | 10 | 500 | 0.0653 | 7,661 |
| 100 | 10 | 1,000 | 0.1632 | 6,129 |
| 500 | 10 | 5,000 | 0.6281 | 7,961 |

### Merge with sort-by + drop-duplicates

| Files | Rows/File | Total Rows | Time (s) | Throughput (rows/s) |
|------:|----------:|-----------:|---------:|--------------------:|
| 10 | 10 | 100 | 0.1054 | 948 |
| 50 | 10 | 500 | 0.0876 | 5,705 |
| 100 | 10 | 1,000 | 0.1172 | 8,529 |
| 500 | 10 | 5,000 | 0.9081 | 5,506 |

### Filter performance

| KM Rows | Ref Variants | Time (s) | Throughput (rows/s) |
|--------:|-------------:|---------:|--------------------:|
| 10 | 5 | 0.0126 | 794 |
| 50 | 5 | 0.0090 | 5,570 |
| 100 | 5 | 0.0076 | 13,210 |
| 500 | 5 | 0.0070 | 71,572 |
| 1,000 | 5 | 0.0136 | 73,590 |

### Notes

- Merge throughput is dominated by file I/O overhead (reading many small files). Larger files yield better throughput.
- Filter throughput scales well with input size since the inner loop is pandas-based.
- The `--sort-by` and `--drop-duplicates` options add modest overhead to merge operations.
