# Developer Guide for kmtools

This document is intended for contributors and developers working on `kmtools`.

---

## Development Setup

```bash
# Clone and enter the repository
git clone https://github.com/trentzz/kmtools.git
cd kmtools

# Install Poetry if not already installed
pip install poetry

# Install dependencies (including dev dependencies)
poetry install

# Run the CLI directly
poetry run kmtools --help
```

To install your development version globally via pipx for testing:

```bash
pipx install --force .
```

---

## Project Layout

```text
src/kmtools/
├── __init__.py          # Package init
├── kmtools.py           # CLI entry point and argument parsing
├── chunk.py             # Chunk class: parallel km find_mutation execution
├── merge.py             # Merge class: combine multiple km output files
├── filter.py            # Filter class: filter km results against reference
├── filter_types.py      # Dataclasses (FilterCondition, FilterResult, etc.)
├── plot.py              # Plot class: generate charts from filtered results
├── exceptions.py        # Custom exception hierarchy (KmtoolsError, etc.)
└── utils.py             # Utility functions (logging, timing)

tests/
├── conftest.py                  # Shared test fixtures
├── data/                        # Realistic test data files (7 files)
│   ├── reference_substitution.tsv
│   ├── reference_deletion.tsv
│   ├── reference_mixed.tsv
│   ├── alt_reference.tsv
│   ├── km_output_substitution.txt
│   ├── km_output_deletion.txt
│   └── km_output_mixed.txt
├── test_chunk.py                # Chunk unit tests
├── test_cli.py                  # CLI argument parsing tests
├── test_data_integration.py     # Tests using data/ files
├── test_filter.py               # Filter calculation tests
├── test_filter_comprehensive.py # Comprehensive filter component tests
├── test_filter_types.py         # Dataclass tests
├── test_integration.py          # End-to-end integration tests (34 tests)
├── test_merge.py                # Merge unit tests
├── test_use_alt.py              # Use-alt mode tests
└── test_utils.py                # Utility function tests

docs/
├── getting-started.md           # Installation and quickstart
├── commands.md                  # Full command reference with argument tables
├── use-alt-mode.md              # --use-alt guide
├── file-formats.md              # All file format specifications
└── testing-and-benchmarking.md  # Test structure and benchmark results

benchmarks/
└── run_benchmarks.py            # Merge and filter performance benchmarks
```

---

## Code Style and Conventions

* Code must be **PEP8 compliant**.
* Use **type hints** for all public functions and class methods.
* Avoid inline lambdas for subcommand functions; use named `run_*` functions for readability.
* All logging and timing should be handled through `Utils.log()` to maintain consistent formatting.
* Custom exceptions inherit from `KmtoolsError` in `exceptions.py`.

---

## Testing

The test suite has **136+ tests** across 11 test files.

```bash
# Run all tests
poetry run pytest -v

# Run with short traceback
poetry run pytest --tb=short

# Run specific test file
poetry run pytest tests/test_integration.py -v

# Run tests matching a keyword
poetry run pytest -k "test_merge" -v
```

### Test categories

| Category | Tests | Description |
|----------|------:|-------------|
| Unit tests | ~60 | Individual component tests (chunk, merge, filter, plot, utils) |
| Use-alt tests | 16 | Alt reference validation, alt filtering, CLI integration |
| Integration tests | 34 | End-to-end pipelines, CLI invocation, error handling |
| Data integration | 10 | Tests using realistic data files from `tests/data/` |

### Test data files

The `tests/data/` directory contains 7 biologically plausible test data files with internally consistent sequences. See [docs/testing-and-benchmarking.md](docs/testing-and-benchmarking.md) for details.

### Benchmarks

```bash
poetry run python benchmarks/run_benchmarks.py
```

---

## Release Process

1. Ensure all tests pass: `poetry run pytest -v`
2. Build the package: `poetry build`
3. Verify the wheel: `pip install dist/kmtools-*.whl && kmtools --version`
4. Tag and push a release:

```bash
git tag -a vX.Y.Z -m "Release vX.Y.Z"
git push origin vX.Y.Z
```

5. Publish to PyPI:

```bash
poetry publish --build
```

---

## Contributing

* Use feature branches: `feature/<short-description>`
* Open a pull request with a clear description of the change.
* Keep commits focused and atomic.
* Run `poetry run pytest -v` before opening a PR.
