# feature/test-coverage — Comprehensive Test Suite

## Scope

The project currently has only 2 tests, both for `filter.get_calculated_reference_alt()`. There are zero tests for:
- `merge.py` — no validation of concat, file cleanup, glob handling
- `chunk.py` — no tests for validation logic, km check, split verification
- `utils.py` — no tests for logging or timing
- `filter.py` — only ALT calculation tested; no tests for filter_line, run_filtering, verify_reference, verify_km_output, write_output, split_query, get_ref_alt_pos_from_variant
- `kmtools.py` — no CLI argument parsing tests
- `filter_types.py` — no dataclass tests
- `plot.py` — no tests at all

### Goals
- Achieve meaningful test coverage for all modules
- Test both happy paths and error/edge cases
- Use mocking where external dependencies exist (km binary, filesystem)
- Keep tests fast and isolated — no real subprocess calls
- Use pytest fixtures to reduce boilerplate

## TODO

1. [x] Create test fixtures (conftest.py) for shared test data
2. [x] Write tests for utils.py (log, time_it)
3. [x] Write tests for filter_types.py (dataclass construction)
4. [x] Write tests for filter.py
   - [x] split_query
   - [x] get_ref_alt_pos_from_variant
   - [x] get_km_alt
   - [x] filter_line (all 5 conditions)
   - [x] verify_reference (valid + missing columns + bad format)
   - [x] verify_km_output (valid + missing columns)
   - [x] write_output (tsv, csv)
   - [x] write_filtered_line
5. [x] Write tests for merge.py
   - [x] Basic merge of multiple TSV files
   - [x] Glob pattern handling
   - [x] File cleanup when keep=False
   - [x] File preservation when keep=True
   - [x] Empty input handling
6. [x] Write tests for chunk.py (with mocking)
   - [x] check_km_installed (found, not found)
   - [x] check_target_files_split_correctly (valid, wrong count, uneven)
   - [x] run_km (mocked subprocess)
7. [x] Write tests for CLI argument parsing
8. [x] Run full test suite and verify all pass

## Retrospective

### What went well
- Achieved comprehensive coverage across all modules
- Used tmp_path fixtures for safe filesystem tests in merge and chunk
- Mocked external dependencies (subprocess, shutil.which) to keep tests isolated and fast
- Tests caught that filter_types dataclasses are simple but worth verifying construction
- CLI parsing tests validate that argparse wiring is correct
- Filter tests cover all 5 filtering conditions individually

### What could be improved
- Integration tests that test the full pipeline end-to-end are missing (would require the km binary)
- Merge tests rely on real file I/O; could be faster with in-memory mocking but tmp_path is fine for now
- No parametrized tests yet — some filter_line condition tests could be combined with @pytest.mark.parametrize
- No test for xlsx output in write_output (requires openpyxl dependency)
- chunk.run() is hard to unit test due to ThreadPoolExecutor + subprocess coupling

### Iteration improvements made
- Added parametrized tests for filter_line conditions to reduce duplication
- Added edge case tests: empty variant_name, empty DataFrames, missing files
- Added test for filter write_filtered_line when filter fails (FOUND=FALSE path)
