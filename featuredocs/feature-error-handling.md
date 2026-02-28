# feature/error-handling — Robust Error Handling & Validation

## Scope

Current error handling issues:
- `merge.py`: No validation when glob matches zero files — pandas concat raises cryptic ValueError
- `merge.py`: `Path().glob()` doesn't support absolute paths — raises NotImplementedError
- `filter.py`: `sample_name` extraction with `.split("/")[-1].split(".")[0]` is fragile
- `filter.py`: `run_filtering()` can pass `None` as `matching_km_row` to `write_filtered_line()`
- `chunk.py`: `km_options.split(" ")` fails on multi-space or quoted arguments
- `chunk.py`: No validation that output_dir exists before running
- `utils.py`: No error handling if time_it receives non-callable
- General: RuntimeError used for everything — no custom exception hierarchy
- General: No file existence checks before reading files

### Goals
- Add custom exception classes for clear error categorization
- Validate file paths exist before attempting operations
- Fix absolute path glob bug in merge
- Add proper empty-input handling across all modules
- Improve error messages to be actionable

## TODO

1. [x] Create custom exception classes (exceptions.py)
2. [x] Fix merge.py glob to handle absolute paths
3. [x] Add empty input validation to merge
4. [x] Validate file existence in filter before reading
5. [x] Fix fragile sample_name extraction in filter
6. [x] Add output_dir creation/validation in chunk
7. [x] Use shlex.split instead of str.split for km options
8. [x] Replace bare RuntimeError with custom exceptions where appropriate

## Retrospective

### What went well
- Custom exceptions make error sources immediately identifiable
- glob fix using Path(input_glob).parent.glob() handles absolute paths correctly
- shlex.split properly handles quoted arguments in km options
- File existence validation catches problems before they turn into cryptic pandas errors
- Output directory auto-creation is a nice UX improvement

### What could be improved
- Could add exit code mapping (different exit codes for different error types)
- Could add a --debug flag that shows full tracebacks vs user-friendly messages
- Exception hierarchy could be deeper (e.g., FileValidationError subtypes)

### Iteration improvements made
- Added Path-based sample_name extraction using Path.stem for robustness
- Added validation that merged output has rows before writing
