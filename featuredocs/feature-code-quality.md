# feature/code-quality — Type Hints, Refactoring & Consistency

## Scope

Current code quality issues:
- `chunk.py`: Unused imports (math, os), dead comment with hardcoded path on line 73, no type hints on __init__ params or run_km
- `filter.py`: Unused import (csv), no type hints on most methods, `get_ref_alt_pos_from_variant` returns a tuple on error but KmVariant on success (inconsistent)
- `merge.py`: No type hints on __init__ params
- `plot.py`: No type hints on __init__ params
- `utils.py`: `time_it` is redundant with verbose check (checks verbose twice)
- `filter_types.py`: No __post_init__ validation
- `kmtools.py`: No type hints on run_* functions
- General: Inconsistent use of str vs Path for file paths

### Goals
- Add type hints to all public method signatures
- Remove unused imports and dead code
- Fix inconsistent return types
- Standardize Path usage where appropriate
- Clean up redundant patterns

## TODO

1. [x] Remove unused imports (math, os from chunk; csv from filter)
2. [x] Remove dead hardcoded path comment from chunk.py
3. [x] Add type hints to all __init__ methods and public methods
4. [x] Fix get_ref_alt_pos_from_variant return type inconsistency
5. [x] Fix time_it redundant verbose check
6. [x] Add return type hints to run methods
7. [x] Remove redundant hasattr check in chunk.run()

## Retrospective

### What went well
- Type hints make the codebase significantly more readable and IDE-friendly
- Removing dead code (hardcoded path comment, unused imports) reduced noise
- Fixing the return type inconsistency in get_ref_alt_pos_from_variant prevents subtle bugs
- Consistent Optional[KmVariant] return is cleaner than returning a tuple of Nones

### What could be improved
- Could add py.typed marker file for PEP 561 compliance
- Could add a ruff or mypy configuration to enforce type hints going forward
- Some internal methods could benefit from type hints too

### Iteration improvements made
- Standardized Path usage in merge and plot constructors
- Cleaned up time_it to avoid double-checking verbose
