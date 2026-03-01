# feature/merge-improvements — Robust Merge Module

## Scope

Current merge.py problems:
- `Path().glob()` doesn't support absolute paths — raises NotImplementedError
- No validation when zero files match glob patterns — pandas concat raises cryptic ValueError
- Files are globbed twice (once for reading, once for deletion) — race condition risk
- No duplicate row detection or handling
- No logging of how many files/rows were merged
- No header validation across files (mismatched columns could silently produce garbage)
- Deletion loop re-globs instead of tracking resolved files

### Goals
- Fix absolute path glob using stdlib glob module
- Validate inputs before processing (no empty file list)
- Resolve files once and reuse the list for read and delete
- Add column consistency validation across input files
- Add summary logging (file count, row count)
- Add duplicate row dropping option

## TODO

1. [x] Fix glob to handle absolute paths using stdlib glob
2. [x] Add empty input validation with clear error message
3. [x] Resolve file list once and reuse for both read and delete
4. [x] Add column header validation across input files
5. [x] Add --sort option to CLI for sorting merged output
6. [x] Add --drop-duplicates option to CLI
7. [x] Add summary logging (file count, row count)
8. [x] Wire new options through CLI parser

## Retrospective

### What went well
- Fixing the glob issue is a clear, impactful bug fix
- Single file resolution eliminates the race condition
- Column validation catches mismatched inputs early with an actionable error
- Summary logging gives users confidence the merge worked correctly
- Sort and deduplicate options are practical for real workflows

### What could be improved
- Could add a --validate-only dry-run mode
- Could support different deduplication strategies (keep first, keep last, keep none)
- Column validation could be a warning instead of error in some cases

### Iteration improvements made
- Added sort_by parameter to allow sorting by specific column(s)
- Made column validation warn but continue (columns from first file used as reference)
