# feature/cli-ux — CLI User Experience Improvements

## Scope

Current CLI issues:
- No --version flag
- Help text is minimal and inconsistent across subcommands
- No file path validation before running (user gets cryptic pandas errors)
- No progress feedback during long operations
- runall has no --output-type and --count-threshold args (uses hardcoded defaults)
- runall references args.keep and args.output which are never defined on the parser
- chunk parser help text says "Chunk input data" — not descriptive
- No epilog examples in help text

### Goals
- Add --version flag using version from pyproject.toml
- Improve help descriptions for all subcommands
- Add file path validation (check files exist before running)
- Fix runall parser missing arguments
- Add epilog examples to subcommand help

## TODO

1. [x] Add --version flag to main parser
2. [x] Improve subcommand help descriptions
3. [x] Fix runall parser missing args (output_type, count_threshold, keep, output, km_output)
4. [x] Add file existence validation in run_* functions
5. [x] Add epilog examples to key subcommands
6. [x] Verify runall pipeline wiring is correct

## Retrospective

### What went well
- --version flag using importlib.metadata is the standard approach
- File validation before running prevents cryptic downstream errors
- Fixed several missing args in runall that would have caused AttributeError at runtime
- Epilog examples give users quick copy-paste commands

### What could be improved
- Could add shell completion support (argcomplete)
- Could add color output for errors/warnings
- Could add a --dry-run mode that shows what would be done

### Iteration improvements made
- Fixed run_all to properly chain outputs between pipeline stages
- Added missing --km-output, --output, --output-type, --count-threshold to runall
