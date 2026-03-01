# feature/use-alt — Direct Alternate Sequence Filtering

## Scope

### Current behaviour (reference mode)
The filter currently works by:
1. Reading a reference file with columns: CHROM, POS, REF, ALT, TYPE
2. For each reference variant, computing the expected full alternate sequence by
   substituting ALT into the km Reference_sequence at position POS
3. Comparing km's reported Sequence against this calculated sequence
4. If they match (plus TYPE, COUNT, INFO conditions), marking as FOUND

### Problem
This requires the user to know the exact genomic position, reference allele, and
alternate allele for each variant. In some workflows, the user already knows the
full expected alternate sequence and just wants to check whether km detected it
— without needing to decompose it into POS/REF/ALT components.

### --use-alt mode
When `--use-alt` is set:
1. The reference file uses a different format with columns: CHROM, ALT_SEQUENCE, TYPE
   - ALT_SEQUENCE is the full expected variant sequence (what you'd see in km's Sequence column)
   - POS, REF, ALT columns are not required
2. The filter skips `get_calculated_reference_alt()` — the ALT_SEQUENCE is used directly
3. Matching: for each alt row, scan km output for rows where:
   - km Sequence (uppercased) matches ALT_SEQUENCE (uppercased)
   - TYPE matches
   - Min_coverage >= count_threshold
   - Info == "vs_ref"
4. Output format adapts: POS/REF/ALT columns show "N/A" since they're not applicable,
   ALT_SEQUENCE is included instead

### Goals
- Add `--use-alt` flag to filter CLI subcommand
- Add `use_alt` parameter to Filter class
- Implement alternate verification for the new reference format
- Implement alternate filtering logic that uses ALT_SEQUENCE directly
- Write adapted output that makes sense for this mode
- Wire through to runall as well
- Tests for all new paths

## TODO

1. [x] Add `use_alt` parameter to Filter.__init__
2. [x] Add `verify_alt_reference()` for the new column format (CHROM, ALT_SEQUENCE, TYPE)
3. [x] Add `run_alt_filtering()` method that compares directly against ALT_SEQUENCE
4. [x] Adapt `write_filtered_line` for use-alt mode (no POS/REF/ALT)
5. [x] Modify `run()` to dispatch to the correct filtering path
6. [x] Add `--use-alt` flag to filter CLI parser
7. [x] Wire `--use-alt` through to runall parser
8. [x] Create test alt reference file
9. [x] Write unit tests for verify_alt_reference
10. [x] Write unit tests for run_alt_filtering
11. [x] Write unit tests for filter_alt_line
12. [x] Write CLI parsing test for --use-alt flag
13. [x] Run all tests

## Retrospective

### What went well
- Clean separation between reference mode and alt mode via use_alt flag
- Reuses most of the existing filter infrastructure (write_output, filter conditions)
- ALT_SEQUENCE comparison is simpler and more direct than the reference calculation
- Test coverage is thorough — tests both modes independently
- Output format adapts sensibly (ALT_SEQUENCE column replaces POS/REF/ALT)

### What could be improved
- Could support a hybrid mode where both ALT_SEQUENCE and POS/REF/ALT are present
- Could add fuzzy matching (e.g., allow N mismatches in sequence comparison)
- Could allow matching against a substring of the km Sequence (partial match mode)
- The two filtering paths share some logic — could be unified with a strategy pattern

### Iteration improvements made
- Added CHROM matching to alt filtering (only compare against km rows for the same chromosome)
- Made output include ALT_SEQUENCE column for traceability
- Added clear error message when wrong reference format is used with --use-alt
