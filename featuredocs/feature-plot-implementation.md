# feature/plot-implementation — Real Matplotlib Charts

## Scope

The plot module is currently a stub that just prints chart names. It needs real implementations for:
- **VAF distribution**: Histogram of variant allele frequencies from filtered results
- **Per-sample summary**: Bar chart showing FOUND vs NOT FOUND counts per sample
- **Per-type summary**: Bar chart showing variant type distribution
- **Overall summary**: Pie chart of overall found/not-found ratio

The filtered output has columns: SAMPLE, CHROM, POS, REF, ALT, TYPE, FOUND, FILTER_NOTES, KMER_VAF, KMER_MIN_COVERAGE, KMER_EXPRESSION, REF_SEQUENCE, VARIANT_SEQUENCE.

### Goals
- Implement 4 chart types using matplotlib
- Read filtered output data using pandas
- Save charts as PNG files to output directory
- Handle edge cases (empty data, missing columns)
- Use non-interactive matplotlib backend for headless environments

## TODO

1. [x] Set matplotlib to Agg backend for headless environments
2. [x] Implement data loading and validation
3. [x] Implement VAF distribution histogram
4. [x] Implement per-sample found/not-found bar chart
5. [x] Implement variant type distribution bar chart
6. [x] Implement overall summary pie chart
7. [x] Handle edge cases (empty data, no VAF values)
8. [x] Wire output filenames correctly

## Retrospective

### What went well
- Using Agg backend ensures it works on servers without a display
- Each chart method is self-contained and can be run independently
- pandas makes data loading and grouping straightforward
- Edge case handling prevents crashes on empty or partial data

### What could be improved
- Could add more chart customization options (colors, figure size, DPI)
- Could support additional output formats (SVG, PDF)
- Could add interactive HTML charts using plotly as an optional dependency
- No multi-sample comparison charts yet

### Iteration improvements made
- Added numeric conversion with error handling for VAF column
- Added output directory auto-creation
- Made chart methods return the file path for logging
