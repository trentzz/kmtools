# Use-Alt Mode

The `--use-alt` flag enables an alternate filtering mode where km sequences are matched directly against a provided `ALT_SEQUENCE` instead of computing the expected variant sequence from genomic coordinates.

## Motivation

In standard mode, the filter computes the expected variant sequence from:
- The reference sequence (`Reference_sequence` from km output)
- The variant position (`POS`) and alleles (`REF`, `ALT`) from the reference file
- The query coordinates to calculate the offset

This works well when you have precise genomic coordinates. However, in some cases:
- You may not have exact POS/REF/ALT information
- You may want to match against a known alternate sequence directly
- The variant representation may differ between your reference and km's output

In these cases, `--use-alt` provides a simpler matching approach: compare the km `Sequence` directly against a provided `ALT_SEQUENCE`.

## Alt reference file format

The alt reference file must be a `.tsv` or `.csv` with these required columns:

| Column | Description |
|--------|-------------|
| `CHROM` | Chromosome identifier (must match the chromosome in the km `Query` field) |
| `ALT_SEQUENCE` | The full expected variant sequence to match against km `Sequence` |
| `TYPE` | Variant type (e.g., `Substitution`, `Deletion`, `Insertion`) |

Example `alt_reference.tsv`:

```
CHROM	ALT_SEQUENCE	TYPE
chr1	TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT	Substitution
chr3	GCAATGCCCTCAGAATCTGTTCAGTTTGTACTTCTGATGACAGAAAGAGCCTCAGAACATCCCCAAACT	Deletion
chr7	TGCTGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCT	Insertion
```

## Matching logic

For each row in the alt reference file, the filter:

1. Iterates through all km output rows
2. Filters to rows where the chromosome (from `Query`) matches `CHROM`
3. Checks all filter conditions:
   - **TYPE match**: km `Type` == alt reference `TYPE`
   - **COUNT threshold**: km `Min_coverage` >= `--count-threshold`
   - **Sequence match**: km `Sequence` (uppercased) == `ALT_SEQUENCE` (uppercased)
   - **vs_ref info**: km `Info` == `"vs_ref"`
4. Reports the first passing row, or the best partial match

## Output columns

Use-alt mode produces different output columns than standard mode:

| Column | Description |
|--------|-------------|
| `SAMPLE` | Sample name (derived from km output filename) |
| `CHROM` | Chromosome from alt reference |
| `ALT_SEQUENCE` | The expected alternate sequence |
| `TYPE` | Variant type |
| `FOUND` | `TRUE` if all conditions passed, `FALSE` otherwise |
| `FILTER_NOTES` | Notes about why filtering failed (e.g., low count) |
| `KMER_VAF` | Variant allele frequency (if found) |
| `KMER_MIN_COVERAGE` | Minimum k-mer coverage (if found) |
| `KMER_EXPRESSION` | Expression level (if found) |
| `REF_SEQUENCE` | Reference sequence from km (if found) |
| `VARIANT_SEQUENCE` | Variant sequence from km (if found) |

## Example workflow

### 1. Prepare an alt reference file

Create a TSV with the sequences you expect to find:

```bash
cat > alt_reference.tsv << 'EOF'
CHROM	ALT_SEQUENCE	TYPE
chr1	TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT	Substitution
chr3	GCAATGCCCTCAGAATCTGTTCAGTTTGTACTTCTGATGACAGAAAGAGCCTCAGAACATCCCCAAACT	Deletion
EOF
```

### 2. Run the filter with --use-alt

```bash
kmtools filter \
    --reference alt_reference.tsv \
    --km-output km_results.txt \
    --output filtered_alt.tsv \
    --use-alt \
    --count-threshold 2 \
    --verbose
```

### 3. Inspect results

```bash
# View found variants
head -1 filtered_alt.tsv && grep TRUE filtered_alt.tsv
```

### 4. Generate plots (optional)

```bash
kmtools plot filtered_alt.tsv --output-dir plots --charts overall,type
```

## When to use standard mode vs use-alt mode

| Scenario | Mode |
|----------|------|
| You have POS/REF/ALT coordinates | Standard mode |
| You have known variant sequences | `--use-alt` mode |
| You want position-based verification | Standard mode |
| You want direct sequence comparison | `--use-alt` mode |
