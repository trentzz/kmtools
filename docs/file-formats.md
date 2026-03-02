# File Formats

## Reference file (standard mode)

Used by `kmtools filter` in standard mode.

**Format:** Tab-separated (`.tsv`) or comma-separated (`.csv`)

**Required columns:**

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `CHROM` | string | Chromosome identifier | `chr1` |
| `POS` | integer | 1-based genomic position of the variant | `114716126` |
| `REF` | string | Reference allele | `C` |
| `ALT` | string | Alternate allele | `T` |
| `TYPE` | string | Variant type | `Substitution`, `Deletion`, `Insertion` |

**Example (`reference.tsv`):**

```
CHROM	POS	REF	ALT	TYPE
chr1	114716126	C	T	Substitution
chr3	178936091	AGG	A	Deletion
chr7	140453136	T	TA	Insertion
```

---

## Alt reference file (use-alt mode)

Used by `kmtools filter --use-alt`.

**Format:** Tab-separated (`.tsv`) or comma-separated (`.csv`)

**Required columns:**

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `CHROM` | string | Chromosome identifier | `chr1` |
| `ALT_SEQUENCE` | string | Full expected variant sequence | `TAGCTGGATT...` |
| `TYPE` | string | Variant type | `Substitution` |

**Example (`alt_reference.tsv`):**

```
CHROM	ALT_SEQUENCE	TYPE
chr1	TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT	Substitution
chr3	GCAATGCCCTCAGAATCTGTTCAGTTTGTACTTCTGATGACAGAAAGAGCCTCAGAACATCCCCAAACT	Deletion
```

---

## km find_mutation output

Produced by `km find_mutation` and consumed by `kmtools filter` and `kmtools merge`.

**Format:** Tab-separated, with optional `#` comment lines

**Required columns:**

| Column | Type | Description |
|--------|------|-------------|
| `Database` | string | Path to the jellyfish database file |
| `Query` | string | Target query in `CHROM_START_END` format |
| `Type` | string | Variant type (`Substitution`, `Deletion`, `Insertion`, `Reference`) |
| `Variant_name` | string | Variant description (`pos:ref/ALT:pos`) |
| `rVAF` | float | Relative variant allele frequency |
| `Expression` | float | Expression level of the variant |
| `Min_coverage` | integer | Minimum k-mer coverage |
| `Start_offset` | integer | Start offset of the variant |
| `Sequence` | string | Observed variant sequence |
| `Reference_expression` | float | Expression level of the reference |
| `Reference_sequence` | string | Expected reference sequence |
| `Info` | string | Status info (`vs_ref`, `reference`, `not_tested`) |

**Example:**

```
Database	Query	Type	Variant_name	rVAF	Expression	Min_coverage	Start_offset	Sequence	Reference_expression	Reference_sequence	Info
file.jf	chr1_114716091_114716161	Substitution	36:c/T:37	0.450	220.5	12	0	TAGCTG...	5033.0	TAGCTG...	vs_ref
```

---

## Filtered output (standard mode)

Produced by `kmtools filter` in standard mode.

**Format:** TSV, CSV, or XLSX (controlled by `--output-type`)

| Column | Description |
|--------|-------------|
| `SAMPLE` | Sample name (derived from km output filename stem) |
| `CHROM` | Chromosome |
| `POS` | Variant position |
| `REF` | Reference allele |
| `ALT` | Alternate allele |
| `TYPE` | Variant type |
| `FOUND` | `TRUE` if variant passed all filter conditions, `FALSE` otherwise |
| `FILTER_NOTES` | Notes about why filtering failed (populated when only the count threshold fails) |
| `KMER_VAF` | rVAF from km output (empty if not found) |
| `KMER_MIN_COVERAGE` | Min_coverage from km output (empty if not found) |
| `KMER_EXPRESSION` | Expression from km output (empty if not found) |
| `REF_SEQUENCE` | Reference_sequence from km output (empty if not found) |
| `VARIANT_SEQUENCE` | Sequence from km output (empty if not found) |

---

## Filtered output (use-alt mode)

Produced by `kmtools filter --use-alt`.

**Format:** TSV, CSV, or XLSX (controlled by `--output-type`)

| Column | Description |
|--------|-------------|
| `SAMPLE` | Sample name (derived from km output filename stem) |
| `CHROM` | Chromosome |
| `ALT_SEQUENCE` | The expected alternate sequence from the alt reference |
| `TYPE` | Variant type |
| `FOUND` | `TRUE` if variant passed all filter conditions, `FALSE` otherwise |
| `FILTER_NOTES` | Notes about why filtering failed |
| `KMER_VAF` | rVAF from km output (empty if not found) |
| `KMER_MIN_COVERAGE` | Min_coverage from km output (empty if not found) |
| `KMER_EXPRESSION` | Expression from km output (empty if not found) |
| `REF_SEQUENCE` | Reference_sequence from km output (empty if not found) |
| `VARIANT_SEQUENCE` | Sequence from km output (empty if not found) |
