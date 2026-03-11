#!/usr/bin/env python3
"""
Benchmarks for kmtools merge and filter operations.

Usage:
    python benchmarks/run_benchmarks.py
"""

import random
import sys
import time
from pathlib import Path

# Add src to path for local development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pandas as pd
from kmtools.merge import Merge
from kmtools.filter import Filter

random.seed(42)

KM_OUTPUT_HEADER = (
    "Database\tQuery\tType\tVariant_name\trVAF\tExpression\t"
    "Min_coverage\tStart_offset\tSequence\tReference_expression\t"
    "Reference_sequence\tInfo"
)

REF_SEQ = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"
VAR_SEQ = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"


def make_km_row(idx: int) -> str:
    """Generate a km output row with slight variation."""
    coverage = random.randint(1, 100)
    vaf = random.uniform(0.0, 1.0)
    expression = random.uniform(0.0, 500.0)
    info = "vs_ref" if random.random() > 0.2 else "reference"
    vtype = random.choice(["Substitution", "Deletion", "Insertion"])
    return (
        f"file.jf\tchr1_114716091_114716161\t{vtype}\t36:c/T:37\t{vaf:.3f}\t{expression:.1f}\t"
        f"{coverage}\t0\t{VAR_SEQ}\t5033.0\t{REF_SEQ}\t{info}"
    )


def write_km_file(path: Path, num_rows: int) -> None:
    """Write a km output file with the given number of rows."""
    lines = [KM_OUTPUT_HEADER]
    for i in range(num_rows):
        lines.append(make_km_row(i))
    path.write_text("\n".join(lines) + "\n")


def write_reference_file(path: Path, num_variants: int) -> None:
    """Write a reference TSV with the given number of variants."""
    lines = ["CHROM\tPOS\tREF\tALT\tTYPE"]
    for i in range(num_variants):
        lines.append(f"chr1\t114716126\tC\tT\tSubstitution")
    path.write_text("\n".join(lines) + "\n")


def benchmark_merge(tmp_dir: Path, num_files: int, rows_per_file: int) -> float:
    """Benchmark merge operation, return elapsed seconds."""
    merge_dir = tmp_dir / f"merge_{num_files}x{rows_per_file}"
    merge_dir.mkdir(parents=True, exist_ok=True)

    for i in range(num_files):
        write_km_file(merge_dir / f"chunk_{i:04d}.txt", rows_per_file)

    output = merge_dir / "merged.txt"
    inputs = [str(merge_dir / f"chunk_{i:04d}.txt") for i in range(num_files)]

    start = time.perf_counter()
    Merge(inputs=inputs, output=str(output), keep=True).run()
    elapsed = time.perf_counter() - start

    # Cleanup
    for f in merge_dir.iterdir():
        f.unlink()
    merge_dir.rmdir()

    return elapsed


def benchmark_merge_with_options(tmp_dir: Path, num_files: int, rows_per_file: int) -> float:
    """Benchmark merge with sort-by and drop-duplicates."""
    merge_dir = tmp_dir / f"merge_opts_{num_files}x{rows_per_file}"
    merge_dir.mkdir(parents=True, exist_ok=True)

    for i in range(num_files):
        write_km_file(merge_dir / f"chunk_{i:04d}.txt", rows_per_file)

    output = merge_dir / "merged.txt"
    inputs = [str(merge_dir / f"chunk_{i:04d}.txt") for i in range(num_files)]

    start = time.perf_counter()
    Merge(inputs=inputs, output=str(output), keep=True, sort_by="Min_coverage", drop_duplicates=True).run()
    elapsed = time.perf_counter() - start

    for f in merge_dir.iterdir():
        f.unlink()
    merge_dir.rmdir()

    return elapsed


def benchmark_filter(tmp_dir: Path, num_km_rows: int, num_ref_variants: int) -> float:
    """Benchmark filter operation, return elapsed seconds."""
    filter_dir = tmp_dir / f"filter_{num_km_rows}x{num_ref_variants}"
    filter_dir.mkdir(parents=True, exist_ok=True)

    ref_file = filter_dir / "reference.tsv"
    km_file = filter_dir / "km_output.txt"
    out_file = filter_dir / "filtered.tsv"

    write_reference_file(ref_file, num_ref_variants)
    write_km_file(km_file, num_km_rows)

    start = time.perf_counter()
    Filter(
        reference=str(ref_file),
        km_output=str(km_file),
        output=str(out_file),
        output_type="tsv",
        count_threshold=2,
    ).run()
    elapsed = time.perf_counter() - start

    for f in filter_dir.iterdir():
        f.unlink()
    filter_dir.rmdir()

    return elapsed


def main():
    import tempfile
    tmp_dir = Path(tempfile.mkdtemp(prefix="kmtools_bench_"))

    print("=" * 70)
    print("kmtools Benchmark Results")
    print("=" * 70)

    # --- Merge benchmarks ---
    print("\n## Merge Benchmarks")
    print(f"{'Files':>8} {'Rows/File':>10} {'Total Rows':>12} {'Time (s)':>10} {'Rows/s':>12}")
    print("-" * 56)

    merge_configs = [
        (10, 10),
        (50, 10),
        (100, 10),
        (500, 10),
    ]

    merge_results = []
    for num_files, rows_per_file in merge_configs:
        elapsed = benchmark_merge(tmp_dir, num_files, rows_per_file)
        total_rows = num_files * rows_per_file
        rows_per_sec = total_rows / elapsed if elapsed > 0 else float('inf')
        print(f"{num_files:>8} {rows_per_file:>10} {total_rows:>12} {elapsed:>10.4f} {rows_per_sec:>12.0f}")
        merge_results.append((num_files, rows_per_file, total_rows, elapsed, rows_per_sec))

    # --- Merge with options ---
    print("\n## Merge with sort-by + drop-duplicates")
    print(f"{'Files':>8} {'Rows/File':>10} {'Total Rows':>12} {'Time (s)':>10} {'Rows/s':>12}")
    print("-" * 56)

    merge_opt_results = []
    for num_files, rows_per_file in merge_configs:
        elapsed = benchmark_merge_with_options(tmp_dir, num_files, rows_per_file)
        total_rows = num_files * rows_per_file
        rows_per_sec = total_rows / elapsed if elapsed > 0 else float('inf')
        print(f"{num_files:>8} {rows_per_file:>10} {total_rows:>12} {elapsed:>10.4f} {rows_per_sec:>12.0f}")
        merge_opt_results.append((num_files, rows_per_file, total_rows, elapsed, rows_per_sec))

    # --- Filter benchmarks ---
    print("\n## Filter Benchmarks")
    print(f"{'KM Rows':>8} {'Ref Vars':>10} {'Time (s)':>10} {'Rows/s':>12}")
    print("-" * 44)

    filter_configs = [
        (10, 5),
        (50, 5),
        (100, 5),
        (500, 5),
        (1000, 5),
    ]

    filter_results = []
    for num_km_rows, num_ref_variants in filter_configs:
        elapsed = benchmark_filter(tmp_dir, num_km_rows, num_ref_variants)
        rows_per_sec = num_km_rows / elapsed if elapsed > 0 else float('inf')
        print(f"{num_km_rows:>8} {num_ref_variants:>10} {elapsed:>10.4f} {rows_per_sec:>12.0f}")
        filter_results.append((num_km_rows, num_ref_variants, elapsed, rows_per_sec))

    # Cleanup
    try:
        tmp_dir.rmdir()
    except OSError:
        pass

    print("\n" + "=" * 70)
    print("Benchmark complete.")

    # Output markdown table for docs
    print("\n\n### Markdown tables for docs/testing-and-benchmarking.md:\n")

    print("#### Merge Performance\n")
    print("| Files | Rows/File | Total Rows | Time (s) | Throughput (rows/s) |")
    print("|------:|----------:|-----------:|---------:|--------------------:|")
    for f, rpf, total, t, rps in merge_results:
        print(f"| {f} | {rpf} | {total} | {t:.4f} | {rps:,.0f} |")

    print("\n#### Merge with sort-by + drop-duplicates\n")
    print("| Files | Rows/File | Total Rows | Time (s) | Throughput (rows/s) |")
    print("|------:|----------:|-----------:|---------:|--------------------:|")
    for f, rpf, total, t, rps in merge_opt_results:
        print(f"| {f} | {rpf} | {total} | {t:.4f} | {rps:,.0f} |")

    print("\n#### Filter Performance\n")
    print("| KM Rows | Ref Variants | Time (s) | Throughput (rows/s) |")
    print("|--------:|-------------:|---------:|--------------------:|")
    for kr, rv, t, rps in filter_results:
        print(f"| {kr} | {rv} | {t:.4f} | {rps:,.0f} |")


if __name__ == "__main__":
    main()
