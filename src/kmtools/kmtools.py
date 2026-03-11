#!/usr/bin/env python3
"""
kmtools: A toolkit for running, filtering, merging, and plotting km find_mutation results.
"""

import argparse
import sys
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path

from kmtools.chunk import Chunk
from kmtools.filter import Filter
from kmtools.merge import Merge
from kmtools.plot import Plot
from kmtools.utils import Utils


def _get_version() -> str:
    try:
        return version("kmtools-bio")
    except PackageNotFoundError:
        return "dev"


# === Subcommand Functions ===


def run_chunk(args: argparse.Namespace) -> None:
    Utils.log(f"Running chunk with args: {args}", args.verbose)
    Chunk(
        threads=args.threads,
        km_find_mutation_options=args.km_find_mutation_options,
        km_target_directory=args.km_target_directory,
        km_jellyfish_file=args.km_jellyfish_file,
        output_dir=args.output_dir,
        prefix=args.prefix,
        merge=args.merge,
        merge_output=args.merge_output,
        merge_keep=args.merge_keep,
        verbose=args.verbose,
    ).run()


def run_filter(args: argparse.Namespace) -> None:
    Utils.log(f"Running filter with args: {args}", args.verbose)

    # Validate input files exist
    for path, label in [(args.reference, "Reference file"), (args.km_output, "KM output file")]:
        if not Path(path).exists():
            print(f"Error: {label} not found: {path}", file=sys.stderr)
            sys.exit(1)

    Filter(
        reference=args.reference,
        km_output=args.km_output,
        output=args.output,
        output_type=args.output_type,
        count_threshold=args.count_threshold,
        use_alt=getattr(args, "use_alt", False),
        verbose=args.verbose,
    ).run()


def run_merge(args: argparse.Namespace) -> None:
    Utils.log(f"Running merge with args: {args}", args.verbose)
    Merge(
        inputs=args.inputs,
        output=args.output,
        keep=args.keep,
        sort_by=getattr(args, "sort_by", None),
        drop_duplicates=getattr(args, "drop_duplicates", False),
        verbose=args.verbose,
    ).run()


def run_plot(args: argparse.Namespace) -> None:
    Utils.log(f"Running plot with args: {args}", args.verbose)

    if not Path(args.file).exists():
        print(f"Error: Input file not found: {args.file}", file=sys.stderr)
        sys.exit(1)

    Plot(
        file=args.file,
        output_dir=args.output_dir,
        charts=args.charts,
        verbose=args.verbose,
    ).run()


def run_all(args: argparse.Namespace) -> None:
    """Runs the entire pipeline: chunk -> merge -> filter -> plot."""
    Utils.log(
        "Running full kmtools pipeline (chunk -> merge -> filter -> plot)", args.verbose
    )

    # --- Step 1: Chunk ---
    chunk = Chunk(
        threads=args.threads,
        km_find_mutation_options=args.km_find_mutation_options,
        km_target_directory=args.km_target_directory,
        km_jellyfish_file=args.km_jellyfish_file,
        output_dir=args.output_dir,
        merge=True,
        merge_output=args.merge_output,
        verbose=args.verbose,
    )
    chunk.run()

    # --- Step 2: Filter ---
    filter_step = Filter(
        reference=args.reference,
        km_output=args.merge_output,
        output=args.filtered_output,
        output_type=args.output_type,
        count_threshold=args.count_threshold,
        verbose=args.verbose,
    )
    filter_step.run()

    # --- Step 3: Plot ---
    plot = Plot(
        file=args.filtered_output,
        output_dir=args.output_dir,
        charts=args.charts,
        verbose=args.verbose,
    )
    plot.run()

    Utils.log("Pipeline complete.", args.verbose)


# === CLI Setup ===


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="kmtools",
        description="kmtools: run, filter, merge, and plot km find_mutation results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {_get_version()}"
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    subparsers = parser.add_subparsers(
        title="subcommands", dest="command", required=True
    )

    # --- chunk ---
    chunk_parser = subparsers.add_parser(
        "chunk",
        help="Run km find_mutation in parallel across chunked target sequences",
        description="Split target sequences into chunks and run km find_mutation in parallel using multiple threads.",
        epilog="Example:\n  kmtools chunk --threads 8 --km-find-mutation-options \"--ratio 0.0001\" "
               "--km-target-directory targets_split --km-jellyfish-file db.jf --merge",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    chunk_parser.add_argument(
        "--threads", type=int, required=True, help="Number of parallel km processes"
    )
    chunk_parser.add_argument(
        "--km-find-mutation-options",
        type=str,
        required=True,
        help="Options passed directly to km find_mutation (e.g., \"--ratio 0.0001\")",
    )
    chunk_parser.add_argument(
        "--km-target-directory",
        type=str,
        required=True,
        help="Directory containing pre-split target subdirectories (one per thread)",
    )
    chunk_parser.add_argument(
        "--km-jellyfish-file",
        type=str,
        required=True,
        help="Path to the jellyfish k-mer counts file (.jf)",
    )
    chunk_parser.add_argument(
        "--output-dir",
        type=str,
        default=".",
        help="Directory to save output files (default: current directory)",
    )
    chunk_parser.add_argument(
        "--merge", action="store_true", help="Automatically merge chunk outputs after processing"
    )
    chunk_parser.add_argument(
        "--prefix",
        type=str,
        default="km_find_mutation_output",
        help="Prefix for output files (default: km_find_mutation_output)",
    )
    chunk_parser.add_argument(
        "--merge-output",
        type=str,
        default="km_find_mutation_merged_output.txt",
        help="Filename for merged output when --merge is used (default: km_find_mutation_merged_output.txt)",
    )
    chunk_parser.add_argument(
        "--merge-keep",
        action="store_true",
        help="Keep intermediate chunk files after merging",
    )
    chunk_parser.set_defaults(func=run_chunk)

    # --- merge ---
    merge_parser = subparsers.add_parser(
        "merge",
        help="Combine multiple km output files into a single result",
        description="Merge multiple km find_mutation output files (TSV) into one combined file.",
        epilog="Example:\n  kmtools merge chunk_*.txt --output merged.txt --keep",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    merge_parser.add_argument(
        "inputs", nargs="+", type=str, help="Input files or glob patterns to merge"
    )
    merge_parser.add_argument(
        "--output", type=str, required=True, help="Output file for merged results"
    )
    merge_parser.add_argument(
        "--keep", action="store_true", help="Keep input files after merging (default: delete them)"
    )
    merge_parser.add_argument(
        "--sort-by",
        type=str,
        default=None,
        help="Sort merged output by this column name",
    )
    merge_parser.add_argument(
        "--drop-duplicates",
        action="store_true",
        help="Remove duplicate rows from merged output",
    )
    merge_parser.set_defaults(func=run_merge)

    # --- filter ---
    filter_parser = subparsers.add_parser(
        "filter",
        help="Filter km results against a reference to identify true variants",
        description="Filter km find_mutation results by comparing against a reference file with known variants.",
        epilog="Example:\n  kmtools filter --reference variants.tsv --km-output merged.txt --output filtered.tsv",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    filter_parser.add_argument(
        "--reference", type=str, required=True,
        help="Reference file with known variants (.csv or .tsv with CHROM, POS, REF, ALT, TYPE columns)",
    )
    filter_parser.add_argument(
        "--km-output", type=str, required=True,
        help="km find_mutation output file to filter",
    )
    filter_parser.add_argument(
        "--output", type=str, required=True, help="Output file for filtered results"
    )
    filter_parser.add_argument(
        "--output-type",
        type=str,
        choices=["tsv", "csv", "xlsx"],
        default="tsv",
        help="Output file format (default: tsv)",
    )
    filter_parser.add_argument(
        "--count-threshold",
        type=int,
        default=2,
        help="Minimum k-mer count threshold for a variant to pass (default: 2)",
    )
    filter_parser.add_argument(
        "--use-alt",
        action="store_true",
        help="Use alternate sequence mode: reference file should have CHROM, ALT_SEQUENCE, TYPE columns "
             "instead of CHROM, POS, REF, ALT, TYPE. Matches km sequences directly against provided ALT_SEQUENCE.",
    )
    filter_parser.set_defaults(func=run_filter)

    # --- plot ---
    plot_parser = subparsers.add_parser(
        "plot",
        help="Generate charts from filtered results",
        description="Generate visualizations (VAF distribution, variant types, sample summary) from filtered km results.",
        epilog="Example:\n  kmtools plot filtered.tsv --output-dir plots --charts vaf,type,sample,overall",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    plot_parser.add_argument("file", type=str, help="Filtered results file to plot")
    plot_parser.add_argument(
        "--output-dir", type=str, default=".", help="Directory to save plot images (default: current directory)"
    )
    plot_parser.add_argument(
        "--charts",
        type=str,
        default="all",
        help="Comma-separated chart types: vaf, type, sample, overall, or all (default: all)",
    )
    plot_parser.set_defaults(func=run_plot)

    # --- runall ---
    runall_parser = subparsers.add_parser(
        "runall",
        help="Run the full pipeline: chunk -> filter -> plot",
        description="Execute the complete kmtools workflow in one command: chunk target sequences, "
                    "run km find_mutation in parallel, merge results, filter against reference, and generate plots.",
        epilog="Example:\n  kmtools runall --threads 8 --km-find-mutation-options \"--ratio 0.0001\" "
               "--km-target-directory targets_split --km-jellyfish-file db.jf "
               "--merge-output merged.txt --reference variants.tsv --filtered-output filtered.tsv",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Chunk step
    runall_parser.add_argument(
        "--threads", type=int, required=True, help="Number of parallel threads for chunking"
    )
    runall_parser.add_argument(
        "--km-find-mutation-options",
        type=str,
        required=True,
        help="Options passed to km find_mutation",
    )
    runall_parser.add_argument(
        "--km-target-directory",
        type=str,
        required=True,
        help="Directory containing pre-split target subdirectories",
    )
    runall_parser.add_argument(
        "--km-jellyfish-file",
        type=str,
        required=True,
        help="Path to the jellyfish k-mer counts file (.jf)",
    )

    # Merge step
    runall_parser.add_argument(
        "--merge-output", type=str, required=True,
        help="Output file for merged chunk results",
    )

    # Filter step
    runall_parser.add_argument(
        "--reference", type=str, required=True,
        help="Reference file with known variants (.csv or .tsv)",
    )
    runall_parser.add_argument(
        "--filtered-output",
        type=str,
        required=True,
        help="Output file for filtered results",
    )
    runall_parser.add_argument(
        "--output-type",
        type=str,
        choices=["tsv", "csv", "xlsx"],
        default="tsv",
        help="Output format for filtered results (default: tsv)",
    )
    runall_parser.add_argument(
        "--count-threshold",
        type=int,
        default=2,
        help="Minimum k-mer count threshold (default: 2)",
    )

    # Plot step
    runall_parser.add_argument(
        "--output-dir", type=str, default=".",
        help="Directory to save plots (default: current directory)",
    )
    runall_parser.add_argument(
        "--charts",
        type=str,
        default="all",
        help="Comma-separated chart types to generate (default: all)",
    )

    runall_parser.set_defaults(func=run_all)

    # --- parse and dispatch ---
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
