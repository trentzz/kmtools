from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from kmtools.utils import Utils


class Plot:
    def __init__(
        self,
        file: str,
        output_dir: str = ".",
        charts: str = "all",
        verbose: bool = False,
    ) -> None:
        self.file = Path(file)
        self.output_dir = Path(output_dir)
        self.charts = (
            charts.split(",")
            if charts != "all"
            else ["vaf", "type", "sample", "overall"]
        )
        self.verbose = verbose
        self.df: Optional[pd.DataFrame] = None

    def run(self) -> None:
        Utils.log(f"Generating plots from {self.file} -> {self.output_dir}", self.verbose)

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._load_data()

        if self.df is None or self.df.empty:
            print("No data to plot.")
            return

        chart_methods = {
            "vaf": self._plot_vaf_distribution,
            "type": self._plot_type_distribution,
            "sample": self._plot_sample_summary,
            "overall": self._plot_overall_summary,
        }

        for chart in self.charts:
            method = chart_methods.get(chart.strip())
            if method is None:
                print(f"Unknown chart type: {chart}")
                continue
            output_path = method()
            if output_path:
                print(f"Generated: {output_path}")

        print("Plotting complete.")

    def _load_data(self) -> None:
        """Load filtered results file."""
        if not self.file.exists():
            raise FileNotFoundError(f"Input file not found: {self.file}")

        if str(self.file).endswith(".csv"):
            self.df = pd.read_csv(self.file, dtype={"FOUND": str})
        else:
            self.df = pd.read_csv(self.file, sep="\t", dtype={"FOUND": str})

        Utils.log(f"Loaded {len(self.df)} rows from {self.file}", self.verbose)

    def _plot_vaf_distribution(self) -> Optional[Path]:
        """Histogram of variant allele frequencies for found variants."""
        if "KMER_VAF" not in self.df.columns:
            Utils.log("Skipping VAF plot: KMER_VAF column not found", self.verbose)
            return None

        found_df = self.df[self.df["FOUND"] == "TRUE"].copy()
        if found_df.empty:
            Utils.log("Skipping VAF plot: no FOUND=TRUE rows", self.verbose)
            return None

        found_df["KMER_VAF"] = pd.to_numeric(found_df["KMER_VAF"], errors="coerce")
        vaf_values = found_df["KMER_VAF"].dropna()

        if vaf_values.empty:
            Utils.log("Skipping VAF plot: no numeric VAF values", self.verbose)
            return None

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(vaf_values, bins=30, color="#4C72B0", edgecolor="white", alpha=0.85)
        ax.set_xlabel("Variant Allele Frequency (rVAF)")
        ax.set_ylabel("Count")
        ax.set_title("VAF Distribution (Found Variants)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()

        output_path = self.output_dir / "vaf_distribution.png"
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        return output_path

    def _plot_type_distribution(self) -> Optional[Path]:
        """Bar chart of variant type distribution."""
        if "TYPE" not in self.df.columns:
            Utils.log("Skipping type plot: TYPE column not found", self.verbose)
            return None

        type_counts = self.df["TYPE"].value_counts()
        if type_counts.empty:
            return None

        fig, ax = plt.subplots(figsize=(8, 5))
        colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974"]
        type_counts.plot(kind="bar", ax=ax, color=colors[: len(type_counts)], edgecolor="white")
        ax.set_xlabel("Variant Type")
        ax.set_ylabel("Count")
        ax.set_title("Variant Type Distribution")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.xticks(rotation=45, ha="right")
        fig.tight_layout()

        output_path = self.output_dir / "type_distribution.png"
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        return output_path

    def _plot_sample_summary(self) -> Optional[Path]:
        """Stacked bar chart showing found vs not-found per sample."""
        if "SAMPLE" not in self.df.columns or "FOUND" not in self.df.columns:
            Utils.log("Skipping sample plot: SAMPLE or FOUND column not found", self.verbose)
            return None

        summary = self.df.groupby(["SAMPLE", "FOUND"]).size().unstack(fill_value=0)

        # Ensure both TRUE and FALSE columns exist
        for col in ["TRUE", "FALSE"]:
            if col not in summary.columns:
                summary[col] = 0

        fig, ax = plt.subplots(figsize=(max(8, len(summary) * 1.2), 5))
        summary[["TRUE", "FALSE"]].plot(
            kind="bar",
            stacked=True,
            ax=ax,
            color=["#55A868", "#C44E52"],
            edgecolor="white",
        )
        ax.set_xlabel("Sample")
        ax.set_ylabel("Variant Count")
        ax.set_title("Found vs Not Found per Sample")
        ax.legend(["Found", "Not Found"])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.xticks(rotation=45, ha="right")
        fig.tight_layout()

        output_path = self.output_dir / "sample_summary.png"
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        return output_path

    def _plot_overall_summary(self) -> Optional[Path]:
        """Pie chart of overall found vs not-found ratio."""
        if "FOUND" not in self.df.columns:
            Utils.log("Skipping overall plot: FOUND column not found", self.verbose)
            return None

        found_counts = self.df["FOUND"].value_counts()
        if found_counts.empty:
            return None

        labels = []
        sizes = []
        colors = []
        color_map = {"TRUE": "#55A868", "FALSE": "#C44E52"}
        label_map = {"TRUE": "Found", "FALSE": "Not Found"}

        for val in ["TRUE", "FALSE"]:
            if val in found_counts.index:
                labels.append(label_map[val])
                sizes.append(found_counts[val])
                colors.append(color_map[val])

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.pie(
            sizes,
            labels=labels,
            colors=colors,
            autopct="%1.1f%%",
            startangle=90,
            textprops={"fontsize": 12},
        )
        ax.set_title("Overall Detection Summary")
        fig.tight_layout()

        output_path = self.output_dir / "overall_summary.png"
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        return output_path
