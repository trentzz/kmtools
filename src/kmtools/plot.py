from __future__ import annotations

from pathlib import Path

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
        self.charts = charts.split(",") if charts != "all" else ["vaf", "patient", "sample", "overall"]
        self.verbose = verbose

    def run(self) -> None:
        Utils.log(f"Generating plots from {self.file} -> {self.output_dir}", self.verbose)
        for chart in self.charts:
            print(f"Generating {chart} chart...")
        print("Plotting complete.")
