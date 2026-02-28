import glob as globmod
from pathlib import Path

import pandas as pd

from kmtools.exceptions import MergeError
from kmtools.utils import Utils


class Merge:
    def __init__(self, inputs, output, keep=False, verbose=False):
        self.inputs = inputs
        self.output = Path(output)
        self.keep = keep
        self.verbose = verbose

    def _resolve_input_files(self) -> list[Path]:
        """Resolve glob patterns to actual file paths, supporting absolute paths."""
        files = []
        for pattern in self.inputs:
            matched = sorted(Path(p) for p in globmod.glob(pattern))
            files.extend(matched)
        return files

    def run(self) -> None:
        input_files = self._resolve_input_files()

        if not input_files:
            raise MergeError(
                f"No input files found matching patterns: {self.inputs}"
            )

        dataframes = []
        for input_file in input_files:
            if not input_file.exists():
                raise MergeError(f"Input file does not exist: {input_file}")
            Utils.log(f"Reading input file: {input_file}", self.verbose)
            df = pd.read_csv(input_file, sep="\t", comment="#")
            dataframes.append(df)

        merged_df = pd.concat(dataframes, ignore_index=True)

        if merged_df.empty:
            Utils.log("Warning: merged result is empty", self.verbose)

        merged_df.to_csv(self.output, index=False, sep="\t")
        Utils.log(f"Merged {len(input_files)} files ({len(merged_df)} rows) -> {self.output}", self.verbose)

        if not self.keep:
            for input_file in input_files:
                input_file.unlink()
                Utils.log(f"Removed {input_file}", self.verbose)
