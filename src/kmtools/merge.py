from __future__ import annotations

from pathlib import Path

import pandas as pd

from kmtools.utils import Utils


class Merge:
    def __init__(
        self,
        inputs: list[str],
        output: str,
        keep: bool = False,
        verbose: bool = False,
    ) -> None:
        self.inputs = inputs
        self.output = Path(output)
        self.keep = keep
        self.verbose = verbose

    def run(self) -> None:
        dataframes = []

        for input_glob in self.inputs:
            for input_file in Path().glob(input_glob):
                Utils.log(f"Reading input file: {input_file}", self.verbose)
                df = pd.read_csv(input_file, sep="\t", comment="#")
                dataframes.append(df)

        merged_df = pd.concat(dataframes, ignore_index=True)
        merged_df.to_csv(self.output, index=False, sep="\t")

        if not self.keep:
            for input_glob in self.inputs:
                for input_file in Path().glob(input_glob):
                    input_file.unlink()
