from __future__ import annotations

import glob as globmod
from pathlib import Path
from typing import Optional

import pandas as pd

from kmtools.exceptions import MergeError
from kmtools.utils import Utils


class Merge:
    def __init__(
        self,
        inputs: list[str],
        output: str,
        keep: bool = False,
        sort_by: Optional[str] = None,
        drop_duplicates: bool = False,
        verbose: bool = False,
    ) -> None:
        self.inputs = inputs
        self.output = Path(output)
        self.keep = keep
        self.sort_by = sort_by
        self.drop_duplicates = drop_duplicates
        self.verbose = verbose

    def _resolve_input_files(self) -> list[Path]:
        """Resolve glob patterns to actual file paths, supporting absolute paths."""
        files = []
        for pattern in self.inputs:
            matched = sorted(Path(p) for p in globmod.glob(pattern))
            files.extend(matched)
        return files

    def _validate_columns(self, dataframes: list[pd.DataFrame], filenames: list[Path]) -> None:
        """Warn if input files have mismatched columns."""
        if len(dataframes) < 2:
            return

        reference_cols = set(dataframes[0].columns)
        for i, df in enumerate(dataframes[1:], start=1):
            current_cols = set(df.columns)
            if current_cols != reference_cols:
                missing = reference_cols - current_cols
                extra = current_cols - reference_cols
                parts = []
                if missing:
                    parts.append(f"missing: {sorted(missing)}")
                if extra:
                    parts.append(f"extra: {sorted(extra)}")
                Utils.log(
                    f"Warning: {filenames[i]} has different columns than {filenames[0]} ({', '.join(parts)})",
                    self.verbose,
                )

    def run(self) -> None:
        input_files = self._resolve_input_files()

        if not input_files:
            raise MergeError(
                f"No input files found matching patterns: {self.inputs}"
            )

        Utils.log(f"Found {len(input_files)} files to merge", self.verbose)

        dataframes = []
        for input_file in input_files:
            if not input_file.exists():
                raise MergeError(f"Input file does not exist: {input_file}")
            Utils.log(f"Reading: {input_file}", self.verbose)
            df = pd.read_csv(input_file, sep="\t", comment="#")
            dataframes.append(df)

        self._validate_columns(dataframes, input_files)

        merged_df = pd.concat(dataframes, ignore_index=True)

        if merged_df.empty:
            Utils.log("Warning: merged result is empty", self.verbose)

        if self.drop_duplicates:
            before = len(merged_df)
            merged_df = merged_df.drop_duplicates()
            dropped = before - len(merged_df)
            if dropped > 0:
                Utils.log(f"Dropped {dropped} duplicate rows", self.verbose)

        if self.sort_by and self.sort_by in merged_df.columns:
            merged_df = merged_df.sort_values(by=self.sort_by).reset_index(drop=True)
            Utils.log(f"Sorted by column: {self.sort_by}", self.verbose)

        merged_df.to_csv(self.output, index=False, sep="\t")
        Utils.log(
            f"Merged {len(input_files)} files ({len(merged_df)} rows) -> {self.output}",
            self.verbose,
        )

        if not self.keep:
            for input_file in input_files:
                input_file.unlink()
                Utils.log(f"Removed: {input_file}", self.verbose)
