"""Tests for kmtools.merge."""

import pytest
import pandas as pd
from pathlib import Path

from kmtools.merge import Merge
from kmtools.exceptions import MergeError


class TestMergeRun:
    def _write_tsv(self, path, rows):
        """Helper to write a TSV file with standard km output headers."""
        header = "Database\tQuery\tType\n"
        lines = [f"{r[0]}\t{r[1]}\t{r[2]}\n" for r in rows]
        path.write_text(header + "".join(lines))

    def test_merge_two_files(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        f1 = tmp_path / "chunk_0.txt"
        f2 = tmp_path / "chunk_1.txt"
        self._write_tsv(f1, [("file.jf", "chr1_100_200", "Substitution")])
        self._write_tsv(f2, [("file.jf", "chr2_300_400", "Deletion")])

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["chunk_*.txt"],
            output=str(output),
            keep=True,
            verbose=False,
        ).run()

        result = pd.read_csv(output, sep="\t")
        assert len(result) == 2
        assert set(result["Query"]) == {"chr1_100_200", "chr2_300_400"}

    def test_merge_deletes_inputs_when_keep_false(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        f1 = tmp_path / "chunk_0.txt"
        self._write_tsv(f1, [("file.jf", "chr1_100_200", "Substitution")])

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["chunk_*.txt"],
            output=str(output),
            keep=False,
            verbose=False,
        ).run()

        assert output.exists()
        assert not f1.exists()

    def test_merge_keeps_inputs_when_keep_true(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        f1 = tmp_path / "chunk_0.txt"
        self._write_tsv(f1, [("file.jf", "chr1_100_200", "Substitution")])

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["chunk_*.txt"],
            output=str(output),
            keep=True,
            verbose=False,
        ).run()

        assert output.exists()
        assert f1.exists()

    def test_merge_multiple_rows_per_file(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        f1 = tmp_path / "chunk_0.txt"
        self._write_tsv(
            f1,
            [
                ("file.jf", "chr1_100_200", "Substitution"),
                ("file.jf", "chr1_300_400", "Deletion"),
            ],
        )
        f2 = tmp_path / "chunk_1.txt"
        self._write_tsv(f2, [("file.jf", "chr2_500_600", "Insertion")])

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["chunk_*.txt"],
            output=str(output),
            keep=True,
            verbose=False,
        ).run()

        result = pd.read_csv(output, sep="\t")
        assert len(result) == 3

    def test_merge_no_matching_files(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        output = tmp_path / "merged.txt"
        with pytest.raises(MergeError):
            Merge(
                inputs=["nonexistent_*.txt"],
                output=str(output),
                keep=True,
                verbose=False,
            ).run()

    def test_merge_multiple_glob_patterns(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        f1 = tmp_path / "a_output.txt"
        f2 = tmp_path / "b_result.txt"
        self._write_tsv(f1, [("file.jf", "chr1_100_200", "Substitution")])
        self._write_tsv(f2, [("file.jf", "chr2_300_400", "Deletion")])

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["a_*.txt", "b_*.txt"],
            output=str(output),
            keep=True,
            verbose=False,
        ).run()

        result = pd.read_csv(output, sep="\t")
        assert len(result) == 2
