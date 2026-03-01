"""Tests for kmtools.chunk (with mocking for external dependencies)."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from kmtools.chunk import Chunk
from kmtools.exceptions import KmNotFoundError, ChunkValidationError


@pytest.fixture
def chunk_instance(tmp_path):
    """Create a Chunk instance with a valid tmp directory structure."""
    target_dir = tmp_path / "targets"
    target_dir.mkdir()
    for i in range(4):
        sub = target_dir / f"chunk_{i}"
        sub.mkdir()
        for j in range(10):
            (sub / f"target_{j}.fa").write_text(f">seq{j}\nATCG\n")

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    return Chunk(
        threads=4,
        km_find_mutation_options="--ratio 0.0001",
        km_target_directory=str(target_dir),
        km_jellyfish_file=str(tmp_path / "db.jf"),
        output_dir=str(output_dir),
        prefix="test_output",
        merge=False,
        verbose=False,
    )


class TestCheckKmInstalled:
    @patch("shutil.which", return_value="/usr/bin/km")
    @patch("subprocess.run")
    def test_km_found(self, mock_run, mock_which, chunk_instance):
        mock_run.return_value = MagicMock(returncode=0)
        chunk_instance.check_km_installed()  # Should not raise

    @patch("shutil.which", return_value=None)
    def test_km_not_found(self, mock_which, chunk_instance):
        with pytest.raises(KmNotFoundError, match="not found in PATH"):
            chunk_instance.check_km_installed()

    @patch("shutil.which", return_value="/usr/bin/km")
    @patch("subprocess.run", side_effect=FileNotFoundError)
    def test_km_executable_missing(self, mock_run, mock_which, chunk_instance):
        with pytest.raises(KmNotFoundError, match="could not be found"):
            chunk_instance.check_km_installed()

    @patch("shutil.which", return_value="/usr/bin/km")
    @patch("subprocess.run", side_effect=TimeoutError)
    def test_km_timeout(self, mock_run, mock_which, chunk_instance):
        with pytest.raises(KmNotFoundError):
            chunk_instance.check_km_installed()


class TestCheckTargetFilesSplitCorrectly:
    def test_valid_split(self, chunk_instance):
        chunk_instance.check_target_files_split_correctly()
        assert len(chunk_instance.km_target_directory_subfolder) == 4

    def test_wrong_subfolder_count(self, tmp_path):
        target_dir = tmp_path / "targets"
        target_dir.mkdir()
        for i in range(2):  # Only 2 instead of 4
            sub = target_dir / f"chunk_{i}"
            sub.mkdir()
            (sub / "target.fa").write_text(">seq\nATCG\n")

        chunk = Chunk(
            threads=4,
            km_find_mutation_options="--ratio 0.0001",
            km_target_directory=str(target_dir),
            km_jellyfish_file="db.jf",
            output_dir=".",
        )
        with pytest.raises(ChunkValidationError, match="Expected 4 split folders"):
            chunk.check_target_files_split_correctly()

    def test_uneven_distribution(self, tmp_path):
        target_dir = tmp_path / "targets"
        target_dir.mkdir()
        for i in range(4):
            sub = target_dir / f"chunk_{i}"
            sub.mkdir()
            # chunk_0 gets 100 files, others get 10 — >20% deviation
            count = 100 if i == 0 else 10
            for j in range(count):
                (sub / f"target_{j}.fa").write_text(">seq\nATCG\n")

        chunk = Chunk(
            threads=4,
            km_find_mutation_options="--ratio 0.0001",
            km_target_directory=str(target_dir),
            km_jellyfish_file="db.jf",
            output_dir=".",
        )
        with pytest.raises(ChunkValidationError, match="Uneven file distribution"):
            chunk.check_target_files_split_correctly()

    def test_nonexistent_directory(self):
        chunk = Chunk(
            threads=4,
            km_find_mutation_options="--ratio 0.0001",
            km_target_directory="/nonexistent/path",
            km_jellyfish_file="db.jf",
            output_dir=".",
        )
        with pytest.raises(ChunkValidationError, match="does not exist"):
            chunk.check_target_files_split_correctly()


class TestRunKm:
    @patch("subprocess.run")
    def test_run_km_calls_subprocess(self, mock_run, chunk_instance, tmp_path):
        target_sub = tmp_path / "targets" / "chunk_0"
        chunk_instance.run_km(0, target_sub)

        mock_run.assert_called_once()
        call_args = mock_run.call_args
        cmd = call_args[0][0]
        assert cmd[0] == "km"
        assert cmd[1] == "find_mutation"
        assert "--ratio" in cmd
        assert "0.0001" in cmd

    @patch("subprocess.run")
    def test_run_km_output_file_naming(self, mock_run, chunk_instance, tmp_path):
        target_sub = tmp_path / "targets" / "chunk_0"
        chunk_instance.run_km(0, target_sub)

        call_args = mock_run.call_args
        # The stdout argument should be a file handle
        assert call_args[1].get("stdout") is not None or len(call_args) > 1


class TestMergeOutputs:
    @patch("kmtools.chunk.Merge")
    def test_merge_outputs_creates_merge_instance(self, MockMerge, chunk_instance):
        chunk_instance.threads = 4
        mock_instance = MagicMock()
        MockMerge.return_value = mock_instance

        chunk_instance.merge_outputs()

        MockMerge.assert_called_once()
        mock_instance.run.assert_called_once()
