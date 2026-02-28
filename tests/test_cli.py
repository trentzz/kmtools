"""Tests for CLI argument parsing in kmtools.kmtools."""

import pytest
from unittest.mock import patch
import sys


class TestChunkParser:
    def test_chunk_required_args(self):
        """Test that chunk subcommand parses required arguments."""
        from kmtools.kmtools import main

        test_args = [
            "kmtools",
            "chunk",
            "--threads",
            "4",
            "--km-find-mutation-options",
            "--ratio 0.0001",
            "--km-target-directory",
            "/tmp/targets",
            "--km-jellyfish-file",
            "/tmp/db.jf",
        ]
        with patch.object(sys, "argv", test_args):
            with patch("kmtools.kmtools.run_chunk") as mock_run:
                main()
                args = mock_run.call_args[0][0]
                assert args.threads == 4
                assert args.km_find_mutation_options == "--ratio 0.0001"
                assert args.km_target_directory == "/tmp/targets"
                assert args.km_jellyfish_file == "/tmp/db.jf"

    def test_chunk_defaults(self):
        from kmtools.kmtools import main

        test_args = [
            "kmtools",
            "chunk",
            "--threads",
            "2",
            "--km-find-mutation-options",
            "opts",
            "--km-target-directory",
            "/tmp",
            "--km-jellyfish-file",
            "/tmp/db.jf",
        ]
        with patch.object(sys, "argv", test_args):
            with patch("kmtools.kmtools.run_chunk") as mock_run:
                main()
                args = mock_run.call_args[0][0]
                assert args.output_dir == "."
                assert args.prefix == "km_find_mutation_output"
                assert args.merge is False
                assert args.merge_keep is False

    def test_chunk_with_merge_flags(self):
        from kmtools.kmtools import main

        test_args = [
            "kmtools",
            "chunk",
            "--threads",
            "2",
            "--km-find-mutation-options",
            "opts",
            "--km-target-directory",
            "/tmp",
            "--km-jellyfish-file",
            "/tmp/db.jf",
            "--merge",
            "--merge-keep",
            "--merge-output",
            "custom_merged.txt",
        ]
        with patch.object(sys, "argv", test_args):
            with patch("kmtools.kmtools.run_chunk") as mock_run:
                main()
                args = mock_run.call_args[0][0]
                assert args.merge is True
                assert args.merge_keep is True
                assert args.merge_output == "custom_merged.txt"


class TestMergeParser:
    def test_merge_required_args(self):
        from kmtools.kmtools import main

        test_args = [
            "kmtools",
            "merge",
            "file1.txt",
            "file2.txt",
            "--output",
            "merged.txt",
        ]
        with patch.object(sys, "argv", test_args):
            with patch("kmtools.kmtools.run_merge") as mock_run:
                main()
                args = mock_run.call_args[0][0]
                assert args.inputs == ["file1.txt", "file2.txt"]
                assert args.output == "merged.txt"
                assert args.keep is False


class TestFilterParser:
    def test_filter_required_args(self):
        from kmtools.kmtools import main

        test_args = [
            "kmtools",
            "filter",
            "--reference",
            "ref.tsv",
            "--km-output",
            "km.txt",
            "--output",
            "filtered.tsv",
        ]
        with patch.object(sys, "argv", test_args):
            with patch("kmtools.kmtools.run_filter") as mock_run:
                main()
                args = mock_run.call_args[0][0]
                assert args.reference == "ref.tsv"
                assert args.km_output == "km.txt"
                assert args.output == "filtered.tsv"
                assert args.output_type == "tsv"
                assert args.count_threshold == 2


class TestPlotParser:
    def test_plot_required_args(self):
        from kmtools.kmtools import main

        test_args = ["kmtools", "plot", "filtered.tsv"]
        with patch.object(sys, "argv", test_args):
            with patch("kmtools.kmtools.run_plot") as mock_run:
                main()
                args = mock_run.call_args[0][0]
                assert args.file == "filtered.tsv"
                assert args.output_dir == "."
                assert args.charts == "all"


class TestNoSubcommand:
    def test_missing_subcommand_exits(self):
        from kmtools.kmtools import main

        test_args = ["kmtools"]
        with patch.object(sys, "argv", test_args):
            with pytest.raises(SystemExit):
                main()
