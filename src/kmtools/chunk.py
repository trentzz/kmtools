from concurrent.futures import ThreadPoolExecutor
import math
import os
from pathlib import Path
import shutil
import subprocess
from kmtools.utils import Utils
from kmtools.merge import Merge

class Chunk:
    def __init__(
        self,
        threads,
        km_find_mutation_options,
        km_target_directory,
        km_jellyfish_file,
        output_dir,
        prefix="km_find_mutation_output",
        merge=False,
        merge_output="km_find_mutation_merged_output.txt",
        merge_keep=False,
        verbose=False,
    ):
        self.threads = threads
        self.km_options = km_find_mutation_options
        self.km_target_directory = km_target_directory
        self.km_jellyfish_file = km_jellyfish_file
        self.output_dir = output_dir
        self.prefix = prefix
        self.merge = merge
        self.merge_output = merge_output
        self.merge_keep = merge_keep
        self.verbose = verbose
        
        self.km_target_directory_subfolder = []

    def run(self) -> None:
        Utils.log(f"Starting Chunk with {self.threads} threads", self.verbose)
        Utils.log(f"KM options: {self.km_options}", self.verbose)

        Utils.log("Checking that km is installed...", self.verbose)

        self.check_km_installed()
        
        self.check_target_files_split_correctly()

        # Run km instances in parallel
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            if not hasattr(self, "km_target_directory_subfolder"):
                raise RuntimeError("Split subfolders not found. Ensure check_target_files_split_correctly() was called")
            # Submit one job per split folder, passing (thread_num, target_subfolder)
            futures = [
                executor.submit(self.run_km, i, sub)
                for i, sub in enumerate(self.km_target_directory_subfolder)
            ]
            # Wait for all to complete and propagate exceptions if any
            for fut in futures:
                fut.result()

        if self.merge:
            self.merge_outputs()
            
            # Optionally remove the output directory if empty
            output_dir_path = Path(self.output_dir)
            if output_dir_path.exists() and output_dir_path.is_dir():
                try:
                    if not any(output_dir_path.iterdir()):
                        Utils.log(f"Removing empty output directory {output_dir_path}", self.verbose)
                        output_dir_path.rmdir()
                except Exception as e:
                    Utils.log(f"Could not remove directory {output_dir_path}: {e}", self.verbose)

        # "--count 2 --ratio 0.00001 /mnt/tzeng-local/tzeng-thesis/titration-target-sequences/70bp-targets /mnt/tzeng-local/tzeng-thesis/RaScALL/output/dedup-humid-v3/TWIST_STDV2_30ng_VAF_2pc_22KVL2LT3_CCGCTACCAA-TTACCTCAGT_L008_HUMID_DEDUPED/TWIST_STDV2_30ng_VAF_2pc_22KVL2LT3_CCGCTACCAA-TTACCTCAGT_L008_HUMID_DEDUPED_countTable31.jf"

    def run_km(self, thread_num, target_subfolder):
        output_file = Path(self.output_dir) / f"{self.prefix}_{thread_num}.txt"

        cmd = (
            [
                "km",
                "find_mutation",
            ]
            + self.km_options.split(" ")
            + [
                str(target_subfolder),
                self.km_jellyfish_file,
            ]
        )

        with open(output_file, "w", encoding="utf-8") as f:
            subprocess.run(cmd, stdout=f, check=True)
            
    def check_target_files_split_correctly(self) -> None:
        # Get the directory and check it exists
        target_dir = Path(self.km_target_directory)
        if not target_dir.exists() or not target_dir.is_dir():
            raise RuntimeError(f"Target directory {target_dir} does not exist or is not a directory")

        # Get all subfolders
        subfolders = [d for d in target_dir.iterdir() if d.is_dir()]

        # Check number of subfolders matches thread count
        if len(subfolders) != self.threads:
            raise RuntimeError(f"Expected {self.threads} split folders in {target_dir}, but found {len(subfolders)}")

        # Check file counts in each subfolder
        file_counts = [len(list(d.glob('*'))) for d in subfolders]
        avg_count = sum(file_counts) / len(file_counts)
        max_deviation = avg_count * 0.2  # Allow 20% deviation

        for count in file_counts:
            if abs(count - avg_count) > max_deviation:
                raise RuntimeError(f"Uneven file distribution detected in split folders. Counts: {file_counts}")

        # Store the subfolders
        self.km_target_directory_subfolder = sorted(subfolders)

    def check_km_installed(self) -> None:
        """
        Verify that the `km` command-line tool is available and runnable.
        Returns the resolved path to `km` if successful.
        Raises RuntimeError if km is not found or not functional.
        """
        km_path = shutil.which("km")
        if km_path is None:
            raise RuntimeError(
                "'km' command not found in PATH. Please install it or ensure it's accessible."
            )

        try:
            result = subprocess.run(
                ["km", "--help"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5,
                check=True,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(
                "'km' executable could not be found after locating its path."
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                "'km --help' timed out. Check that km is installed correctly."
            ) from exc
        except Exception as e:
            raise RuntimeError(f"Unexpected error when checking 'km': {e}") from e

        if result.returncode != 0:
            raise RuntimeError(
                f"'km' command returned nonzero exit code ({result.returncode}). "
                "Please verify your installation."
            )

        Utils.log(f"'km' found and verified at {km_path}", self.verbose)

    def merge_outputs(self) -> None:
        Utils.log("Merging output files...", self.verbose)
        
        Merge(
            inputs=[str(Path(self.output_dir) / f"{self.prefix}_{i}.txt") for i in range(self.threads)],
            output=self.merge_output,
            keep=self.merge_keep,
            verbose=self.verbose,
        ).run()
        