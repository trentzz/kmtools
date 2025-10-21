from concurrent.futures import ThreadPoolExecutor
import math
import os
from pathlib import Path
import shutil
import subprocess
from kmtools.utils import Utils


class Chunk:
    def __init__(
        self,
        threads,
        km_find_mutation_options,
        km_target_directory,
        km_jellyfish_file,
        output_dir,
        merge=False,
        verbose=False,
    ):
        self.threads = threads
        self.km_options = km_find_mutation_options
        self.km_target_directory = km_target_directory
        self.km_jellyfish_file = km_jellyfish_file
        self.output_dir = output_dir
        self.merge = merge
        self.verbose = verbose
        
        self.splits_dir = Path(self.output_dir) / "splits"

    def run(self) -> None:
        Utils.log(f"Starting Chunk with {self.threads} threads", self.verbose)
        Utils.log(f"KM options: {self.km_options}", self.verbose)

        Utils.log("Checking that km is installed...", self.verbose)

        self.check_km_installed()

        # Create directories for splits
        os.makedirs(self.splits_dir, exist_ok=True)

        # Get list of files and calculate files per thread
        files = list(Path(self.km_target_directory).glob("*"))
        files_per_thread = math.ceil(len(files) / self.threads)

        # Split files into directories
        for i in range(self.threads):
            thread_dir = self.splits_dir / f"split_{i}"
            os.makedirs(thread_dir, exist_ok=True)
            start_idx = i * files_per_thread
            end_idx = min((i + 1) * files_per_thread, len(files))

            for file in files[start_idx:end_idx]:
                shutil.copy2(file, thread_dir)

        # Run km instances in parallel
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            executor.map(self.run_km, range(self.threads))

        # # Cleanup split directories
        # if not self.merge:
        #     shutil.rmtree(self.splits_dir)

        # "--count 2 --ratio 0.00001 /mnt/tzeng-local/tzeng-thesis/titration-target-sequences/70bp-targets /mnt/tzeng-local/tzeng-thesis/RaScALL/output/dedup-humid-v3/TWIST_STDV2_30ng_VAF_2pc_22KVL2LT3_CCGCTACCAA-TTACCTCAGT_L008_HUMID_DEDUPED/TWIST_STDV2_30ng_VAF_2pc_22KVL2LT3_CCGCTACCAA-TTACCTCAGT_L008_HUMID_DEDUPED_countTable31.jf"

    def run_km(self, thread_num):
        thread_dir = self.splits_dir / f"split_{thread_num}"
        output_file = Path(self.output_dir) / f"output_{thread_num}.txt"

        cmd = (
            [
                "km",
                "find_mutation",
            ]
            + self.km_options.split(" ")
            + [
                str(thread_dir),
                self.km_jellyfish_file,
            ]
        )

        with open(output_file, "w", encoding="utf-8") as f:
            subprocess.run(cmd, stdout=f, check=True)

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
