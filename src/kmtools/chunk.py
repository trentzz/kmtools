from concurrent.futures import ThreadPoolExecutor
import shlex
from pathlib import Path
import shutil
import subprocess

from kmtools.exceptions import ChunkValidationError, KmNotFoundError
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

        self.check_km_installed()
        self._ensure_output_dir()
        self.check_target_files_split_correctly()

        # Validate jellyfish file exists
        jf_path = Path(self.km_jellyfish_file)
        if not jf_path.exists():
            raise ChunkValidationError(
                f"Jellyfish file not found: {self.km_jellyfish_file}"
            )

        # Run km instances in parallel
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = [
                executor.submit(self.run_km, i, sub)
                for i, sub in enumerate(self.km_target_directory_subfolder)
            ]
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

    def _ensure_output_dir(self) -> None:
        """Create output directory if it doesn't exist."""
        output_path = Path(self.output_dir)
        if not output_path.exists():
            Utils.log(f"Creating output directory: {output_path}", self.verbose)
            output_path.mkdir(parents=True, exist_ok=True)

    def run_km(self, thread_num, target_subfolder):
        output_file = Path(self.output_dir) / f"{self.prefix}_{thread_num}.txt"

        cmd = (
            ["km", "find_mutation"]
            + shlex.split(self.km_options)
            + [str(target_subfolder), self.km_jellyfish_file]
        )

        Utils.log(f"Thread {thread_num}: running {' '.join(cmd)}", self.verbose)

        with open(output_file, "w", encoding="utf-8") as f:
            subprocess.run(cmd, stdout=f, check=True)

    def check_target_files_split_correctly(self) -> None:
        target_dir = Path(self.km_target_directory)
        if not target_dir.exists() or not target_dir.is_dir():
            raise ChunkValidationError(
                f"Target directory {target_dir} does not exist or is not a directory"
            )

        subfolders = [d for d in target_dir.iterdir() if d.is_dir()]

        if len(subfolders) != self.threads:
            raise ChunkValidationError(
                f"Expected {self.threads} split folders in {target_dir}, "
                f"but found {len(subfolders)}"
            )

        file_counts = [len(list(d.glob("*"))) for d in subfolders]

        if any(c == 0 for c in file_counts):
            raise ChunkValidationError(
                f"Empty subfolder(s) detected. File counts: {file_counts}"
            )

        avg_count = sum(file_counts) / len(file_counts)
        max_deviation = avg_count * 0.2

        for count in file_counts:
            if abs(count - avg_count) > max_deviation:
                raise ChunkValidationError(
                    f"Uneven file distribution detected in split folders. "
                    f"Counts: {file_counts}"
                )

        self.km_target_directory_subfolder = sorted(subfolders)

    def check_km_installed(self) -> None:
        km_path = shutil.which("km")
        if km_path is None:
            raise KmNotFoundError(
                "'km' command not found in PATH. "
                "Install it with: pipx install km-walk"
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
            raise KmNotFoundError(
                "'km' executable could not be found after locating its path."
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise KmNotFoundError(
                "'km --help' timed out. Check that km is installed correctly."
            ) from exc
        except Exception as e:
            raise KmNotFoundError(
                f"Unexpected error when checking 'km': {e}"
            ) from e

        if result.returncode != 0:
            raise KmNotFoundError(
                f"'km' command returned nonzero exit code ({result.returncode}). "
                "Please verify your installation."
            )

        Utils.log(f"'km' found and verified at {km_path}", self.verbose)

    def merge_outputs(self) -> None:
        Utils.log("Merging output files...", self.verbose)

        Merge(
            inputs=[
                str(Path(self.output_dir) / f"{self.prefix}_{i}.txt")
                for i in range(self.threads)
            ],
            output=self.merge_output,
            keep=self.merge_keep,
            verbose=self.verbose,
        ).run()
