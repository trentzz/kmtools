from __future__ import annotations

import sys
import time
from typing import Any, Callable


class Utils:
    """Utility functions for logging, timing, and other shared operations."""

    @staticmethod
    def log(message: str, verbose: bool = False) -> None:
        """Print a message if verbose mode is enabled."""
        if verbose:
            print(message, file=sys.stderr, flush=True)

    @staticmethod
    def time_it(label: str, func: Callable, *args: Any, verbose: bool = False, **kwargs: Any) -> Any:
        """Time execution of a function and optionally log the duration."""
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        Utils.log(f"{label} completed in {elapsed:.2f}s", verbose)
        return result
