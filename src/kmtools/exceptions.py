"""Custom exceptions for kmtools."""


class KmtoolsError(Exception):
    """Base exception for all kmtools errors."""


class FileValidationError(KmtoolsError):
    """Raised when a required file is missing or has invalid format."""


class KmNotFoundError(KmtoolsError):
    """Raised when the km binary is not found or not functional."""


class ChunkValidationError(KmtoolsError):
    """Raised when chunk input validation fails (directory structure, file counts)."""


class MergeError(KmtoolsError):
    """Raised when merge operations fail (no input files, empty results)."""


class FilterError(KmtoolsError):
    """Raised when filtering fails (missing columns, invalid data)."""
