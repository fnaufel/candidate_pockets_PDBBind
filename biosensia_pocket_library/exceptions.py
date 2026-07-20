"""Domain exceptions used by the pocket-library pipeline."""


class PocketLibraryError(Exception):
    """Base class for expected pipeline failures."""


class ConfigurationError(PocketLibraryError):
    """Raised for invalid or incompatible configuration."""


class SourceIntegrityError(PocketLibraryError):
    """Raised when source identity or integrity cannot be established."""


class ParseError(PocketLibraryError):
    """Raised when a required molecular source cannot be parsed."""


class ValidationError(PocketLibraryError):
    """Raised when produced artifacts violate their contract."""
