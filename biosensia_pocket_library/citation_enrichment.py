"""Citation enrichment entry point (mmCIF-first and geometry-independent)."""

from .rcsb import enrich_from_mmcif

__all__ = ["enrich_from_mmcif"]
