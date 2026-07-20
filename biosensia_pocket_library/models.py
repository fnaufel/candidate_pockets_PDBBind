"""In-memory domain models."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(slots=True)
class ProcessingIssue:
    stage: str
    issue_code: str
    severity: str
    message: str
    complex_id: str | None = None
    pocket_instance_id: str | None = None
    source_file_id: str | None = None
    exception_type: str | None = None
    details: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class BindingMeasurement:
    measurement_id: str
    complex_id: str
    pdb_id: str
    measurement_type_raw: str
    measurement_type_normalized: str | None
    relation_raw: str | None
    relation_normalized: str | None
    value_raw: str | None
    value_numeric: float | None
    unit_raw: str | None
    unit_normalized: str | None
    value_molar: float | None
    value_inverse_molar: float | None
    p_measurement_name: str | None
    p_relation: str | None
    p_value: float | None
    normalization_kind: str | None
    measurement_raw: str
    parse_status: str
    parse_warning_codes: list[str]
    source_index_line_number: int


@dataclass(slots=True)
class IndexRecord:
    complex_id: str
    distribution_id: str
    pdb_id: str
    resolution_raw: str
    resolution_angstrom: float | None
    experimental_method_hint: str | None
    release_year: int
    ligand_label: str | None
    index_comment: str | None
    index_line_redacted: str
    source_line_sha256: str
    primary_index_line_number: int
    occurrence_line_numbers: list[int]
    measurement: BindingMeasurement


@dataclass(slots=True)
class LigandGeometry:
    ligand_instance_id: str
    source_format: str
    source_path: Path
    source_sha256: str
    atomic_numbers: np.ndarray
    elements: list[str]
    coordinates: np.ndarray
    component_ids: np.ndarray
    content_hash: str
    derivation_hash: str
    chemistry: dict[str, Any]
    comparison: dict[str, Any]
    warnings: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ProteinAtom:
    pdbbind_atom_key: str
    source_order: int
    model_id: int
    record_type: str
    auth_chain_id: str
    auth_residue_number: int
    insertion_code: str
    residue_name: str
    atom_name: str
    altloc: str
    element: str
    occupancy: float | None
    b_factor: float | None
    x: float
    y: float
    z: float
    serial_number: int | None
    locally_classified_protein: bool
    exclusion_reason: str | None = None

    @property
    def coordinate(self) -> np.ndarray:
        return np.asarray((self.x, self.y, self.z), dtype=np.float64)

    @property
    def residue_key(self) -> tuple[str, int, str, str]:
        return (
            self.auth_chain_id,
            self.auth_residue_number,
            self.insertion_code,
            self.residue_name,
        )


@dataclass(slots=True)
class ParsedProtein:
    atoms: list[ProteinAtom]
    excluded_atoms: list[ProteinAtom]
    model_count: int
    selected_model_id: int
    discarded_altloc_count: int
    inferred_element_count: int
    warnings: list[str] = field(default_factory=list)
    missing_atom_residues: set[tuple[str, int, str, str]] = field(default_factory=set)


@dataclass(slots=True)
class ExtractedPocket:
    pocket_instance_id: str
    complex_id: str
    ligand_instance_id: str
    pdb_id: str
    content_hash: str
    derivation_hash: str
    contact_atoms: list[ProteinAtom]
    residue_expanded_atoms: list[ProteinAtom]
    exported_atoms: list[ProteinAtom]
    minimum_distances: dict[str, float]
    contact_residues: list[tuple[str, int, str, str]]
    crop_applied: bool
    max_retained_distance: float | None
    min_discarded_distance: float | None
    geometry_quality_tier: str
    warning_codes: list[str]
    error_codes: list[str]
