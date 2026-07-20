"""Computable, independent quality classification."""

from __future__ import annotations

import tomllib
from pathlib import Path

import numpy as np
from rdkit import Chem

from .config import BuildConfig
from .models import ExtractedPocket, LigandGeometry, ParsedProtein

DEFAULT_EFFECTS = {
    "DETERMINISTIC_CROP_APPLIED": "B", "LIGAND_MOL2_FALLBACK": "B",
    "SDF_MOL2_DISAGREEMENT": "B", "MODIFIED_RESIDUE_INCLUDED": "B",
    "EXPLICIT_ELEMENT_MAPPING_APPLIED": "B",
    "PROBABLE_COVALENT_CONTACT": "C", "EXCLUDED_COMPONENT_BRIDGES_CONTACT": "C",
    "VERY_SMALL_POCKET": "C", "UNSUPPORTED_ATOM_EXCLUDED": "C",
    "SPATIALLY_SEPARATED_LIGAND_COMPONENTS": "C", "LOCAL_MISSING_ATOM_RECORD": "C",
    "LIGAND_CHEMISTRY_UNUSUAL": "C",
    "MISSING_INDEX_FILE": "none", "MISSING_COMPLEX_DIRECTORY": "none",
    "DUPLICATE_COMPLEX_DIRECTORY": "none", "SOURCE_FILENAME_CASE_MISMATCH": "none",
    "MISSING_EXPECTED_FILE": "none", "EMPTY_SOURCE_FILE": "none",
    "EXTRA_COMPLEX_FILE": "none", "GEOMETRY_PROCESSING_FAILED": "rejected",
    "LIGAND_BOTH_FORMATS_FAILED": "rejected", "PROTEIN_PARSE_FAILED": "rejected",
    "POCKET_EXTRACTION_FAILED": "rejected",
    "RCSB_ENRICHMENT_FAILED": "none",
    "INVALID_REFERENCE_OVERRIDE": "none", "CONFLICTING_REFERENCE_OVERRIDES": "none",
    "PDBBIND_POCKET_COMPARISON_UNAVAILABLE": "none", "ALTERNATE_LOCATIONS_DISCARDED": "none",
    "PROTEIN_ELEMENTS_INFERRED": "none",
}


def classify_geometry(
    pocket: ExtractedPocket, ligand: LigandGeometry, protein: ParsedProtein, config: BuildConfig
) -> str:
    codes = set(pocket.warning_codes)
    ligand_heavy = ligand.coordinates[ligand.atomic_numbers > 1]
    for atom in pocket.contact_atoms:
        if atom.element == "H":
            continue
        distances = np.linalg.norm(ligand_heavy - atom.coordinate, axis=1)
        ligand_atomic = ligand.atomic_numbers[ligand.atomic_numbers > 1]
        radii = np.asarray([Chem.GetPeriodicTable().GetRcovalent(int(number))
                            + Chem.GetPeriodicTable().GetRcovalent(atom.element)
                            + config.quality.covalent_radius_margin_angstrom for number in ligand_atomic])
        if np.any(distances <= radii):
            codes.add("PROBABLE_COVALENT_CONTACT")
            break
    component_count = int(ligand.component_ids.max()) + 1 if len(ligand.component_ids) else 0
    for first in range(component_count):
        for second in range(first + 1, component_count):
            left = ligand.coordinates[(ligand.component_ids == first) & (ligand.atomic_numbers > 1)]
            right = ligand.coordinates[(ligand.component_ids == second) & (ligand.atomic_numbers > 1)]
            if len(left) and len(right) and np.min(np.linalg.norm(left[:, None] - right[None, :], axis=2)) > config.quality.separated_component_cutoff_angstrom:
                codes.add("SPATIALLY_SEPARATED_LIGAND_COMPONENTS")
    contact_coords = np.asarray([atom.coordinate for atom in pocket.contact_atoms if atom.element != "H"])
    for atom in protein.excluded_atoms:
        if atom.element == "H" or atom.residue_name in {"HOH", "WAT", "DOD"} or not len(contact_coords):
            continue
        ligand_distance = float(np.min(np.linalg.norm(ligand_heavy - atom.coordinate, axis=1)))
        protein_distance = float(np.min(np.linalg.norm(contact_coords - atom.coordinate, axis=1)))
        if ligand_distance <= config.quality.excluded_component_bridge_cutoff_angstrom and protein_distance <= config.quality.excluded_component_bridge_cutoff_angstrom:
            codes.add("EXCLUDED_COMPONENT_BRIDGES_CONTACT")
            break
    effects = load_quality_effects(config.quality.rules_file)
    unknown = codes - set(effects)
    if unknown:
        raise ValueError(f"Unknown geometry issue codes: {sorted(unknown)}")
    pocket.warning_codes = sorted(codes)
    pocket.geometry_quality_tier = "C" if any(effects[code] == "C" for code in codes) else (
        "B" if any(effects[code] == "B" for code in codes) else "A"
    )
    return pocket.geometry_quality_tier


def load_quality_effects(path: Path) -> dict[str, str]:
    if not path.is_file():
        return dict(DEFAULT_EFFECTS)
    raw = tomllib.loads(path.read_text(encoding="utf-8"))
    rules = raw.get("rules", {})
    effects = {key: value["tier"] for key, value in rules.items()}
    return effects
