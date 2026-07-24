"""One-to-one comparison of re-extracted and PDBbind-provided pocket views."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from .config import BuildConfig
from .models import ExtractedPocket, LigandGeometry, ProteinAtom
from .protein_parser import parse_protein


def compare_pdbbind_pocket(
    pocket: ExtractedPocket, ligand: LigandGeometry, path: Path | None,
    source_file_id: str | None, config: BuildConfig,
) -> tuple[list[dict], list[dict], str]:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return [_unavailable_row(pocket, source_file_id, view) for view in _views()], [], "unavailable"
    try:
        parsed = parse_protein(path, config, pdb_id=pocket.pdb_id, complex_id=pocket.complex_id)
        provided = [atom for atom in parsed.atoms if atom.element != "H"]
    except Exception:
        return [_unavailable_row(pocket, source_file_id, view, "parse_failed") for view in _views()], [], "unavailable"
    provided_min = _distances(provided, ligand)
    provided_contact = [atom for atom in provided if provided_min[atom.pdbbind_atom_key] <= config.pocket.distance_cutoff_angstrom]
    rows: list[dict] = []
    differences: list[dict] = []
    qualities: list[str] = []
    for view_name, extracted in _views(pocket):
        provided_view = provided_contact if view_name == "contact_atom" else provided
        provided_view_min = {atom.pdbbind_atom_key: provided_min[atom.pdbbind_atom_key] for atom in provided_view}
        row, atom_rows = _compare_view(pocket, view_name, extracted, provided_view, provided_view_min, source_file_id, config)
        rows.append(row)
        differences.extend(atom_rows)
        qualities.append(_quality(row, config))
    # The residue-expanded view is primary; contact-atom is a boundary diagnostic.
    overall = qualities[1]
    return rows, differences, overall


def _views(pocket: ExtractedPocket | None = None):
    if pocket is None:
        return ("contact_atom", "residue_expanded_atom")
    return (("contact_atom", pocket.contact_atoms),
            ("residue_expanded_atom", pocket.residue_expanded_atoms))


def _compare_view(pocket, view, extracted, provided, provided_min, source_file_id, config):
    extracted_by_key = {atom.pdbbind_atom_key: atom for atom in extracted}
    provided_by_key = {atom.pdbbind_atom_key: atom for atom in provided}
    exact_keys = set(extracted_by_key) & set(provided_by_key)
    unmatched_left = [atom for key, atom in extracted_by_key.items() if key not in exact_keys]
    unmatched_right = [atom for key, atom in provided_by_key.items() if key not in exact_keys]
    fallback: list[tuple[ProteinAtom, ProteinAtom, float]] = []
    used_right: set[str] = set()
    ambiguous_assignment = False
    left_groups: dict[tuple, list[ProteinAtom]] = defaultdict(list)
    right_groups: dict[tuple, list[ProteinAtom]] = defaultdict(list)
    for atom in unmatched_left:
        left_groups[_fallback_group(atom)].append(atom)
    for atom in unmatched_right:
        right_groups[_fallback_group(atom)].append(atom)
    tolerance = config.structure.coordinate_match_tolerance_angstrom
    for group in sorted(set(left_groups) & set(right_groups)):
        left_values = sorted(left_groups[group], key=lambda atom: atom.pdbbind_atom_key)
        right_values = sorted(right_groups[group], key=lambda atom: atom.pdbbind_atom_key)
        costs = np.asarray([[np.linalg.norm(left.coordinate - right.coordinate)
                             for right in right_values] for left in left_values], dtype=np.float64)
        # A deterministic infinitesimal lexicographic term breaks equal-cost assignments.
        tie_break = np.arange(costs.size, dtype=np.float64).reshape(costs.shape) * 1e-12
        row_indices, column_indices = linear_sum_assignment(costs + tie_break)
        for row_index, column_index in zip(row_indices, column_indices, strict=True):
            distance = float(costs[row_index, column_index])
            if distance > tolerance:
                continue
            left, right = left_values[row_index], right_values[column_index]
            fallback.append((left, right, distance))
            used_right.add(right.pdbbind_atom_key)
            if (np.count_nonzero(np.isclose(costs[row_index], distance, atol=1e-12)) > 1 or
                    np.count_nonzero(np.isclose(costs[:, column_index], distance, atol=1e-12)) > 1):
                ambiguous_assignment = True
    matched_left = exact_keys | {left.pdbbind_atom_key for left, _, _ in fallback}
    matched_right = exact_keys | used_right
    only_left = set(extracted_by_key) - matched_left
    only_right = set(provided_by_key) - matched_right
    exact_distances = [float(np.linalg.norm(extracted_by_key[key].coordinate - provided_by_key[key].coordinate)) for key in exact_keys]
    coordinate_distances = exact_distances + [distance for _, _, distance in fallback]
    left_residues = {atom.residue_key for atom in extracted}
    right_residues = {atom.residue_key for atom in provided}
    union_atoms = len(matched_left) + len(only_left) + len(only_right)
    union_residues = len(left_residues | right_residues)
    re_distances = [pocket.minimum_distances[atom.pdbbind_atom_key] for atom in extracted
                    if atom.pdbbind_atom_key in pocket.minimum_distances]
    pd_distances = list(provided_min.values())
    warnings: list[str] = []
    if coordinate_distances and max(coordinate_distances) > config.comparison.coordinate_rmsd_warning_angstrom:
        warnings.append("COMMON_ATOM_COORDINATE_DIFFERENCE")
    if ambiguous_assignment:
        warnings.append("AMBIGUOUS_FALLBACK_ATOM_ASSIGNMENT")
    row = {
        "pocket_instance_id": pocket.pocket_instance_id, "comparison_view": view,
        "pdbbind_pocket_file_id": source_file_id, "comparison_status": "compared",
        "left_geometry_role": "pdbbind_reextracted", "right_geometry_role": "pdbbind_provided_pocket",
        "reextracted_atom_count": len(extracted), "pdbbind_atom_count": len(provided),
        "reextracted_heavy_atom_count": sum(a.element != "H" for a in extracted),
        "pdbbind_heavy_atom_count": sum(a.element != "H" for a in provided),
        "common_atom_exact_count": len(exact_keys), "common_atom_fallback_count": len(fallback),
        "only_reextracted_atom_count": len(only_left), "only_pdbbind_atom_count": len(only_right),
        "atom_jaccard": len(matched_left) / union_atoms if union_atoms else 1.0,
        "reextracted_residue_count": len(left_residues), "pdbbind_residue_count": len(right_residues),
        "common_residue_count": len(left_residues & right_residues),
        "only_reextracted_residue_count": len(left_residues - right_residues),
        "only_pdbbind_residue_count": len(right_residues - left_residues),
        "residue_jaccard": len(left_residues & right_residues) / union_residues if union_residues else 1.0,
        "reextracted_chain_ids": sorted({a.auth_chain_id for a in extracted}),
        "pdbbind_chain_ids": sorted({a.auth_chain_id for a in provided}),
        "chain_sets_equal": {a.auth_chain_id for a in extracted} == {a.auth_chain_id for a in provided},
        "reextracted_subset_of_pdbbind": not only_left, "pdbbind_subset_of_reextracted": not only_right,
        "common_atom_coordinate_rmsd": float(np.sqrt(np.mean(np.square(coordinate_distances)))) if coordinate_distances else None,
        "common_atom_max_coordinate_difference": max(coordinate_distances, default=None),
        **_distance_summary("reextracted", re_distances), **_distance_summary("pdbbind", pd_distances),
        "warning_codes": warnings,
    }
    atom_rows: list[dict] = []
    for key in sorted(exact_keys):
        atom_rows.append(_difference_row(pocket, view, "matched_exact_identity", extracted_by_key[key], provided_by_key[key], "exact_key",
                                         float(np.linalg.norm(extracted_by_key[key].coordinate - provided_by_key[key].coordinate))))
    for left, right, distance in fallback:
        atom_rows.append(_difference_row(pocket, view, "matched_fallback", left, right, "coordinate_element", distance))
    for key in sorted(only_left):
        atom_rows.append(_difference_row(pocket, view, "only_reextracted", extracted_by_key[key], None, None, None))
    for key in sorted(only_right):
        atom_rows.append(_difference_row(pocket, view, "only_pdbbind", None, provided_by_key[key], None, None))
    return row, atom_rows


def _distance_summary(prefix: str, values: list[float]) -> dict:
    array = np.asarray(values)
    return {f"{prefix}_maximum_ligand_distance": float(np.max(array)) if len(array) else None,
            f"{prefix}_mean_ligand_distance": float(np.mean(array)) if len(array) else None,
            f"{prefix}_median_ligand_distance": float(np.median(array)) if len(array) else None,
            f"{prefix}_p95_ligand_distance": float(np.percentile(array, 95)) if len(array) else None}


def _difference_row(pocket, view, classification, left, right, method, distance):
    atom = left or right
    return {"pocket_instance_id": pocket.pocket_instance_id, "comparison_view": view,
            "comparison_class": classification,
            "reextracted_atom_key": left.pdbbind_atom_key if left else None,
            "pdbbind_atom_key": right.pdbbind_atom_key if right else None,
            "match_method": method, "coordinate_distance": distance,
            "auth_chain_id": atom.auth_chain_id, "auth_residue_number": atom.auth_residue_number,
            "insertion_code": atom.insertion_code, "residue_name": atom.residue_name,
            "atom_name": atom.atom_name, "element": atom.element}


def _distances(atoms: list[ProteinAtom], ligand: LigandGeometry) -> dict[str, float]:
    heavy = ligand.coordinates[ligand.atomic_numbers > 1]
    return {atom.pdbbind_atom_key: float(np.min(np.linalg.norm(heavy - atom.coordinate, axis=1))) for atom in atoms}


def _fallback_group(atom: ProteinAtom) -> tuple:
    return (atom.auth_chain_id, atom.auth_residue_number, atom.insertion_code,
            atom.residue_name, atom.atom_name, atom.element)


def _quality(row: dict, config: BuildConfig) -> str:
    if row["atom_jaccard"] < config.comparison.atom_jaccard_severe_minimum or row["residue_jaccard"] < config.comparison.residue_jaccard_severe_minimum:
        return "severe_difference"
    if row["atom_jaccard"] < config.comparison.atom_jaccard_moderate_minimum or row["residue_jaccard"] < config.comparison.residue_jaccard_moderate_minimum:
        return "moderate_difference"
    if (row["common_atom_max_coordinate_difference"] or 0.0) > config.comparison.coordinate_rmsd_warning_angstrom:
        return "moderate_difference"
    return "concordant"


def _unavailable_row(pocket, source_file_id, view, status="unavailable"):
    names = ["reextracted_atom_count", "pdbbind_atom_count", "reextracted_heavy_atom_count",
             "pdbbind_heavy_atom_count", "common_atom_exact_count", "common_atom_fallback_count",
             "only_reextracted_atom_count", "only_pdbbind_atom_count", "atom_jaccard",
             "reextracted_residue_count", "pdbbind_residue_count", "common_residue_count",
             "only_reextracted_residue_count", "only_pdbbind_residue_count", "residue_jaccard",
             "chain_sets_equal", "reextracted_subset_of_pdbbind", "pdbbind_subset_of_reextracted",
             "common_atom_coordinate_rmsd", "common_atom_max_coordinate_difference",
             "reextracted_maximum_ligand_distance", "pdbbind_maximum_ligand_distance",
             "reextracted_mean_ligand_distance", "pdbbind_mean_ligand_distance",
             "reextracted_median_ligand_distance", "pdbbind_median_ligand_distance",
             "reextracted_p95_ligand_distance", "pdbbind_p95_ligand_distance"]
    row = {name: None for name in names}
    row.update({"pocket_instance_id": pocket.pocket_instance_id, "comparison_view": view,
                "pdbbind_pocket_file_id": source_file_id, "comparison_status": status,
                "left_geometry_role": "pdbbind_reextracted", "right_geometry_role": "pdbbind_provided_pocket",
                "reextracted_chain_ids": [], "pdbbind_chain_ids": [], "warning_codes": ["PDBBIND_POCKET_COMPARISON_UNAVAILABLE"]})
    return row
