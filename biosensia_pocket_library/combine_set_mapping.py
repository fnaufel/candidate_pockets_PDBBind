"""Map authoritative DrugCLIP pickle atoms to neighboring PDB metadata."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .config import BuildConfig
from .models import ProteinAtom
from .protein_parser import parse_protein
from .trusted_pickle import drugclip_pocket_token


@dataclass(frozen=True, slots=True)
class AtomMapping:
    source_role: str | None
    source_file_id: str | None
    atoms: tuple[ProteinAtom | None, ...]
    statuses: tuple[str, ...]
    matched_count: int
    maximum_distance: float | None


def map_pickle_atoms(
    raw_atoms: tuple[str, ...], coordinates: np.ndarray, candidates: list[tuple[str, Path, str]],
    config: BuildConfig, *, pdb_id: str, complex_id: str,
) -> AtomMapping:
    mappings = [
        _map_candidate(raw_atoms, coordinates, role, path, file_id, config,
                       pdb_id=pdb_id, complex_id=complex_id)
        for role, path, file_id in candidates if path.is_file() and path.stat().st_size
    ]
    if not mappings:
        return AtomMapping(None, None, (None,) * len(raw_atoms), ("unresolved",) * len(raw_atoms), 0, None)
    return sorted(
        mappings,
        key=lambda item: (-item.matched_count,
                          item.maximum_distance if item.maximum_distance is not None else float("inf"),
                          item.source_role or ""),
    )[0]


def _map_candidate(raw_atoms, coordinates, role, path, file_id, config, *, pdb_id, complex_id):
    try:
        parsed = parse_protein(
            path, config, pdb_id=pdb_id, complex_id=complex_id,
            distribution_id=config.combine_set.distribution_id,
        )
        pdb_atoms = parsed.atoms
    except Exception:
        return AtomMapping(role, file_id, (None,) * len(raw_atoms), ("unresolved",) * len(raw_atoms), 0, None)
    normalized = tuple(drugclip_pocket_token(atom) for atom in raw_atoms)
    if len(pdb_atoms) == len(normalized):
        pdb_tokens = tuple(atom.element for atom in pdb_atoms)
        pdb_xyz = np.asarray([atom.coordinate for atom in pdb_atoms], dtype=np.float64)
        errors = np.linalg.norm(coordinates - pdb_xyz, axis=1)
        if pdb_tokens == normalized and np.all(errors <= config.combine_set.compare_coordinate_tolerance_angstrom):
            return AtomMapping(role, file_id, tuple(pdb_atoms), ("exact_ordered",) * len(raw_atoms),
                               len(raw_atoms), float(errors.max()) if len(errors) else None)

    mapped: list[ProteinAtom | None] = [None] * len(raw_atoms)
    statuses = ["unresolved"] * len(raw_atoms)
    distances: list[float] = []
    tolerance = config.structure.coordinate_match_tolerance_angstrom
    for element in sorted(set(normalized)):
        left = np.asarray([index for index, token in enumerate(normalized) if token == element], dtype=np.int64)
        right = [atom for atom in pdb_atoms if atom.element == element]
        if not len(left) or not right:
            continue
        tree = cKDTree(np.asarray([atom.coordinate for atom in right], dtype=np.float64))
        nearest_distance, nearest_index = tree.query(coordinates[left], k=min(2, len(right)))
        nearest_distance = np.asarray(nearest_distance).reshape((len(left), -1))
        nearest_index = np.asarray(nearest_index).reshape((len(left), -1))
        proposals: dict[int, list[tuple[int, float, bool]]] = {}
        for row, source_index in enumerate(left):
            distance = float(nearest_distance[row, 0])
            if distance > tolerance:
                continue
            target_index = int(nearest_index[row, 0])
            ambiguous = nearest_distance.shape[1] > 1 and np.isclose(
                nearest_distance[row, 0], nearest_distance[row, 1], atol=1e-12, rtol=0
            )
            proposals.setdefault(target_index, []).append((int(source_index), distance, bool(ambiguous)))
        for target_index, values in proposals.items():
            if len(values) != 1:
                for source_index, _, _ in values:
                    statuses[source_index] = "ambiguous"
                continue
            source_index, distance, ambiguous = values[0]
            if ambiguous:
                statuses[source_index] = "ambiguous"
                continue
            mapped[source_index] = right[target_index]
            statuses[source_index] = "spatial_unique"
            distances.append(distance)
    return AtomMapping(role, file_id, tuple(mapped), tuple(statuses),
                       sum(atom is not None for atom in mapped), max(distances, default=None))
