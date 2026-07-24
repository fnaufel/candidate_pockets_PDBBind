"""Validation for explicitly trusted DrugCLIP pickle source records."""

from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .exceptions import ConfigurationError, ParseError


@dataclass(frozen=True, slots=True)
class CombineSetPickle:
    pocket_name: str
    pocket_atoms: tuple[str, ...]
    pocket_coordinates: np.ndarray
    ligand_atoms: tuple[str, ...]
    ligand_coordinate_sets: tuple[np.ndarray, ...]
    label: Any


def load_trusted_combine_set_pickle(
    path: Path, *, trusted: bool, dictionary: set[str]
) -> CombineSetPickle:
    """Load one vendored pickle after an explicit trust decision."""
    if not trusted:
        raise ConfigurationError(
            "Refusing to deserialize data.pkl without combine_set.trusted_pickles=true"
        )
    try:
        with path.open("rb") as handle:
            record = pickle.load(handle)
    except Exception as error:
        raise ParseError(f"Could not load trusted pickle: {type(error).__name__}: {error}") from error
    if not isinstance(record, dict):
        raise ParseError(f"Expected pickle dictionary, found {type(record).__name__}")
    required = {"pocket", "pocket_atoms", "pocket_coordinates"}
    missing = required - set(record)
    if missing:
        raise ParseError(f"Pickle lacks required fields: {sorted(missing)}")

    pocket_name = str(record["pocket"]).strip()
    if not pocket_name:
        raise ParseError("Pickle pocket field is empty")
    pocket_atoms = _atom_tokens(record["pocket_atoms"], field="pocket_atoms")
    pocket_coordinates = _single_coordinate_array(
        record["pocket_coordinates"], field="pocket_coordinates"
    )
    if len(pocket_atoms) != len(pocket_coordinates):
        raise ParseError("pocket_atoms and pocket_coordinates have different lengths")
    transformed = tuple(drugclip_pocket_token(atom) for atom in pocket_atoms)
    unsupported = set(transformed) - dictionary
    if unsupported:
        raise ParseError(f"DrugCLIP dictionary lacks transformed tokens: {sorted(unsupported)}")

    ligand_atoms = _atom_tokens(record.get("atoms", []), field="atoms", allow_empty=True)
    ligand_coordinate_sets = _coordinate_arrays(record.get("coordinates", []), field="coordinates")
    if ligand_atoms and any(len(item) != len(ligand_atoms) for item in ligand_coordinate_sets):
        raise ParseError("atoms and coordinates have different lengths")
    return CombineSetPickle(
        pocket_name=pocket_name,
        pocket_atoms=pocket_atoms,
        pocket_coordinates=pocket_coordinates,
        ligand_atoms=ligand_atoms,
        ligand_coordinate_sets=ligand_coordinate_sets,
        label=record.get("label"),
    )


def drugclip_pocket_token(atom: str) -> str:
    """Apply ``AffinityPocketDataset.pocket_atom`` without importing Uni-Core."""
    if not atom:
        raise ParseError("Pocket atom token is empty")
    return atom[1] if atom[0].isdigit() and len(atom) > 1 else atom[0]


def _atom_tokens(value: Any, *, field: str, allow_empty: bool = False) -> tuple[str, ...]:
    try:
        atoms = tuple(str(item).strip() for item in value)
    except TypeError as error:
        raise ParseError(f"{field} is not an atom sequence") from error
    if (not atoms and not allow_empty) or any(not atom for atom in atoms):
        raise ParseError(f"{field} is empty or contains an empty token")
    return atoms


def _single_coordinate_array(value: Any, *, field: str) -> np.ndarray:
    arrays = _coordinate_arrays(value, field=field)
    if len(arrays) != 1:
        raise ParseError(f"{field} must contain exactly one coordinate set")
    return arrays[0]


def _coordinate_arrays(value: Any, *, field: str) -> tuple[np.ndarray, ...]:
    array = np.asarray(value)
    try:
        if array.dtype != object and array.ndim == 2 and array.shape[1] == 3:
            arrays = (np.asarray(array, dtype=np.float64, order="C"),)
        elif array.dtype != object and array.ndim == 3 and array.shape[2] == 3:
            arrays = tuple(np.asarray(item, dtype=np.float64, order="C") for item in array)
        else:
            arrays = tuple(np.asarray(item, dtype=np.float64, order="C") for item in value)
    except (TypeError, ValueError) as error:
        raise ParseError(f"{field} is not numeric coordinate data") from error
    if not arrays:
        return ()
    if any(item.ndim != 2 or item.shape[1] != 3 for item in arrays):
        raise ParseError(f"{field} is not one or more (n, 3) arrays")
    if not all(np.isfinite(item).all() for item in arrays):
        raise ParseError(f"{field} contains non-finite coordinates")
    return arrays
