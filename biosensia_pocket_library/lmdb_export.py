"""Deterministic DrugCLIP-compatible single-file LMDB export."""

from __future__ import annotations

import hashlib
import json
import os
import pickle
from pathlib import Path

import lmdb
import numpy as np

from .config import BuildConfig
from .drugclip_contract import read_dictionary
from .hashing import atomic_write_bytes, canonical_json_bytes, length_frame, sha256_bytes, sha256_file
from .progress import track
from .sidecars import read_sidecar


PROFILE_FILTERS = {
    "default": {"tiers": ("A", "B"), "filename": "candidate_pockets.lmdb"},
    "tier-a": {"tiers": ("A",), "filename": "candidate_pockets_tier_a.lmdb"},
    "tiers-ab": {"tiers": ("A", "B"), "filename": "candidate_pockets_tiers_ab.lmdb"},
    "all-usable": {"tiers": ("A", "B", "C"), "filename": "candidate_pockets_all_usable.lmdb"},
}


def export_lmdb(
    run_dir: Path, config: BuildConfig, contract: dict, profile: str = "default",
    *, overwrite: bool = False, progress: bool = True,
) -> tuple[dict, list[dict]]:
    if profile not in PROFILE_FILTERS:
        raise ValueError(f"Unknown LMDB profile {profile!r}")
    definition = dict(PROFILE_FILTERS[profile])
    if profile == "default":
        definition["tiers"] = tuple(config.lmdb.include_geometry_quality_tiers)
    sidecars = run_dir / "sidecars"
    pocket_rows = {row["pocket_instance_id"]: row for row in read_sidecar(sidecars, "pockets")
                   if row["processing_status"] in {"accepted", "accepted_with_warnings"}
                   and row["geometry_quality_tier"] in definition["tiers"]}
    atoms_by_pocket: dict[str, list[dict]] = {key: [] for key in pocket_rows}
    for atom in read_sidecar(sidecars, "pocket_atoms"):
        if atom["pocket_instance_id"] in atoms_by_pocket and atom["retained_after_crop"]:
            atoms_by_pocket[atom["pocket_instance_id"]].append(atom)
    payloads: list[tuple[str, bytes]] = []
    dictionary = read_dictionary(config.paths.drugclip_dictionary)
    for pocket_id in sorted(pocket_rows):
        atom_rows = sorted(atoms_by_pocket[pocket_id], key=lambda row: row["export_order"])
        tokens = [row["element"] for row in atom_rows]
        coords = np.asarray([[row["x"], row["y"], row["z"]] for row in atom_rows], dtype="<f4", order="C")
        _validate_record_fields(pocket_id, tokens, coords, dictionary, config.pocket.max_pocket_atoms)
        value = pickle.dumps({"pocket": pocket_id, "pocket_atoms": tokens,
                              "pocket_coordinates": coords}, protocol=4)
        payloads.append((pocket_id, value))
    map_size = _map_size(payloads, config)
    lmdb_dir = run_dir / "lmdb"
    lmdb_dir.mkdir(parents=True, exist_ok=True)
    destination = lmdb_dir / str(definition["filename"])
    if destination.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite {destination}")
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    temporary.unlink(missing_ok=True)
    rows: list[dict] = []
    logical = hashlib.sha256()
    for index, (pocket_id, value) in enumerate(payloads):
        key = str(index).encode("ascii")
        framed = length_frame((key, value))
        logical.update(framed)
        pocket = pocket_rows[pocket_id]
        rows.append({"library_profile": profile, "lmdb_path": destination.relative_to(run_dir).as_posix(),
                     "record_index": index, "lmdb_key": key.decode("ascii"),
                     "pocket_instance_id": pocket_id,
                     "pocket_geometry_content_hash": pocket["pocket_geometry_content_hash"],
                     "pocket_derivation_hash": pocket["pocket_derivation_hash"],
                     "atom_count": len(pickle.loads(value)["pocket_atoms"]),
                     "serialized_record_sha256": sha256_bytes(value),
                     "logical_record_sha256": sha256_bytes(framed)})
    try:
        map_size = _write_lmdb(temporary, payloads, map_size, auto=config.lmdb.map_size == "auto",
                               profile=profile, progress=progress)
        _validate_lmdb_file(temporary, rows, dictionary, config.pocket.max_pocket_atoms, progress=progress)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    metadata = {
        "profile_name": profile, "filter_ast": {"column": "geometry_quality_tier", "operator": "in",
                                                  "values": list(definition["tiers"])},
        "schema_versions": {"sidecar": "1.0.0", "extraction": config.pipeline.extraction_version,
                            "lmdb_record": "1.0.0"},
        "source_sidecar_hashes": {
            "pockets.parquet": sha256_file(sidecars / "pockets.parquet"),
            "pocket_atoms.parquet": sha256_file(sidecars / "pocket_atoms.parquet"),
        },
        "record_count": len(rows), "pickle_protocol": 4, "lmdb_physical_sha256": _hash_output(destination, progress),
        "lmdb_logical_sha256": logical.hexdigest(), "map_size": map_size,
        "dictionary_sha256": contract.get("dictionary_sha256"), "task_sha256": contract.get("task_sha256"),
        "loader_sha256": contract.get("loader_sha256"), "helper_sha256": contract.get("helper_sha256"),
        "drugclip_library_contract_fingerprint": contract.get("drugclip_library_contract_fingerprint",
                                                               contract.get("drugclip_contract_fingerprint")),
        "cache_action": f"Namespace or invalidate BioSensIA-DC embedding cache by {logical.hexdigest()}",
    }
    profile_path = destination.with_suffix(destination.suffix + ".profile.json")
    atomic_write_bytes(profile_path, canonical_json_bytes(metadata) + b"\n")
    metadata["profile_metadata_sha256"] = sha256_file(profile_path)
    metadata["path"] = destination.as_posix()
    return metadata, rows


def validate_lmdb(run_dir: Path, profile: str, config: BuildConfig, *, progress: bool = True) -> list[str]:
    definition = PROFILE_FILTERS[profile]
    path = run_dir / "lmdb" / str(definition["filename"])
    rows = [row for row in read_sidecar(run_dir / "sidecars", "lmdb_records") if row["library_profile"] == profile]
    dictionary = read_dictionary(config.paths.drugclip_dictionary)
    return _validate_lmdb_file(path, rows, dictionary, config.pocket.max_pocket_atoms, progress=progress)


def _validate_lmdb_file(path, rows, dictionary, max_atoms, *, progress=True):
    errors: list[str] = []
    env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
    try:
        with env.begin() as transaction:
            count = transaction.stat()["entries"]
            if count != len(rows):
                errors.append(f"LMDB count {count} != sidecar count {len(rows)}")
            for index in track(range(count), description="Validating LMDB", total=count, enabled=progress):
                key = str(index).encode("ascii")
                value = transaction.get(key)
                if value is None:
                    errors.append(f"Missing dense LMDB key {index}")
                    continue
                try:
                    record = pickle.loads(value)
                    if set(record) != {"pocket", "pocket_atoms", "pocket_coordinates"}:
                        raise ValueError("wrong record keys")
                    _validate_record_fields(record["pocket"], record["pocket_atoms"],
                                            record["pocket_coordinates"], dictionary, max_atoms)
                    if rows and record["pocket"] != rows[index]["pocket_instance_id"]:
                        raise ValueError("record order differs from sidecar")
                except Exception as error:
                    errors.append(f"Invalid LMDB record {index}: {error}")
    finally:
        env.close()
    if errors:
        raise ValueError("; ".join(errors))
    return errors


def _validate_record_fields(name, atoms, coordinates, dictionary, max_atoms):
    if not isinstance(name, str) or not isinstance(atoms, list) or not isinstance(coordinates, np.ndarray):
        raise ValueError("wrong record field type")
    if coordinates.dtype != np.float32 or coordinates.dtype.str != "<f4" or not coordinates.flags.c_contiguous:
        raise ValueError("coordinates must be contiguous little-endian float32")
    if coordinates.shape != (len(atoms), 3) or not 1 <= len(atoms) <= max_atoms:
        raise ValueError("invalid atom/coordinate shape or count")
    if not np.isfinite(coordinates).all() or "H" in atoms:
        raise ValueError("nonfinite coordinates or hydrogen token")
    missing = set(atoms) - dictionary
    if missing:
        raise ValueError(f"dictionary lacks tokens {sorted(missing)}")
    transformed = [token[1] if token and token[0].isdigit() and len(token) > 1 else token[0] for token in atoms]
    if transformed != atoms:
        raise ValueError("DrugCLIP token transformation is lossy")


def _map_size(payloads, config):
    if isinstance(config.lmdb.map_size, int):
        return config.lmdb.map_size
    page = 4096
    # LMDB stores large values in page-aligned overflow pages; counting only
    # serialized bytes substantially underestimates the physical map.
    page_footprint = sum(
        ((len(value) + len(str(index)) + 32 + page - 1) // page) * page
        for index, (_, value) in enumerate(payloads)
    )
    branch_pages = ((len(payloads) + 99) // 100) * page
    raw = page_footprint + branch_pages + 1_048_576
    size = int(raw * (1 + config.lmdb.map_size_headroom_fraction))
    return max(16 * page, ((size + page - 1) // page) * page)


def _write_lmdb(path, payloads, map_size, *, auto, profile, progress):
    """Write atomically, growing and restarting only an automatically sized map."""
    attempt = 1
    while True:
        path.unlink(missing_ok=True)
        env = lmdb.open(str(path), subdir=False, map_size=map_size, lock=False,
                        sync=True, metasync=True, max_dbs=1)
        try:
            description = f"Exporting LMDB {profile}"
            if attempt > 1:
                description += f" (retry {attempt}, {map_size // (1024 * 1024)} MiB map)"
            with env.begin(write=True) as transaction:
                for index, (_, value) in enumerate(track(payloads, description=description,
                                                         total=len(payloads), enabled=progress)):
                    key = str(index).encode("ascii")
                    if not transaction.put(key, value, overwrite=False):
                        raise RuntimeError(f"Duplicate LMDB key {key!r}")
            env.sync()
            return map_size
        except lmdb.MapFullError as error:
            if not auto:
                raise lmdb.MapFullError(
                    f"Configured LMDB map_size={map_size} is too small; increase lmdb.map_size"
                ) from error
            map_size = max(map_size * 2, map_size + 16 * 1024 * 1024)
            attempt += 1
        finally:
            env.close()


def _hash_output(path: Path, progress: bool) -> str:
    from .progress import file_progress
    with file_progress(path, description=f"Hashing {path.name}",
                       enabled=progress and path.stat().st_size > 50_000_000) as bar:
        return sha256_file(path, progress=bar)
