#!/usr/bin/env python3
"""Verify a produced LMDB through the actual linked DrugCLIP wrapper chain.

Run this with the BioSensIA-DC Python environment, which supplies PyTorch and
the compiled Uni-Core installation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from tqdm.auto import tqdm


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--biosensia-root", type=Path, required=True)
    parser.add_argument("--max-pocket-atoms", type=int, default=256)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    biosensia = args.biosensia_root.resolve()
    sys.path.insert(0, str(biosensia / "external/Uni-Core"))
    sys.path.insert(0, str(biosensia / "external/DrugCLIP"))
    bar = tqdm(total=5, desc="Verifying linked DrugCLIP loader", unit="stage")
    try:
        from unimol.data import (  # noqa: PLC0415 - paths are runtime inputs
            AffinityPocketDataset,
            CroppingPocketDataset,
            LMDBDataset,
            NormalizeDataset,
            RemoveHydrogenPocketDataset,
        )
        bar.update(1)

        lmdb_path = run_dir / "lmdb/candidate_pockets.lmdb"
        raw_dataset = LMDBDataset(str(lmdb_path))
        if len(raw_dataset) < 1:
            raise ValueError("The default LMDB contains no records")
        raw = raw_dataset[0]
        bar.update(1)
        affinity = AffinityPocketDataset(raw_dataset, 1, "pocket_atoms", "pocket_coordinates", False, "pocket")
        transformed = affinity[0]
        hydrogen_removed = RemoveHydrogenPocketDataset(
            affinity, "pocket_atoms", "pocket_coordinates", True, True
        )
        after_hydrogen = hydrogen_removed[0]
        cropped = CroppingPocketDataset(
            hydrogen_removed, 1, "pocket_atoms", "pocket_coordinates", args.max_pocket_atoms
        )
        after_crop = cropped[0]
        normalized = NormalizeDataset(cropped, "pocket_coordinates")[0]
        bar.update(1)

        original_atoms = np.asarray(raw["pocket_atoms"])
        original_coordinates = np.asarray(raw["pocket_coordinates"], dtype=np.float32)
        checks = {
            "loader_atom_transform_is_identity": bool(
                np.array_equal(original_atoms, transformed["pocket_atoms"])
            ),
            "hydrogen_wrapper_is_identity": bool(
                np.array_equal(transformed["pocket_atoms"], after_hydrogen["pocket_atoms"])
                and np.array_equal(
                    transformed["pocket_coordinates"], after_hydrogen["pocket_coordinates"]
                )
            ),
            "crop_wrapper_is_identity": bool(
                np.array_equal(after_hydrogen["pocket_atoms"], after_crop["pocket_atoms"])
                and np.array_equal(
                    after_hydrogen["pocket_coordinates"], after_crop["pocket_coordinates"]
                )
            ),
            "input_dtype_float32": original_coordinates.dtype == np.float32,
            "normalized_center_is_zero": bool(
                np.allclose(normalized["pocket_coordinates"].mean(axis=0), 0.0, atol=5e-5)
            ),
            "normalization_preserves_pairwise_distances": bool(
                np.allclose(
                    _pairwise(after_crop["pocket_coordinates"]),
                    _pairwise(normalized["pocket_coordinates"]),
                    rtol=1e-5,
                    atol=1e-5,
                )
            ),
            "pocket_name_preserved": transformed["pocket"] == raw["pocket"],
        }
        pockets = pq.read_table(run_dir / "sidecars/pockets.parquet", columns=["pocket_instance_id"]).column(0).to_pylist()
        checks["pocket_name_joins_sidecar"] = transformed["pocket"] in set(pockets)
        if not all(checks.values()):
            failed = [name for name, passed in checks.items() if not passed]
            raise AssertionError(f"DrugCLIP loader compatibility failed: {failed}")
        bar.update(1)

        manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
        artifact = {
            "status": "passed",
            "record_index": 0,
            "pocket_instance_id": transformed["pocket"],
            "atom_count": len(after_crop["pocket_atoms"]),
            "checks": checks,
            "biosensia_commit": manifest["drugclip_contract"].get("biosensia_commit"),
            "drugclip_contract_fingerprint": manifest["drugclip_contract_fingerprint"],
            "task_sha256": manifest["drugclip_contract"].get("task_sha256"),
            "loader_sha256": manifest["drugclip_contract"].get("loader_sha256"),
            "helper_sha256": manifest["drugclip_contract"].get("helper_sha256"),
            "dictionary_sha256": manifest["drugclip_contract"].get("dictionary_sha256"),
            "checkpoint_sha256": manifest["drugclip_contract"].get("checkpoint_sha256"),
            "checkpoint_encoding_status": "not_run_by_loader_contract_test",
        }
        destination = run_dir / "reports/drugclip_loader_integration.json"
        _atomic_text(destination, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
        digest = hashlib.sha256(destination.read_bytes()).hexdigest()
        manifest.setdefault("integration_tests", {})["drugclip_loader"] = {
            "status": "passed",
            "artifact": destination.relative_to(run_dir).as_posix(),
            "sha256": digest,
        }
        manifest["output_files"] = [item for item in manifest.get("output_files", [])
                                    if item["path"] != destination.relative_to(run_dir).as_posix()]
        manifest["output_files"].append({"path": destination.relative_to(run_dir).as_posix(),
                                         "size_bytes": destination.stat().st_size, "sha256": digest})
        _atomic_text(run_dir / "manifest.json", json.dumps(manifest, sort_keys=True, separators=(",", ":")) + "\n")
        bar.update(1)
        print(json.dumps(artifact, indent=2, sort_keys=True))
        return 0
    finally:
        bar.close()


def _pairwise(coordinates):
    coordinates = np.asarray(coordinates, dtype=np.float64)
    return np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=2)


def _atomic_text(path: Path, value: str) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(value, encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
