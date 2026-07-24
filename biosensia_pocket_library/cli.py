"""Command-line interface for reproducible candidate-pocket library builds."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .config import BuildConfig, load_config
from .combine_set_pipeline import build_combine_set_library
from .combine_set_source import discover_combine_set, inventory_combine_set, select_combine_set
from .drugclip_contract import verify_drugclip_contract
from .exceptions import ConfigurationError, SourceIntegrityError
from .finalization import finalize_run
from .index_parser import parse_index
from .inventory import inventory_sources
from .lmdb_export import export_lmdb
from .manifest import write_manifest
from .pipeline import DISTRIBUTION_ID, _artifact_inventory, build_library
from .rcsb import download_mmcif_files
from .reporting import generate_reports
from .sidecars import read_sidecar, write_sidecars
from .validation import validate_run


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        return _dispatch(args)
    except (ConfigurationError, argparse.ArgumentError) as error:
        print(f"Configuration error: {error}", file=sys.stderr)
        return 2
    except SourceIntegrityError as error:
        print(f"Source-integrity error: {error}", file=sys.stderr)
        return 3
    except (ValueError, FileNotFoundError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    except Exception as error:
        print(f"Infrastructure error ({type(error).__name__}): {error}", file=sys.stderr)
        return 3


def _dispatch(args) -> int:
    if args.command in {"export-lmdb", "finalize", "validate", "report"}:
        run_dir = args.run_dir.resolve()
        config_path = args.config or run_dir / "config.resolved.toml"
        config = load_config(config_path)
    else:
        overrides = _overrides(args)
        config = load_config(args.config, overrides=overrides)
    progress = config.pipeline.progress and not getattr(args, "no_progress", False)
    if getattr(args, "offline", False) and getattr(args, "refresh_cache", False):
        raise ConfigurationError("--offline and --refresh-cache are incompatible")
    if args.command == "check-drugclip-contract":
        print(json.dumps(verify_drugclip_contract(config, progress=progress), indent=2, sort_keys=True))
        return 0
    if args.command == "parse-index":
        records, _, summary = parse_index(config.paths.index_dir / "INDEX_general_PL.2020R1.lst", DISTRIBUTION_ID)
        selected = _cli_select(records, args)
        summary["selected_count"] = len(selected)
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
    if args.command == "inventory":
        records, _, summary = parse_index(config.paths.index_dir / "INDEX_general_PL.2020R1.lst", DISTRIBUTION_ID)
        selected = _cli_select(records, args)
        files, directories, issues = inventory_sources(config.paths.index_dir, config.paths.complex_root,
                                                        selected, progress=progress)
        print(json.dumps({**summary, "selected_count": len(selected), "file_count": len(files),
                          "discovered_count": sum(bool(item["complex_directory"]) for item in directories.values()),
                          "issue_count": len(issues)}, indent=2, sort_keys=True))
        return 0
    if args.command == "inventory-combine-set":
        root = config.paths.combine_set_root
        assert root is not None
        discovered = discover_combine_set(root)
        selected = select_combine_set(
            discovered, args.pdb_id or _read_ids(args.pdb_ids_file), args.limit
        )
        files, directories, issues = inventory_combine_set(selected, config, progress=progress)
        print(json.dumps({"discovered_count": len(discovered), "selected_count": len(selected),
                          "file_count": len(files), "bundle_count": len(directories),
                          "issue_count": len(issues)}, indent=2, sort_keys=True))
        return 0
    if args.command == "download-rcsb":
        records, _, _ = parse_index(config.paths.index_dir / "INDEX_general_PL.2020R1.lst", DISTRIBUTION_ID)
        selected = _cli_select(records, args)
        paths, _, failures = download_mmcif_files([record.pdb_id for record in selected], config,
                                                  refresh=args.refresh_cache, progress=progress)
        print(json.dumps({"cached_count": len(paths), "failure_count": len(failures)}, indent=2))
        return 0
    if args.command in {"build", "build-sidecars"}:
        run_dir = build_library(config, pdb_ids=args.pdb_id or _read_ids(args.pdb_ids_file), limit=args.limit,
                                year_from=args.year_from, year_to=args.year_to, resume=args.resume,
                                overwrite_run=args.overwrite_run, export=args.command == "build", progress=progress)
        print(run_dir)
        return 0
    if args.command in {"build-combine-set", "build-combine-set-sidecars"}:
        run_dir = build_combine_set_library(
            config, pdb_ids=args.pdb_id or _read_ids(args.pdb_ids_file), limit=args.limit,
            resume=args.resume, overwrite_run=args.overwrite_run,
            export=args.command == "build-combine-set", progress=progress,
        )
        print(run_dir)
        return 0
    if args.command == "export-lmdb":
        contract = verify_drugclip_contract(config, progress=progress)
        metadata, lmdb_rows = export_lmdb(run_dir, config, contract, args.profile,
                                          overwrite=args.overwrite, progress=progress)
        rows = {name: read_sidecar(run_dir / "sidecars", name) for name in _table_names()}
        rows["lmdb_records"] = [row for row in rows["lmdb_records"] if row["library_profile"] != args.profile] + lmdb_rows
        sidecar_results = write_sidecars(run_dir / "sidecars", rows, progress=progress)
        manifest_path = run_dir / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest.setdefault("lmdb_profiles", {})[args.profile] = {
            **metadata, "path": Path(metadata["path"]).relative_to(run_dir).as_posix()
        }
        manifest["sidecar_artifacts"] = {
            name: {**item, "path": Path(item["path"]).relative_to(run_dir).as_posix()}
            for name, item in sidecar_results.items()
        }
        manifest["output_files"] = _artifact_inventory(run_dir, progress=progress)
        write_manifest(run_dir, manifest)
        print(json.dumps(metadata, indent=2, sort_keys=True))
        return 0
    if args.command == "finalize":
        manifest = finalize_run(run_dir, config, progress=progress)
        print(json.dumps({"run_id": manifest["run_id"], "status": manifest["status"],
                          "completed_at_utc": manifest["completed_at_utc"],
                          "counts": manifest["counts"]}, indent=2, sort_keys=True))
        return 0
    if args.command == "validate":
        errors = validate_run(run_dir, config, progress=progress)
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        print("Validation passed")
        return 0
    if args.command == "report":
        manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
        outputs = generate_reports(run_dir, manifest)
        manifest["output_files"] = _artifact_inventory(run_dir, progress=progress)
        write_manifest(run_dir, manifest)
        print(json.dumps(outputs, indent=2, sort_keys=True))
        return 0
    raise ConfigurationError(f"Unknown command {args.command}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="biosensia-pocket-library")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("check-drugclip-contract", "inventory", "parse-index", "download-rcsb",
                    "build-sidecars", "build", "inventory-combine-set",
                    "build-combine-set-sidecars", "build-combine-set"):
        item = subparsers.add_parser(command)
        _common(item, selection=command != "check-drugclip-contract",
                combine_set="combine-set" in command)
        if command in {"build", "build-sidecars", "build-combine-set", "build-combine-set-sidecars"}:
            item.add_argument("--resume", action="store_true")
            item.add_argument("--overwrite-run", action="store_true")
        if command == "download-rcsb":
            item.add_argument("--refresh-cache", action="store_true")
    export_parser = subparsers.add_parser("export-lmdb")
    export_parser.add_argument("--run-dir", type=Path, required=True)
    export_parser.add_argument("--config", type=Path)
    export_parser.add_argument("--profile", choices=("default", "tier-a", "tiers-ab", "all-usable"), default="default")
    export_parser.add_argument("--overwrite", action="store_true")
    export_parser.add_argument("--no-progress", action="store_true")
    for command in ("finalize", "validate", "report"):
        item = subparsers.add_parser(command)
        item.add_argument("--run-dir", type=Path, required=True)
        item.add_argument("--config", type=Path)
        item.add_argument("--no-progress", action="store_true")
    return parser


def _common(parser, *, selection, combine_set=False):
    parser.add_argument("--config", type=Path)
    if combine_set:
        parser.add_argument("--combine-set-root", type=Path)
        parser.add_argument("--trust-pickles", action="store_true")
    else:
        parser.add_argument("--index-dir", type=Path)
        parser.add_argument("--complex-root", type=Path)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument("--no-progress", action="store_true")
    if selection:
        parser.add_argument("--pdb-id", action="append")
        parser.add_argument("--pdb-ids-file", type=Path)
        parser.add_argument("--limit", type=int)
        if not combine_set:
            parser.add_argument("--year-from", type=int)
            parser.add_argument("--year-to", type=int)


def _overrides(args):
    values = {}
    for attribute, key in (("index_dir", "paths.index_dir"), ("complex_root", "paths.complex_root"),
                           ("combine_set_root", "paths.combine_set_root"),
                           ("workers", "pipeline.workers")):
        value = getattr(args, attribute, None)
        if value is not None:
            values[key] = value
    if getattr(args, "offline", False):
        values["pipeline.offline"] = True
    if getattr(args, "fail_fast", False):
        values["pipeline.fail_fast"] = True
    if getattr(args, "no_progress", False):
        values["pipeline.progress"] = False
    if getattr(args, "trust_pickles", False):
        values["combine_set.trusted_pickles"] = True
    return values


def _read_ids(path):
    return [line.strip().lower() for line in path.read_text().splitlines() if line.strip()] if path else None


def _cli_select(records, args):
    ids = args.pdb_id or _read_ids(args.pdb_ids_file)
    wanted = set(ids) if ids else None
    selected = [record for record in records if (wanted is None or record.pdb_id in wanted)
                and (args.year_from is None or record.release_year >= args.year_from)
                and (args.year_to is None or record.release_year <= args.year_to)]
    selected.sort(key=lambda record: record.pdb_id)
    return selected[:args.limit] if args.limit is not None else selected


def _table_names():
    from .schemas import TABLES
    return TABLES


if __name__ == "__main__":
    raise SystemExit(main())
