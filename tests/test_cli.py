from pathlib import Path

from biosensia_pocket_library.cli import _parser


def test_export_lmdb_defaults_to_required_default_profile():
    args = _parser().parse_args(["export-lmdb", "--run-dir", str(Path("run"))])
    assert args.profile == "default"


def test_finalize_accepts_an_explicit_run_directory():
    args = _parser().parse_args(["finalize", "--run-dir", str(Path("run"))])
    assert args.command == "finalize"


def test_combine_set_build_accepts_explicit_pickle_trust():
    args = _parser().parse_args(["build-combine-set", "--limit", "1", "--trust-pickles"])
    assert args.command == "build-combine-set"
    assert args.trust_pickles is True


def test_post_build_rcsb_workflow_commands_are_source_neutral():
    parser = _parser()
    plan = parser.parse_args([
        "plan-rcsb", "--run-dir", "run", "--output", "request.tsv",
    ])
    prefetch = parser.parse_args([
        "prefetch-rcsb", "--request", "request.tsv", "--cache-dir", "cache",
    ])
    enrich = parser.parse_args([
        "enrich-from-cache", "--run-dir", "run", "--cache-dir", "cache",
    ])
    assert plan.command == "plan-rcsb"
    assert prefetch.command == "prefetch-rcsb"
    assert enrich.command == "enrich-from-cache"
    assert not hasattr(enrich, "combine_set_root")
