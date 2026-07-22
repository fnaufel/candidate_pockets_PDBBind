from pathlib import Path

from biosensia_pocket_library.cli import _parser


def test_export_lmdb_defaults_to_required_default_profile():
    args = _parser().parse_args(["export-lmdb", "--run-dir", str(Path("run"))])
    assert args.profile == "default"


def test_finalize_accepts_an_explicit_run_directory():
    args = _parser().parse_args(["finalize", "--run-dir", str(Path("run"))])
    assert args.command == "finalize"
