"""Configuration loading, validation, and fingerprinting."""

from __future__ import annotations

import dataclasses
import os
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, get_args, get_origin, get_type_hints

from .constants import DEFAULT_MODIFIED_AMINO_ACIDS
from .exceptions import ConfigurationError
from .hashing import canonical_json_hash, sha256_file


@dataclass(slots=True)
class PipelineConfig:
    schema_version: str = "1.0.0"
    extraction_version: str = "1"
    random_seed: int = 1
    offline: bool = False
    workers: int = 1
    fail_fast: bool = False
    progress: bool = True


@dataclass(slots=True)
class PathsConfig:
    index_dir: Path = Path("data/raw/index")
    complex_root: Path = Path("data/raw/P-L")
    external_cache_dir: Path = Path("data/cache/external")
    output_root: Path = Path("data/processed/pdbbind_2020_v2024p_20250804")
    drugclip_dir: Path | None = None
    biosensia_root: Path | None = None
    drugclip_dictionary: Path | None = None
    drugclip_checkpoint: Path | None = None


@dataclass(slots=True)
class PocketConfig:
    distance_cutoff_angstrom: float = 6.0
    distance_uses_ligand_heavy_atoms: bool = True
    include_protein_hydrogens: bool = False
    include_allowlisted_polymer_hetatm: bool = True
    polymer_classification_policy: str = "pdbbind_local_v1"
    modified_residue_allowlist: tuple[str, ...] = DEFAULT_MODIFIED_AMINO_ACIDS
    max_pocket_atoms: int = 256
    deterministic_crop: bool = True
    minimum_pocket_atoms_hard: int = 1
    minimum_pocket_atoms_warning: int = 20


@dataclass(slots=True)
class StructureConfig:
    model_policy: str = "first"
    altloc_policy: str = "highest_occupancy"
    coordinate_match_tolerance_angstrom: float = 0.50
    strict_atom_match_tolerance_angstrom: float = 0.10
    alignment_minimum_ca_atoms: int = 10
    alignment_minimum_residue_coverage: float = 0.50
    alignment_maximum_rmsd_angstrom: float = 2.0
    probable_mapping_minimum_residue_coverage: float = 0.25
    mapping_candidate_tie_margin: float = 0.01


@dataclass(slots=True)
class ComparisonConfig:
    atom_jaccard_moderate_minimum: float = 0.70
    atom_jaccard_severe_minimum: float = 0.40
    residue_jaccard_moderate_minimum: float = 0.80
    residue_jaccard_severe_minimum: float = 0.50
    coordinate_rmsd_warning_angstrom: float = 0.10


@dataclass(slots=True)
class QualityConfig:
    rules_file: Path = Path("config/pocket-quality-rules.toml")
    covalent_radius_margin_angstrom: float = 0.40
    excluded_component_bridge_cutoff_angstrom: float = 3.0
    separated_component_cutoff_angstrom: float = 8.0
    nearby_nonprotein_cutoff_angstrom: float = 8.0


@dataclass(slots=True)
class ElementsConfig:
    explicit_mappings: dict[str, str] = field(default_factory=dict)
    unsupported_policy: str = "reject"


@dataclass(slots=True)
class LigandConfig:
    primary_format: str = "sdf"
    fallback_format: str = "mol2"
    sanitize: bool = True
    allow_multiple_components: bool = True
    multiple_sdf_record_policy: str = "reject"
    pocket_defining_component_policy: str = "all_components"
    require_3d_coordinates: bool = True


@dataclass(slots=True)
class RcsbConfig:
    download_mmcif: bool = True
    download_compressed: bool = True
    timeout_seconds: float = 60.0
    maximum_retries: int = 5
    requests_per_second: float = 2.0
    use_data_api: bool = True
    use_graphql: bool = True
    cache_mode: str = "content_addressed"


@dataclass(slots=True)
class BibliographyConfig:
    external_enrichment_enabled: bool = False
    extract_mmcif_citations: bool = True
    query_crossref: bool = False
    query_pubmed: bool = False
    allow_pdbbind_page_lookup: bool = False
    contact_email_env: str = "BIOSENSIA_BIBLIOGRAPHY_EMAIL"
    pubmed_api_key_env: str = "NCBI_API_KEY"
    manual_overrides_path: Path = Path("data/config/affinity_reference_overrides.csv")


@dataclass(slots=True)
class LmdbConfig:
    map_size: int | str = "auto"
    map_size_headroom_fraction: float = 0.25
    overwrite: bool = False
    pickle_protocol: int = 4
    include_geometry_quality_tiers: tuple[str, ...] = ("A", "B")


@dataclass(slots=True)
class BuildConfig:
    pipeline: PipelineConfig = field(default_factory=PipelineConfig)
    paths: PathsConfig = field(default_factory=PathsConfig)
    pocket: PocketConfig = field(default_factory=PocketConfig)
    structure: StructureConfig = field(default_factory=StructureConfig)
    comparison: ComparisonConfig = field(default_factory=ComparisonConfig)
    quality: QualityConfig = field(default_factory=QualityConfig)
    elements: ElementsConfig = field(default_factory=ElementsConfig)
    ligand: LigandConfig = field(default_factory=LigandConfig)
    rcsb: RcsbConfig = field(default_factory=RcsbConfig)
    bibliography: BibliographyConfig = field(default_factory=BibliographyConfig)
    lmdb: LmdbConfig = field(default_factory=LmdbConfig)
    project_root: Path = field(default=Path("."), init=False, repr=False)

    def as_dict(self) -> dict[str, Any]:
        value = dataclasses.asdict(self)
        value.pop("project_root", None)
        return _convert_paths(value, self.project_root)

    @property
    def semantic_hash(self) -> str:
        value = self.as_dict()
        value["pipeline"].pop("workers", None)
        value["pipeline"].pop("fail_fast", None)
        value["pipeline"].pop("progress", None)
        value.pop("paths", None)  # source and DrugCLIP contents have independent fingerprints
        value["quality"]["rules_sha256"] = (
            sha256_file(self.quality.rules_file) if self.quality.rules_file.is_file()
            else canonical_json_hash({"builtin_quality_rules": "v1"})
        )
        value["quality"].pop("rules_file", None)
        value["bibliography"]["manual_overrides_sha256"] = (
            sha256_file(self.bibliography.manual_overrides_path)
            if self.bibliography.manual_overrides_path.is_file() else None
        )
        value["bibliography"].pop("manual_overrides_path", None)
        for key in ("timeout_seconds", "maximum_retries", "requests_per_second", "cache_mode"):
            value["rcsb"].pop(key, None)
        for key in ("map_size", "map_size_headroom_fraction", "overwrite"):
            value["lmdb"].pop(key, None)
        return canonical_json_hash(value)

    @property
    def operational_hash(self) -> str:
        paths = self.as_dict()["paths"]
        # Encoder selection is recorded configuration, not a library operation.
        paths.pop("drugclip_checkpoint", None)
        return canonical_json_hash({
            "workers": self.pipeline.workers, "fail_fast": self.pipeline.fail_fast,
            "progress": self.pipeline.progress, "paths": paths,
            "rcsb_network": {"timeout_seconds": self.rcsb.timeout_seconds,
                             "maximum_retries": self.rcsb.maximum_retries,
                             "requests_per_second": self.rcsb.requests_per_second,
                             "cache_mode": self.rcsb.cache_mode},
            "lmdb_physical": {"map_size": self.lmdb.map_size,
                              "map_size_headroom_fraction": self.lmdb.map_size_headroom_fraction,
                              "overwrite": self.lmdb.overwrite},
        })


_SECTIONS = {item.name: get_type_hints(BuildConfig)[item.name]
             for item in dataclasses.fields(BuildConfig) if item.init}


def load_config(
    path: str | Path | None = None,
    overrides: dict[str, Any] | None = None,
    project_root: str | Path = ".",
) -> BuildConfig:
    root = Path(project_root).resolve()
    raw: dict[str, Any] = {}
    if path is not None:
        with Path(path).open("rb") as handle:
            raw = tomllib.load(handle)
    if overrides:
        _merge_dotted_overrides(raw, overrides)
    unknown = set(raw) - set(_SECTIONS)
    if unknown:
        raise ConfigurationError(f"Unknown configuration sections: {sorted(unknown)}")

    kwargs: dict[str, Any] = {}
    for section, section_type in _SECTIONS.items():
        values = raw.get(section, {})
        if not isinstance(values, dict):
            raise ConfigurationError(f"Configuration section {section!r} must be a table")
        kwargs[section] = _construct_dataclass(section_type, values, section)
    config = BuildConfig(**kwargs)
    config.project_root = root
    _resolve_paths(config, root)
    _validate(config)
    return config


def _construct_dataclass(cls: type[Any], values: dict[str, Any], section: str) -> Any:
    fields = {item.name: item for item in dataclasses.fields(cls)}
    unknown = set(values) - set(fields)
    if unknown:
        raise ConfigurationError(f"Unknown keys in [{section}]: {sorted(unknown)}")
    hints = get_type_hints(cls)
    converted: dict[str, Any] = {}
    for key, value in values.items():
        converted[key] = _coerce(value, hints[key])
    return cls(**converted)


def _coerce(value: Any, annotation: Any) -> Any:
    origin = get_origin(annotation)
    args = get_args(annotation)
    if annotation is Path:
        return Path(value)
    if origin is tuple:
        return tuple(value)
    if origin is None:
        return value
    if origin is dict:
        return dict(value)
    if type(None) in args:
        if value in (None, "auto"):
            return None
        non_none = next(arg for arg in args if arg is not type(None))
        return _coerce(value, non_none)
    if Path in args and isinstance(value, str):
        return Path(value)
    return value


def _resolve_paths(config: BuildConfig, root: Path) -> None:
    for item in dataclasses.fields(config.paths):
        value = getattr(config.paths, item.name)
        if isinstance(value, Path) and not value.is_absolute():
            setattr(config.paths, item.name, _resolve(root, value))

    if config.paths.drugclip_dir is None:
        candidates = (
            root / "data/DrugCLIP",
            root / "BioSensIA-DC/external/DrugCLIP",
            root / "external/DrugCLIP",
        )
        config.paths.drugclip_dir = next(
            (_absolute(candidate) for candidate in candidates if candidate.exists()),
            candidates[0],
        )
    else:
        config.paths.drugclip_dir = _resolve(root, config.paths.drugclip_dir)

    drugclip = config.paths.drugclip_dir
    if config.paths.biosensia_root is None:
        if drugclip.parent.name == "external":
            config.paths.biosensia_root = drugclip.parent.parent
        else:
            candidate = root / "BioSensIA-DC"
            config.paths.biosensia_root = _absolute(candidate) if candidate.exists() else None
    if config.paths.drugclip_dictionary is None:
        config.paths.drugclip_dictionary = drugclip / "data/dict_pkt.txt"
    elif not config.paths.drugclip_dictionary.is_absolute():
        config.paths.drugclip_dictionary = _resolve(root, config.paths.drugclip_dictionary)
    if config.paths.drugclip_checkpoint is None:
        config.paths.drugclip_checkpoint = drugclip / "checkpoint_best.pt"
    elif not config.paths.drugclip_checkpoint.is_absolute():
        config.paths.drugclip_checkpoint = _resolve(root, config.paths.drugclip_checkpoint)

    config.quality.rules_file = _resolve(root, config.quality.rules_file)
    config.bibliography.manual_overrides_path = _resolve(root, config.bibliography.manual_overrides_path)


def _resolve(root: Path, value: Path) -> Path:
    return _absolute(value if value.is_absolute() else root / value)


def _absolute(value: Path) -> Path:
    """Normalize dot segments without dereferencing configured symlinks."""
    return Path(os.path.abspath(value))


def _validate(config: BuildConfig) -> None:
    if config.pipeline.workers < 1:
        raise ConfigurationError("pipeline.workers must be at least 1")
    if config.pocket.distance_cutoff_angstrom <= 0:
        raise ConfigurationError("pocket.distance_cutoff_angstrom must be positive")
    if config.pocket.max_pocket_atoms < 1:
        raise ConfigurationError("pocket.max_pocket_atoms must be positive")
    if not (1 <= config.pocket.minimum_pocket_atoms_hard <= config.pocket.minimum_pocket_atoms_warning
            <= config.pocket.max_pocket_atoms):
        raise ConfigurationError("Pocket minimum atom thresholds must satisfy 1 <= hard <= warning <= max")
    if config.pocket.include_protein_hydrogens:
        raise ConfigurationError("Extraction schema version 1 requires protein hydrogens to be excluded")
    if not config.pocket.deterministic_crop:
        raise ConfigurationError("Extraction schema version 1 requires deterministic cropping")
    if config.structure.model_policy != "first":
        raise ConfigurationError("Extraction schema version 1 supports only structure.model_policy='first'")
    if config.structure.altloc_policy != "highest_occupancy":
        raise ConfigurationError("Extraction schema version 1 supports only highest_occupancy altloc policy")
    if config.ligand.multiple_sdf_record_policy not in {"reject", "first"}:
        raise ConfigurationError("Unsupported multiple_sdf_record_policy")
    if config.ligand.primary_format not in {"sdf", "mol2"} or config.ligand.fallback_format not in {"sdf", "mol2"}:
        raise ConfigurationError("Ligand primary/fallback formats must be sdf or mol2")
    if config.ligand.primary_format == config.ligand.fallback_format:
        raise ConfigurationError("Ligand primary and fallback formats must differ")
    if config.ligand.pocket_defining_component_policy != "all_components":
        raise ConfigurationError("Extraction schema version 1 requires all_components ligand policy")
    if not config.ligand.require_3d_coordinates:
        raise ConfigurationError("Extraction schema version 1 requires bound 3D coordinates")
    if config.elements.unsupported_policy not in {"reject", "exclude"}:
        raise ConfigurationError("elements.unsupported_policy must be reject or exclude")
    for name in ("atom_jaccard_moderate_minimum", "atom_jaccard_severe_minimum",
                 "residue_jaccard_moderate_minimum", "residue_jaccard_severe_minimum"):
        if not 0 <= getattr(config.comparison, name) <= 1:
            raise ConfigurationError(f"comparison.{name} must lie in [0,1]")
    if config.comparison.atom_jaccard_severe_minimum > config.comparison.atom_jaccard_moderate_minimum:
        raise ConfigurationError("Atom Jaccard severe threshold cannot exceed moderate threshold")
    if config.comparison.residue_jaccard_severe_minimum > config.comparison.residue_jaccard_moderate_minimum:
        raise ConfigurationError("Residue Jaccard severe threshold cannot exceed moderate threshold")
    if config.lmdb.pickle_protocol != 4:
        raise ConfigurationError("Schema version 1 requires pickle protocol 4")
    if not config.lmdb.include_geometry_quality_tiers or not set(config.lmdb.include_geometry_quality_tiers) <= {"A", "B", "C"}:
        raise ConfigurationError("LMDB geometry tiers must be a nonempty subset of A, B, C")
    if not 0 <= config.lmdb.map_size_headroom_fraction <= 10:
        raise ConfigurationError("LMDB map-size headroom must lie in [0,10]")
    if config.rcsb.maximum_retries < 1 or config.rcsb.requests_per_second < 0:
        raise ConfigurationError("RCSB retries must be positive and request rate nonnegative")
    if config.rcsb.cache_mode != "content_addressed":
        raise ConfigurationError("Only content_addressed RCSB cache mode is supported")
    if not config.bibliography.external_enrichment_enabled and any((
        config.bibliography.query_crossref, config.bibliography.query_pubmed,
        config.bibliography.allow_pdbbind_page_lookup,
    )):
        raise ConfigurationError("External bibliography sources require external_enrichment_enabled=true")


def _merge_dotted_overrides(raw: dict[str, Any], overrides: dict[str, Any]) -> None:
    for dotted, value in overrides.items():
        section, separator, key = dotted.partition(".")
        if not separator:
            raise ConfigurationError(f"Override must be section.key: {dotted}")
        raw.setdefault(section, {})[key] = value


def _convert_paths(value: Any, root: Path) -> Any:
    if isinstance(value, dict):
        return {key: _convert_paths(item, root) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_convert_paths(item, root) for item in value]
    if isinstance(value, Path):
        try:
            return value.relative_to(root).as_posix()
        except ValueError:
            return value.as_posix()
    return value
