# Specification: Build a BioSensIA-DC Candidate Pocket Library from the PDBbind 2020/v2024-Reprocessed Special Distribution

## 1. Objective

Implement a reproducible pipeline that builds a candidate protein-pocket library from the special PDBbind distribution updated on 2025-08-04. This distribution contains the 19,037 protein-ligand complexes from PDBbind v2020 that remain in PDBbind v2024, uses structures reprocessed with the PDBbind v2024 workflow, and includes later index corrections. It must not be described as an unmodified copy of the original PDBbind 2020R1 release.

The pipeline must:

1. Discover and validate the PDBBind input files.
2. Parse the PDBBind protein–ligand index and binding measurements.
3. Parse and validate each bound ligand.
4. Re-extract a pocket reproducibly from the PDBBind processed protein and ligand coordinates.
5. Compare the re-extracted pocket with the PDBBind-provided `*_pocket.pdb` file.
6. Identify which protein chains contribute atoms to the re-extracted pocket.
7. Download and cache the corresponding RCSB PDBx/mmCIF entry.
8. Map PDBBind chain and residue identifiers to RCSB chains, polymer entities, and UniProt accessions.
9. Enrich citations and bibliographic metadata without retaining the legacy `*.pdf` reference token from the PDBBind index.
10. Store normalized metadata in Parquet sidecars.
11. Export one or more minimal DrugCLIP-compatible pocket LMDB files.
12. Produce machine-readable validation reports, processing statuses, warnings, and checksums.

The provenance chain for every accepted pocket must be traceable:

```text
PDBbind distribution identity and source-file checksums
  → PDBBind index record
  → PDBBind complex directory
  → selected ligand file and ligand coordinates
  → PDBBind processed protein coordinates
  → extraction configuration
  → selected pocket atoms
  → contributing PDB chains
  → RCSB chain instances and polymer entities
  → UniProt accessions
  → DrugCLIP LMDB record
```

The PDBbind v2024-reprocessed files in this special distribution are the sole source of the geometry supplied to DrugCLIP. RCSB data is used only for identity mapping, annotations, alignment, validation, and bibliographic enrichment. RCSB availability or improved mapping must never change selected PDBbind atoms, coordinates, geometry-content hashes, or pocket identifiers.

---

## 2. Scope

### 2.1 Dataset identity and included input set

Record these distinct identity fields:

```text
dataset_family = "PDBbind"
nominal_complex_set_version = "2020"
structure_processing_version = "2024"
index_revision_date = "2025-08-04"
distribution_label = "PDBbind v2020/v2024-reprocessed special distribution"
```

Parse the date and descriptive wording from `data/raw/index/README` and the index headers when possible. Source-derived values and checksums take precedence over display defaults. Derive a stable `distribution_id` from the normalized identity fields and SHA-256 checksums of `README` and all four index files. Use a short form of `distribution_id` in identifier namespaces; `pdbbind-2020r1` alone is insufficient to distinguish this package from the original 2020R1 distribution.

The candidate library is built from the protein–ligand index:

```text
data/raw/index/INDEX_general_PL.2020R1.lst
```

The other index files must be inventoried and included in the manifest, but they are not used to build the initial protein-pocket library:

```text
data/raw/index/INDEX_general_NL.2020R1.lst
data/raw/index/INDEX_general_PN.2020R1.lst
data/raw/index/INDEX_general_PP.2020R1.lst
```

The complex directories are below:

```text
data/raw/P-L/1981-2000/<pdb_id>/
data/raw/P-L/2001-2010/<pdb_id>/
data/raw/P-L/2011-2019/<pdb_id>/
```

A normal complex directory contains:

```text
<pdb_id>_ligand.sdf
<pdb_id>_ligand.mol2
<pdb_id>_pocket.pdb
<pdb_id>_protein.pdb
```

### 2.2 Non-goals for the first implementation

The initial implementation must not:

* modify or overwrite any file under `data/raw`;
* regenerate or optimize ligand conformations;
* perform molecular docking;
* predict binding affinity;
* use the reference ligand’s binding measurement to modify the DrugCLIP ranking;
* aggregate pockets into proteins during library construction;
* remove structural redundancy automatically;
* use RCSB data to decide which PDBbind atoms belong to the exported geometry;
* use RCSB coordinates instead of the PDBbind processed geometry;
* assume that the primary PDB citation is the source of the PDBBind affinity measurement;
* store the PDBBind legacy PDF filename token;
* commit PDBbind raw files, derived coordinates, coordinate-bearing sidecars, LMDB files, caches, or run outputs to Git because their distribution may be license-restricted or they may be too large.

---

## 3. Grounding constraints and verified DrugCLIP contract

Ground DrugCLIP and BioSensIA-DC behavior in the checkout at `https://github.com/fnaufel/BioSensIA-DC`. The current project exposes its `external/DrugCLIP` directory through the local symbolic link `data/DrugCLIP`. Permit a configured equivalent, but record both the lexical link path and resolved target.

Before extraction or export, a `check-drugclip-contract` stage must resolve the link and record:

* the resolved checkout path;
* the checkout Git commit and dirty state;
* SHA-256 checksums of the active task, dataset wrappers, dictionary, and BioSensIA LMDB helper;
* the configured maximum pocket size, dictionary path, optional checkpoint path, and coordinate-normalization behavior.

A missing or broken DrugCLIP link blocks compatibility validation and LMDB export, but not PDBbind-only inventory and geometry sidecars. Derive the BioSensIA-DC root from the resolved target when possible and verify `lmdb_helpers.py`; otherwise require `biosensia_root`. The checkpoint is not read or hashed during candidate-library construction: it cannot change pocket geometry, sidecars, or LMDB bytes. A missing checkpoint blocks only the separate encoder integration test or embedding generation.

Separate the **library contract** from the **encoder contract**. The versioned library contract covers only the task, loader, helper, and dictionary files that define or validate the exported interface. Record the linked checkout commit and dirty state as provenance, but do not include them in the library fingerprint because unrelated checkout changes cannot affect the library when the relevant files have identical hashes. The encoder contract combines a completed library or LMDB logical fingerprint with the exact checkpoint hash and must be created only by a command that actually loads, validates, or uses that checkpoint.

The configured checkpoint path may remain in resolved configuration for a later encoder command, but it participates in neither the library semantic hash nor the library operational hash.

The specification-review baseline is BioSensIA-DC commit `01c79cedb37379cf4d70beb9eb309fdf75518bd5`. The local linked dictionary at that revision contains special tokens plus `C`, `N`, `O`, `S`, and `H`. Runtime validation must still use the actual linked revision and hashes rather than assuming this baseline remains current.

Inspection established the current contract:

* retrieval records require `pocket`, `pocket_atoms`, and `pocket_coordinates`;
* values are pickled dictionaries under numeric ASCII LMDB keys;
* `AffinityPocketDataset` converts each supplied atom token to its first character, or its second character when it starts with a digit;
* `RemoveHydrogenPocketDataset` removes tokens equal to `H`;
* `CroppingPocketDataset` performs seeded, sampled, center-weighted cropping above `max_pocket_atoms`, whose default is 256;
* `NormalizeDataset` subtracts the pocket centroid in the loader;
* the loader casts coordinates to `float32`.

The compatibility stage must verify these facts against the linked checkout. Contract drift is fatal for LMDB export until a versioned adapter or specification update handles it.

Export canonical one-letter element tokens supported by the active dictionary. Reject by default any element that cannot be represented losslessly by the first-character loader behavior; for example, never pass `Se` and let the loader silently reinterpret it as `S`.

Remove hydrogens and pre-crop to at most the active maximum so the stochastic loader crop is never reached. Store uncentered PDBbind coordinates because centering occurs in the loader. The integration test must prove that hydrogen removal and cropping do not change an exported record.

The current helper uses dense numeric keys, `subdir=False`, `lock=False`, and pickled dictionaries. Reuse it only where it matches this specification; otherwise preserve the convention while adding atomic output. Fix serialization at pickle protocol 4 for schema version 1.

RCSB distinguishes PDB entries, polymer entities, non-polymer entities, polymer entity instances or chains, and chemical components. The Data API follows this hierarchy and supplements mmCIF data with external annotations.

RCSB distinguishes archive-assigned `label_asym_id` chain identifiers from author-supplied `auth_asym_id` identifiers. Legacy PDB files generally expose author-style chain identifiers, while mmCIF represents both systems.

RCSB distributes PDBx/mmCIF files through scripted HTTPS downloads.

---

## 4. Recommended package structure

Create a dedicated package rather than adding the entire implementation to `biosensia_target_fishing.py`.

```text
biosensia_pocket_library/
├── __init__.py
├── cli.py
├── config.py
├── constants.py
├── manifest.py
├── index_parser.py
├── source_inventory.py
├── ligand_parser.py
├── protein_parser.py
├── pocket_extractor.py
├── pocket_comparison.py
├── drugclip_compatibility.py
├── rcsb_download.py
├── rcsb_metadata.py
├── structure_mapping.py
├── citation_enrichment.py
├── quality.py
├── sidecars.py
├── schemas.py
├── lmdb_export.py
├── validation.py
├── hashing.py
├── models.py
└── exceptions.py
```

Add an executable entry point or script:

```text
scripts/build_pdbbind_pocket_library.py
```

The pipeline must expose both:

* a Python API suitable for tests and notebooks;
* a command-line interface suitable for batch jobs and Codex-driven execution.

---

## 5. Dependencies

Retain the project’s Python constraint:

```text
Python >= 3.11, < 3.12
```

Reuse existing dependencies where possible:

```text
numpy
polars
pyarrow
rdkit
lmdb
tqdm
pytest
```

Add:

```text
gemmi
httpx
```

Optional additions are permitted only when their use is justified and tested:

```text
tenacity
```

Use `gemmi` for legacy PDB and PDBx/mmCIF parsing. Use RDKit for ligand chemistry. Use `httpx` for cached HTTP requests. Use Polars and PyArrow for sidecar generation.

All new dependencies must be pinned or bounded in `pyproject.toml`. Record exact resolved versions in the manifest. Do not add `rapidfuzz` unless an implemented, reviewed matching rule actually requires it; identifier and atom matching must not depend on unconstrained fuzzy text similarity.

---

## 6. Directory layout

The pipeline must use the following default layout:

```text
data/
├── raw/
│   ├── index/
│   │   ├── INDEX_general_NL.2020R1.lst
│   │   ├── INDEX_general_PL.2020R1.lst
│   │   ├── INDEX_general_PN.2020R1.lst
│   │   ├── INDEX_general_PP.2020R1.lst
│   │   └── README
│   └── P-L/
│       ├── 1981-2000/
│       ├── 2001-2010/
│       └── 2011-2019/
├── cache/
│   └── external/
│       ├── objects/
│       │   └── sha256/<first_two_hex>/<sha256>
│       ├── request_index/
│       │   ├── rcsb_mmcif/
│       │   ├── rcsb_data_api/
│       │   ├── rcsb_graphql/
│       │   ├── crossref/
│       │   ├── pubmed/
│       │   └── pdbbind_pages/
│       └── metadata/
│           └── <request_hash>.json
└── processed/
    └── pdbbind_2020_v2024p_20250804/
        └── <run_id>/
            ├── manifest.json
            ├── config.resolved.toml
            ├── sidecars/
            ├── lmdb/
            ├── reports/
            ├── logs/
            └── checkpoints/
```

The `<run_id>` must include short semantic-configuration, source-fingerprint, selection-fingerprint, and DrugCLIP-library-contract hashes:

```text
pb20-v24p-20250804-v1-<semantic_config_hash_8>-<source_hash_8>-<selection_hash_8>-<contract_hash_8>
```

Do not use only a timestamp. A timestamp may be appended for readability. The run must refuse to reuse an existing directory if its recorded full hashes differ.

`data/raw/`, `data/cache/`, `data/processed/`, and the machine-specific `data/DrugCLIP` symlink must be ignored by Git. Documentation must state that derived coordinate-bearing artifacts inherit distribution concerns from the licensed source and must not be redistributed without authorization.

---

## 7. Configuration

Support an optional TOML configuration file and command-line overrides.

Example:

```toml
[pipeline]
schema_version = "1.0.0"
extraction_version = "1"
random_seed = 1
offline = false
workers = 8
fail_fast = false

[paths]
index_dir = "data/raw/index"
complex_root = "data/raw/P-L"
external_cache_dir = "data/cache/external"
output_root = "data/processed/pdbbind_2020_v2024p_20250804"
drugclip_dir = "data/DrugCLIP"
biosensia_root = "auto"
drugclip_dictionary = "data/DrugCLIP/data/dict_pkt.txt"
drugclip_checkpoint = "data/DrugCLIP/checkpoint_best.pt"

[pocket]
distance_cutoff_angstrom = 6.0
distance_uses_ligand_heavy_atoms = true
include_protein_hydrogens = false
include_allowlisted_polymer_hetatm = true
polymer_classification_policy = "pdbbind_local_v1"
modified_residue_allowlist = ["MSE"]
max_pocket_atoms = 256
deterministic_crop = true
minimum_pocket_atoms_hard = 1
minimum_pocket_atoms_warning = 20

[structure]
model_policy = "first"
altloc_policy = "highest_occupancy"
coordinate_match_tolerance_angstrom = 0.50
strict_atom_match_tolerance_angstrom = 0.10
alignment_minimum_ca_atoms = 10
alignment_minimum_residue_coverage = 0.50
alignment_maximum_rmsd_angstrom = 2.0
probable_mapping_minimum_residue_coverage = 0.25
mapping_candidate_tie_margin = 0.01

[comparison]
atom_jaccard_moderate_minimum = 0.70
atom_jaccard_severe_minimum = 0.40
residue_jaccard_moderate_minimum = 0.80
residue_jaccard_severe_minimum = 0.50
coordinate_rmsd_warning_angstrom = 0.10

[quality]
rules_file = "config/pocket-quality-rules.toml"
covalent_radius_margin_angstrom = 0.40
excluded_component_bridge_cutoff_angstrom = 3.0
separated_component_cutoff_angstrom = 8.0
nearby_nonprotein_cutoff_angstrom = 8.0

[elements]
explicit_mappings = {}
unsupported_policy = "reject"

[ligand]
primary_format = "sdf"
fallback_format = "mol2"
sanitize = true
allow_multiple_components = true
multiple_sdf_record_policy = "reject"
pocket_defining_component_policy = "all_components"
require_3d_coordinates = true

[rcsb]
download_mmcif = true
download_compressed = true
timeout_seconds = 60
maximum_retries = 5
requests_per_second = 2.0
use_data_api = true
use_graphql = true
cache_mode = "content_addressed"

[bibliography]
external_enrichment_enabled = false
extract_mmcif_citations = true
query_crossref = false
query_pubmed = false
allow_pdbbind_page_lookup = false
contact_email_env = "BIOSENSIA_BIBLIOGRAPHY_EMAIL"
pubmed_api_key_env = "NCBI_API_KEY"

[lmdb]
map_size = "auto"
map_size_headroom_fraction = 0.25
overwrite = false
pickle_protocol = 4
include_geometry_quality_tiers = ["A", "B"]
```

Configuration precedence is defaults, then TOML, then CLI overrides. Unknown keys are errors. Semantic configuration, operational configuration, and secrets must be distinguished. Worker count and progress display are operational; cutoff, classification, crop, filters, and offline/enrichment policies are semantic. Environment-provided secret values must never be written to resolved configuration, manifests, logs, hashes, or reports; record only the environment-variable names and whether values were available.

Compute hashes from canonical UTF-8 JSON with recursively sorted object keys, preserved list order, normalized POSIX-style project-relative paths, finite JSON numbers, and no insignificant whitespace. Reject NaN, infinity, duplicate TOML keys, and paths that cannot be normalized safely. Define:

* `semantic_config_hash` over all nonsecret values capable of changing scientific or exported contents, including enrichment enabled/offline policy;
* `operational_config_hash` over execution-only values such as workers and progress display;
* `source_fingerprint` over `distribution_id` plus every selected local source-file ID and checksum;
* `selection_fingerprint` over the canonical sorted selected PDB-ID list and the filter specification that produced it;
* `drugclip_library_contract_fingerprint` over a version tag and the task/loader/helper/dictionary hashes. The legacy alias `drugclip_contract_fingerprint` may be retained for readers, but neither the checkpoint nor whole-checkout revision/dirty state participates.

Changing only the configured checkpoint, or making unrelated changes elsewhere in the linked BioSensIA-DC checkout, must not change the candidate-library run ID. Relevant task/loader/helper/dictionary changes must change the versioned library-contract fingerprint. Checkpoint hashes belong to downstream encoder/embedding artifacts together with the LMDB logical checksum.

Backward compatibility is mandatory. A directory whose manifest predates the versioned library contract may contain only `drugclip_contract_fingerprint` and may embed `checkpoint_sha256` in `drugclip_contract` or LMDB profile metadata. Treat those fields as an opaque legacy run identity and retained provenance: do not rename the directory, rewrite its identity, require the checkpoint to validate it, or reject it because the new field is absent. Commands taking an explicit `--run-dir` (`validate`, `report`, loader verification, and `export-lmdb`) must continue to process such directories. New or re-exported profile metadata may use the new library-contract fields without mutating the legacy top-level identity. New runs must never recreate a legacy checkpoint-coupled identity.

Operational changes may reuse scientific checkpoints when all stage dependency hashes match. Semantic, source, or selection changes may not.

Every nonsecret resolved configuration value must be copied into:

```text
config.resolved.toml
manifest.json
```

---

## 8. Pipeline execution model

The pipeline must consist of resumable stages with an explicit dependency graph.

```text
bootstrap-identity
check-drugclip-contract
inventory
parse-index
parse-structures
extract-pockets
compare-pockets
download-rcsb
map-structures
enrich-citations
quality-control
write-sidecars
export-lmdb
validate
report
```

`bootstrap-identity` is a read-only preflight that resolves configuration and the DrugCLIP link, parses enough index metadata to determine selection, inventories and hashes the selected local source files, computes the four full fingerprints, and only then creates the run directory and initial manifest. The resumable `inventory` and `parse-index` stages reproduce and validate those bootstrap results as authoritative sidecars. Bootstrap failure must not leave a run directory.

Each stage must:

1. Read immutable inputs or outputs from preceding stages.
2. Write its result atomically.
3. Produce a completion marker containing the stage name and version, relevant semantic-configuration hash, input and upstream artifact hashes, source fingerprint, code commit and dirty-state fingerprint, and output hashes.
4. Be safely restartable.
5. Refuse to reuse a checkpoint when any marker field is incompatible or an output is missing or fails its recorded logical checksum.
6. Record failures as data instead of terminating the entire run, unless `--fail-fast` is specified.

Intermediate checkpoints may use Parquet files under:

```text
checkpoints/<stage-name>/
```

Temporary files must be written under the run directory and renamed atomically only after successful completion. Directory-valued stages must build a sibling temporary directory and atomically rename it. A failure after partial per-complex work must still produce isolated issue records without publishing an incomplete stage marker.

The geometry stages must depend only on PDBbind inputs, local geometry policies, and the verified DrugCLIP token contract. RCSB and bibliographic stages may enrich or assess geometry but must not feed back into extraction, content hashes, derivation identifiers, or geometry acceptance.

---

## 9. Source inventory and manifest

### 9.1 Manifest requirements

Create `manifest.json` before processing begins and update it after every completed stage.

It must contain:

```json
{
  "manifest_schema_version": "1.0.0",
  "pipeline_name": "biosensia-pdbbind-pocket-library",
  "pipeline_version": "...",
  "git_commit": "...",
  "git_dirty": false,
  "run_id": "...",
  "semantic_config_hash": "...",
  "operational_config_hash": "...",
  "source_fingerprint": "...",
  "selection_fingerprint": "...",
  "drugclip_contract_fingerprint": "...",
  "drugclip_library_contract_fingerprint": "...",
  "started_at_utc": "...",
  "completed_at_utc": null,
  "status": "running",
  "python_version": "...",
  "platform": "...",
  "dependency_versions": {},
  "drugclip_contract": {
    "link_path": "data/DrugCLIP",
    "resolved_revision": "...",
    "dirty": false,
    "task_sha256": "...",
    "loader_sha256": "...",
    "helper_sha256": "...",
    "dictionary_sha256": "...",
    "encoder_checkpoint": {
      "configured_path": "data/DrugCLIP/checkpoint_best.pt",
      "sha256": null,
      "verification_status": "not_evaluated"
    }
  },
  "dataset": {
    "name": "PDBbind",
    "distribution_id": "...",
    "distribution_label": "PDBbind v2020/v2024-reprocessed special distribution",
    "nominal_complex_set_version": "2020",
    "structure_processing_version": "2024",
    "index_revision_date": "2025-08-04",
    "source_readme_sha256": "...",
    "index_file": "data/raw/index/INDEX_general_PL.2020R1.lst",
    "index_declared_complex_count": null,
    "index_parsed_complex_count": null,
    "discovered_complex_directory_count": null
  },
  "configuration": {},
  "stage_statuses": {},
  "counts": {},
  "output_files": []
}
```

Identity fields should be read from `README` and index headers when possible, not hardcoded. Record disagreements as fatal source-identity issues. `git_commit` and `git_dirty` may be null when the project is not in a Git checkout, but that condition must be explicit. Record the linked BioSensIA-DC/DrugCLIP revision separately.

### 9.2 Source file inventory

Create:

```text
sidecars/source_files.parquet
```

One row per relevant source or downloaded file.

Required columns:

| Column               | Type               | Description                                                                                                                             |
| -------------------- | ------------------ | --------------------------------------------------------------------------------------------------------------------------------------- |
| `source_file_id`     | string             | Stable identifier derived from `source_kind`, normalized relative path, and SHA-256; a changed file receives a new ID                   |
| `source_kind`        | string             | `pdbbind_index`, `pdbbind_protein`, `pdbbind_ligand_sdf`, `pdbbind_ligand_mol2`, `pdbbind_pocket`, `rcsb_mmcif`, `rcsb_api_cache`, etc. |
| `pdb_id`             | nullable string    | Associated PDB ID                                                                                                                       |
| `path`               | string             | Path relative to the project root when possible                                                                                         |
| `size_bytes`         | integer            | File size                                                                                                                               |
| `sha256`             | string             | SHA-256 digest                                                                                                                          |
| `modified_time_utc`  | nullable timestamp | Filesystem modification time                                                                                                            |
| `download_url`       | nullable string    | Source URL for downloaded files                                                                                                         |
| `downloaded_at_utc`  | nullable timestamp | Download time                                                                                                                           |
| `http_etag`          | nullable string    | HTTP ETag                                                                                                                               |
| `http_last_modified` | nullable string    | HTTP Last-Modified                                                                                                                      |
| `validation_status`  | string             | `valid`, `missing`, `empty`, `checksum_mismatch`, `unreadable`                                                                          |
| `warning_codes`      | list of strings    | Nonfatal issues                                                                                                                         |

Inventory every expected PDBbind file for indexed complexes, including missing-file rows, plus discovered extras and every downloaded response used by the run. Do not follow arbitrary symlinks while inventorying raw data. Normalize paths relative to the project root without recording machine-specific absolute paths.

### 9.3 Checksum boundary

The checksum claim applies to:

* `README`, all four index files, and every expected or actually consumed complex file;
* every cached external response actually used by the run;
* `config.resolved.toml`, final Parquet sidecars, final report data, LMDB data files, and adjacent LMDB profile metadata.

It excludes the checksum inventory itself, `manifest.json`, mutable logs, checkpoints, temporary files, directory entries, and transient lock files. The manifest may refer to `source_files.parquet` instead of embedding all complex-file hashes. Final artifact checksums must be written only after artifacts are closed and immutable.

For LMDB, record both the physical file SHA-256 and an authoritative logical SHA-256 over the ordered sequence of length-framed `(key, serialized_value)` byte pairs. Reproducibility assertions use the logical digest when LMDB library or filesystem details make physical bytes differ.

---

## 10. PDBBind index parsing

### 10.1 Input

Parse:

```text
data/raw/index/INDEX_general_PL.2020R1.lst
```

Ignore blank lines and lines whose first non-space character is `#`.

### 10.2 Parsed fields

For each data line, extract:

* PDB ID;
* resolution or method token;
* release year;
* raw binding measurement;
* ligand label or ligand name when present;
* curator comments after the ligand label;
* original line number;
* a redacted reconstruction of the original line;
* SHA-256 of the original source line bytes for provenance without retaining its text in a table.

Before parsing the ligand label or comment, remove every PDF-looking token in the text after `//`, case-insensitively. For version 1, split into maximal non-whitespace tokens, peel surrounding punctuation for detection, and remove the entire original token when the basename after the last `/` or `\` ends in `.pdf`; also remove an immediately adjacent separator punctuation mark left with no content. Collapse whitespace left by removal. No original unredacted line may be persisted in a table, manifest field, report, exception, or log. The parser and global structured-logging scrubber may recognize such tokens transiently only to remove them. Tests must cover multiple tokens, mixed case, path separators, punctuation, and a PDF-looking token later in the curator comment.

After redaction, treat the first balanced parenthesized group as `ligand_label` when present; retain its interior verbatim after trimming. The remaining post-label text is `index_comment`. If parentheses are unbalanced, retain all redacted text as the comment and emit a warning. Measurements mentioned in comments are not additional authoritative measurements unless a separately versioned parser and evidence rule explicitly enables them.

### 10.3 Complex identifier

Use:

```text
complex_id = "pb20v24p-<distribution_hash_8>:<lowercase_pdb_id>"
```

Example:

```text
pb20v24p-a1b2c3d4:2l3r
```

The parser must reject conflicting duplicate PDB IDs. Byte-for-byte-equivalent duplicates are collapsed into one complex, with all source line numbers retained in a separate `index_record_occurrences.parquet`; counts must distinguish physical data lines from unique complexes.

### 10.4 Resolution fields

Store:

* `resolution_raw`;
* `resolution_angstrom`, when numeric;
* `experimental_method_hint`, such as `NMR` when indicated by the index.

Do not interpret `NMR` as a numeric resolution.

### 10.5 Binding-measurement parser

Parse forms including:

```text
Kd=49uM
Ki=0.43uM
Kd<10uM
IC50=5nM
Ka=1.2e6M-1
```

Support:

* measurement types such as `Kd`, `Ki`, `IC50`, `EC50`, `Ka`, and unknown types;
* relations `=`, `<`, `>`, `<=`, `>=`, and approximate values;
* concentration units `M`, `mM`, `uM`, `µM`, `nM`, `pM`, and `fM`;
* inverse-concentration `Ka` units such as `M-1`, `M^-1`, and `uM^-1`;
* whitespace variants;
* scientific notation;
* unparseable raw values.

Normalize dissociation, inhibition, and response concentration measurements to molar units:

$$
v_{\mathrm{M}} = v \times f_{\mathrm{unit}}
$$

When a positive normalized molar value exists, calculate:

$$
pX = -\log_{10}(v_{\mathrm{M}})
$$

The relation must be reversed when transformed to the `-log10` concentration scale:

```text
Kd < 10 µM  →  pKd > 5
Kd > 10 µM  →  pKd < 5
```

Never replace the raw measurement string with the normalized value.

`Ka` is an association constant and must not be normalized as a concentration. Normalize it to `M^-1` and, when positive, calculate `pKa_association = log10(Ka_M_inverse)` without reversing the relation. Keep this field semantically distinct from acid-dissociation pKa. Unknown types may be lexed into raw fields but must not receive a normalized value or p-value. Approximate markers include `~`, `≈`, and a configured textual prefix; uncertainty ranges and values embedded only in comments remain raw comment text in version 1.

`measurement_id` is a stable hash of `complex_id`, measurement ordinal, and the exact primary measurement token. Floating-point overflow, zero, negative values, malformed exponents, and unsupported dimensional units must produce explicit parse statuses rather than exceptions.

---

## 11. Binding-measurement sidecar

Create:

```text
sidecars/binding_measurements.parquet
```

Use one row per measurement rather than embedding the measurement directly into `complexes.parquet`. The current index may provide one measurement per complex, but the data model must support future multiple measurements.

Required columns:

| Column                        | Type             |
| ----------------------------- | ---------------- |
| `measurement_id`              | string           |
| `complex_id`                  | string           |
| `pdb_id`                      | string           |
| `measurement_type_raw`        | string           |
| `measurement_type_normalized` | nullable string  |
| `relation_raw`                | nullable string  |
| `relation_normalized`         | nullable string  |
| `value_raw`                   | nullable string  |
| `value_numeric`               | nullable float64 |
| `unit_raw`                    | nullable string  |
| `unit_normalized`             | nullable string  |
| `value_molar`                 | nullable float64 |
| `value_inverse_molar`         | nullable float64 |
| `p_measurement_name`          | nullable string  |
| `p_relation`                  | nullable string  |
| `p_value`                     | nullable float64 |
| `normalization_kind`          | nullable string  |
| `measurement_raw`             | string           |
| `parse_status`                | string           |
| `parse_warning_codes`         | list of strings  |
| `source_index_line_number`    | integer          |

Allowed `parse_status` values:

```text
parsed_exact
parsed_censored
parsed_approximate
unsupported_measurement_type
unsupported_unit
malformed
missing
```

The binding measurement is metadata about the deposited PDBBind ligand–protein pair. It must not be added to the default DrugCLIP target-fishing score.

Create `index_record_occurrences.parquet` with `complex_id`, `source_index_line_number`, and `source_line_sha256`. It contains no unredacted line text.

---

## 12. Complex discovery and validation

For every index PDB ID, locate exactly one directory under the configured year-range directories.

Expected files:

```text
<pdb_id>_ligand.sdf
<pdb_id>_ligand.mol2
<pdb_id>_pocket.pdb
<pdb_id>_protein.pdb
```

Record:

* missing directory;
* duplicate directory;
* missing file;
* empty file;
* unexpected filename;
* extra files;
* case mismatch;
* unreadable file.

An extra file is not fatal. A missing protein file is fatal. A missing `pocket.pdb` file prevents comparison but does not necessarily prevent re-extraction. At least one usable ligand format is required.

---

## 13. Ligand parsing and validation

### 13.1 Selection policy

Try the ligand formats in this order:

1. `*_ligand.sdf`
2. `*_ligand.mol2`

The SDF file is the preferred source when it is valid. The MOL2 file is the fallback.

Never generate a new ligand conformation. Never perform geometry optimization. Preserve the PDBBind-provided bound coordinates.

### 13.2 RDKit loading

Attempt staged loading:

1. Parse without sanitization.
2. Verify that a molecule object exists.
3. Verify that a conformer with three-dimensional coordinates exists.
4. Attempt RDKit sanitization.
5. Record the sanitization result and exception details.
6. Calculate chemical identifiers only when the corresponding operation succeeds.

An SDF containing more than one nonempty record is distinct from one molecule containing disconnected components. Under `multiple_sdf_record_policy = "reject"`, reject that format with an explicit issue rather than silently choosing the first record; MOL2 fallback may still succeed. Future selection policies must be versioned and retain every record's ordinal and validation result.

A usable ligand requires:

* at least one atom;
* at least one heavy atom;
* finite coordinates shaped `(N, 3)` from an explicit conformer;
* matching atom and coordinate counts.

Do not rely only on RDKit's `Is3D` flag or nonzero Z variance: a legitimate bound ligand can be planar. Record the conformer's dimensional flag and geometric rank as diagnostics. Sanitization failure is a warning if geometry and element identities remain usable. It becomes fatal only if neither format supplies a usable ligand.

### 13.3 Compare SDF and MOL2

When both files are usable, compare:

* atom count;
* heavy-atom count;
* element multiset;
* formal charge;
* connected-component count;
* bond count;
* graph isomorphism when feasible;
* canonical SMILES;
* isomeric SMILES;
* coordinate RMSD after deterministic atom mapping when feasible.

Record disagreement rather than silently choosing one representation.

Default selection:

```text
valid SDF                    → select SDF
invalid SDF, valid MOL2      → select MOL2
both valid and consistent    → select SDF
both valid but inconsistent  → select SDF and emit a warning
both invalid                 → reject geometry
```

Make the selection policy configurable.

Graph comparison must use a deterministic one-to-one mapping. Enumerate symmetry-equivalent isomorphisms up to a configured cap, choose the mapping with minimum direct-frame RMSD, and break ties by lexicographically ordered source atom indices. Do not rigidly align SDF and MOL2 coordinates for the primary RMSD because both should describe the same bound frame; an aligned RMSD may be stored separately as a diagnostic. Report mapping ambiguity and cap exhaustion.

### 13.4 Ligand chemistry fields

Store, when available:

* canonical SMILES;
* isomeric SMILES;
* InChI;
* InChIKey;
* molecular formula;
* formal charge;
* molecular weight;
* heavy-atom count;
* atom-element counts;
* component count;
* RDKit sanitization status;
* stereochemistry status;
* source format;
* source-file SHA-256.

Do not remove salts or disconnected components from stored bound geometry. Record one row per component in `ligand_components.parquet`, including atom indices, composition, centroid, pairwise component separation, and `is_pocket_defining`. Version 1 uses all components for extraction under `pocket_defining_component_policy = "all_components"`; flag spatially separated components and make alternative future policies explicit and derivation-affecting. A standardized parent identifier may be derived separately but must never replace bound geometry.

### 13.5 Ligand content and derivation hashes

Create two hashes:

* `ligand_geometry_content_hash`: atomic numbers, component boundaries, and selected coordinates in a source-order-independent canonical permutation;
* `ligand_derivation_hash`: parser schema version, selected source-file SHA-256, record and conformer ordinal, selection/sanitization/component policies, and the content hash.

For the content hash, order atoms within each component lexicographically by atomic number and exact normalized `(x,y,z)` float64 byte values, then order components by their canonical atom sequence; use source index only to make indistinguishable duplicate rows stable without serializing it. The derivation hash retains original atom order and mapping. Schema version 1 serializes coordinates as C-contiguous little-endian float64 after normalizing negative zero, with length-prefixed UTF-8 strings and arrays, explicit null markers, and fixed field order. NaN and infinity are rejected before hashing. Do not round source coordinates for the content hash. These rules apply to all canonical hashes unless a subsection explicitly defines another representation.

---

## 14. Protein parsing

### 14.1 Source

Use:

```text
<pdb_id>_protein.pdb
```

as the only source of protein geometry for re-extraction.

### 14.2 Model policy

Default:

```text
model_policy = "first"
```

If no `MODEL` records exist, treat the structure as model 1.

Record:

* number of models;
* selected model;
* whether the structure appears to be NMR;
* whether additional models were ignored.

Future support for one pocket per NMR model may be added, but it is outside the first implementation.

### 14.3 Alternate-location policy

Resolve alternate locations deterministically per atom identity defined without the alternate-location field.

Default precedence:

1. highest numeric occupancy, treating missing occupancy as lower than any numeric value;
2. blank alternate-location identifier;
3. alternate location `A`;
4. lexicographically smallest alternate-location identifier.

Record the number of discarded records, missing occupancies, zero-occupancy selections, and exact ties. Blank and nonblank conformers with the same identity participate in the same choice. Retain the chosen `altloc` in the canonical atom key.

### 14.4 Element determination

Use the PDB element field when present.

When absent, infer the element conservatively from the atom name and residue context. Record every inferred element.

Do not confuse atom names such as `CA` for an alpha carbon with calcium.

### 14.5 Protein versus non-protein atoms

Geometry classification must be deterministic from the PDBbind file and local versioned policy only:

* retain non-hydrogen `ATOM` records only when their normalized residue name is in the versioned canonical amino-acid classification table;
* retain non-hydrogen `HETATM` residues only when their normalized residue name is in the configured, versioned modified-amino-acid allowlist;
* exclude water, non-polymer ligands, ions, cofactors, and all other `HETATM` records from the DrugCLIP tensor;
* retain excluded nearby components in sidecars;
* reject or explicitly map any retained atom whose element cannot be represented losslessly by DrugCLIP.

The canonical table must explicitly list the 20 standard amino acids and any accepted PDB residue aliases; nucleic-acid `ATOM` records are excluded. Unknown residues are retained in parsed sidecars but excluded from geometry and flagged. RCSB mapping may confirm, dispute, or annotate this local classification but must not change geometry. A dispute produces a mapping/quality issue. The canonical table, modified-residue allowlist, element policy, and their hashes are derivation inputs.

### 14.6 Canonical atom identity

Represent a PDBBind protein atom using:

```text
pdb_id
model_id
record_type
auth_chain_id
auth_residue_number
insertion_code
residue_name
atom_name
altloc
element
```

Also retain:

```text
serial_number
occupancy
b_factor
x
y
z
source_order
```

Define a canonical `pdbbind_atom_key` from `distribution_id`, `complex_id`, and the stable identity fields. Do not use only atom serial numbers.

### 14.7 Canonical residue identity

Use:

```text
pdb_id
model_id
auth_chain_id
auth_residue_number
insertion_code
residue_name
```

Prefix the residue key with `distribution_id` and `complex_id`. Do not discard insertion codes.

---

## 15. Pocket extraction

### 15.1 Comparable pocket representations

Retain three explicit, uncropped representations:

1. `contact_atom`: locally classified protein heavy atoms individually within the cutoff;
2. `contact_residue`: residues containing at least one `contact_atom`;
3. `residue_expanded_atom`: all locally classified protein heavy atoms belonging to a `contact_residue`, including atoms beyond the cutoff.

A residue belongs to `contact_residue` when at least one of its heavy atoms lies within the configured cutoff of at least one pocket-defining ligand heavy atom.

For residue (r), include it when:

$$
\min_{a \in r,\ \ell \in L}
\lVert x_a-x_\ell\rVert
\leq d
$$

where:

* (L) is the selected ligand heavy-atom set;
* (d) is `distance_cutoff_angstrom`.

### 15.2 DrugCLIP export view

A protein heavy atom belongs to the DrugCLIP pocket when:

$$
\min_{\ell \in L}
\lVert x_a-x_\ell\rVert
\leq d
$$

Version 1 exports the `contact_atom` view. Store this choice as `drugclip_export_view = "contact_atom"` and include it in derivation hashes. The compatibility characterization must compare this view and the `residue_expanded_atom` view against representative DrugCLIP training records; changing the export view requires a new extraction schema version. All three views remain available in sidecars.

### 15.3 Distance computation

Requirements:

* use ligand heavy atoms only by default;
* use protein heavy atoms only for the LMDB view;
* use vectorized NumPy computation or a spatial index;
* avoid constructing an unnecessarily large full distance matrix for extreme structures;
* calculate each selected protein atom’s minimum ligand distance;
* preserve unrounded float64 distances internally;
* store distances in Ångström;
* reject NaN or infinite coordinates.

Use Euclidean Cartesian distance with inclusive `<=` cutoff, no periodic boundary conditions, and Ångström units. A blocked vectorized calculation or deterministic spatial index may be used; tests must prove identical selections at cutoff boundaries. Distance values stored for reproducible hashing must normalize negative zero.

### 15.4 Chain set

The contributing-chain set is the sorted unique set of `auth_chain_id` values among the selected protein atoms.

A valid pocket may involve multiple chains.

Do not automatically collapse:

```text
["A", "B"]
```

into one protein assignment.

### 15.5 Deterministic atom order

Before hashing and LMDB export, order selected atoms deterministically by:

1. minimum distance to a ligand heavy atom;
2. `auth_chain_id`;
3. residue number;
4. insertion code;
5. residue name;
6. atom name;
7. alternate-location identifier;
8. original source order.

Normalize blank chain and insertion identifiers to the empty string for ordering; numeric residue number precedes insertion code; strings use Unicode code-point order after no locale-dependent transformation. Use a stable sort. The exact rules are part of the extraction schema.

### 15.6 Deterministic cropping

When the raw re-extracted atom set contains more than `max_pocket_atoms` atoms:

1. Sort atoms by the deterministic ordering above.
2. Retain the nearest `max_pocket_atoms`.
3. Do not sample randomly.
4. Record the original and retained counts.
5. Record the greatest minimum-ligand distance among retained atoms.
6. Record the smallest minimum-ligand distance among discarded atoms.
7. Preserve the complete uncropped atom and residue metadata in sidecars.
8. Export only the deterministic cropped set to LMDB.

Ties in minimum distance use the remaining deterministic order fields. This prevents the verified loader from performing its sampled center-weighted crop because the exported record is already at or below the active maximum.

### 15.7 Pocket content and derivation hashes

Create:

* `pocket_geometry_content_hash` from the ordered exported element tokens and uncentered coordinates only, using little-endian float32 because those exact values are serialized to LMDB;
* `pocket_derivation_hash` from the extraction schema, distribution ID, PDBbind protein SHA-256, ligand derivation hash, all model/altloc/classification/cutoff/export-view/crop/element-mapping policies, ordered atom identities, float64 source coordinates and minimum ligand distances, and the content hash.

The derivation hash changes whenever any derivation input changes, even if the resulting exported geometry happens to be identical. The content hash changes only when exported tokens or float32 coordinates change. Record both. Use length framing, explicit type tags, little-endian arrays, fixed field order, normalized negative zero, and no locale-dependent serialization.

### 15.8 Pocket instance identifier

Use:

```text
pocket_instance_id =
    "pb20v24p-<distribution_hash_8>:<pdb_id>:<first_16_hex_chars_of_derivation_hash>"
```

Example:

```text
pb20v24p-a1b2c3d4:2l3r:8f3a25b19c173e62
```

Do not include an RCSB ligand-instance mapping in this identifier because enrichment may be unresolved or may improve later. The identifier must depend only on source geometry and versioned extraction rules.

---

## 16. PDBBind-provided pocket comparison

### 16.1 Source

Parse:

```text
<pdb_id>_pocket.pdb
```

using the same model, alternate-location, element, local classification, and hydrogen policies used for extraction. Parse excluded water, ions, and other components for diagnostics, but compare like with like after filtering.

Derive the provided pocket's views explicitly: its `residue_expanded_atom` view is every locally classified protein heavy atom in the file, and its `contact_atom` view is the subset individually within the configured ligand cutoff. Residues are derived from the filtered protein atoms. Thus whole-residue atoms beyond the cutoff do not create a false contact-view disagreement.

Do not crop either pocket before performing the primary comparison.

### 16.2 Comparison levels

Perform separate comparisons for `contact_atom` and `residue_expanded_atom`, identified by `comparison_view`, at:

1. chain level;
2. residue level;
3. atom-identity level;
4. coordinate level;
5. ligand-distance level.

### 16.3 Required summary metrics

Record for each comparison view:

* re-extracted atom count;
* PDBBind-provided atom count;
* re-extracted heavy-atom count;
* PDBBind-provided heavy-atom count;
* common atom-identity count;
* atoms only in re-extracted pocket;
* atoms only in PDBBind pocket;
* atom-set Jaccard similarity;
* re-extracted residue count;
* PDBBind-provided residue count;
* common residue count;
* residues only in re-extracted pocket;
* residues only in PDBBind pocket;
* residue-set Jaccard similarity;
* re-extracted chain set;
* PDBBind-provided chain set;
* chain-set equality;
* common-atom coordinate RMSD;
* common-atom maximum coordinate difference;
* maximum minimum-ligand distance in each pocket;
* mean minimum-ligand distance in each pocket;
* median minimum-ligand distance in each pocket;
* 95th percentile of minimum-ligand distance in each pocket;
* whether one atom set is a subset of the other;
* whether one residue set is a subset of the other.

The primary PDBbind-provided-pocket comparison is `residue_expanded_atom`, because provided pocket files commonly contain whole residues. `contact_atom` is a secondary extraction-boundary diagnostic. Never mix the two views in one Jaccard calculation.

The atom-set Jaccard similarity is:

$$
J(A,B)=\frac{|A\cap B|}{|A\cup B|}
$$

Define `J(empty, empty) = 1` and `J(empty, nonempty) = 0`, while also emitting the relevant empty-pocket issue.

### 16.4 Atom matching

Primary matching must use the canonical atom identity.

When atom identities cannot be matched because PDBBind processing altered identifiers, permit a fallback match using:

* chain;
* residue number;
* insertion code;
* residue name;
* atom name;
* element;
* coordinate tolerance.

Fallback matching must be a deterministic one-to-one assignment. Within each chain/residue/name group, choose the minimum-coordinate-distance assignment under tolerance; break equal-cost ties lexicographically by canonical atom key. Never count one atom more than once. Mark ambiguous equal-cost assignments and cap exhaustion. Coordinate RMSD is calculated directly in the shared PDBbind frame without rigid alignment.

### 16.5 Difference details

Store every atom present in only one representation in a normalized sidecar, not only as counts.

Required classifications:

```text
only_reextracted
only_pdbbind
matched_exact_identity
matched_fallback
ambiguous_unmatched
```

Do not embed potentially large atom-difference lists in `pockets.parquet`.

### 16.6 Comparison warnings

Suggested warning codes:

```text
PDBBIND_POCKET_MISSING
PDBBIND_POCKET_PARSE_FAILED
POCKET_CHAIN_SET_DIFFERS
POCKET_ATOM_JACCARD_LOW
POCKET_RESIDUE_JACCARD_LOW
POCKET_COORDINATE_MISMATCH
PDBBIND_POCKET_HAS_DISTANT_RESIDUE_ATOMS
REEXTRACTED_POCKET_EMPTY
```

Thresholds must be configurable, versioned, and used only by pocket-comparison quality. They must not affect geometry selection, geometry tier, content hashes, or LMDB eligibility.

`PDBBIND_POCKET_HAS_DISTANT_RESIDUE_ATOMS` is informational because whole-residue pocket files are expected to contain such atoms; it never changes comparison quality by itself.

For each primary metric, values at or above the moderate minimum are concordant, values below the severe minimum are severe, and intervening values are moderate. Any severe primary metric makes the view severe; otherwise any moderate metric makes it moderate. Coordinate mismatch uses its own maximum threshold. The primary `residue_expanded_atom` result determines the pocket-level comparison quality, while secondary-view results remain visible.

---

## 17. RCSB mmCIF caching

### 17.1 Download

For each distinct PDB ID, index the canonical request under:

```text
data/cache/external/request_index/rcsb_mmcif/<lowercase_pdb_id>.json
```

The index points to an immutable compressed payload in the content-addressed object store; it is not a second mutable copy of the mmCIF bytes.

Use the RCSB endpoint pattern:

```text
https://files.rcsb.org/download/<PDB_ID>.cif.gz
```

RCSB documents compressed and uncompressed PDBx/mmCIF downloads for PDB entries.

### 17.2 Cache validation

For each downloaded file:

* require a successful HTTP response;
* reject empty content;
* verify gzip readability;
* verify that the mmCIF data block identifies the expected PDB entry;
* calculate SHA-256;
* record ETag and Last-Modified when available;
* use atomic writes;
* never mutate content already referenced by a run.

Cache payloads by content SHA-256 and maintain a replaceable request-to-content index. A run records and thereafter reuses the exact content hash it consumed. `--refresh-cache` may update the request index for future runs but must not alter a prior run's snapshot or allow a resumed stage to drift.

Every HTTP/API cache entry must have a metadata envelope containing canonical request method, URL, normalized parameters or body hash, selected nonsecret request headers, response status, response headers, retrieval time, payload SHA-256, parser/schema version, and error classification. Do not use a single ambiguous `<pdb_id>.json` for several Data API resources. Cache negative responses with bounded expiry, and apply one process-wide rate limiter across workers.

### 17.3 Failure handling

Classify RCSB status as:

```text
current
obsolete
replaced
removed
not_found
download_failed
parse_failed
```

If an entry is obsolete or replaced:

* retain the original PDBBind PDB ID and geometry;
* record the replacement PDB ID when known;
* do not substitute replacement coordinates;
* record replacement metadata only in explicitly marked supplemental fields;
* never use replacement chain, entity, UniProt, ligand, or citation mappings as mappings of the original entry without a separate manual override and evidence trail.

### 17.4 Offline mode

With `--offline`:

* make no network requests;
* use cached files and API responses only;
* mark missing enrichment as unresolved;
* allow geometry processing to continue.

`--offline` and `--refresh-cache` are mutually exclusive. Offline cache misses must not be retried as network failures.

---

## 18. RCSB structural identity and chain mapping

### 18.1 Sources

Use the cached mmCIF as the primary source for:

* `auth_asym_id`;
* `label_asym_id`;
* author and label residue numbering;
* polymer entity IDs;
* polymer versus non-polymer classification;
* model and alternate-location identifiers;
* deposition citations.

Use RCSB Data API or GraphQL as the preferred source for:

* enriched polymer-entity metadata;
* UniProt mappings;
* organism and taxonomy;
* names and descriptions;
* entry-level experimental metadata;
* external database annotations.

The RCSB Data API exposes separate core objects for entries, polymer entities, non-polymer entities, assemblies, and polymer entity instances.

### 18.2 Chain mapping

For every PDBbind contributing `auth_chain_id`:

1. Find matching polymer atom records in mmCIF using `auth_asym_id`.
2. Determine corresponding `label_asym_id` values.
3. Determine the polymer entity ID for each `label_asym_id`.
4. Determine whether the entity is a protein.
5. Retrieve UniProt mappings for the polymer entity.
6. Write every candidate mapping to `chain_mapping_candidates.parquet` with method, evidence, score components, rank, and selected/ambiguous flags.

`protein_chains.parquet` contains one stable PDBbind chain row only. It stores the final overall status but no scalar RCSB identity that would erase multiplicity.

### 18.3 Ambiguity resolution

If a direct `auth_asym_id` mapping is ambiguous:

1. Match PDBBind residues to mmCIF residues using author residue number, insertion code, residue name, and atom name.
2. Compare coordinates under the configured tolerance.
3. Use sequence/order evidence before coordinates when duplicated author chain IDs or symmetry-related copies exist.
4. If necessary, align the PDBbind processed protein to each plausible RCSB asymmetric-unit candidate using shared alpha-carbon atoms and the Kabsch algorithm.
5. Reattempt residue and atom matching in the aligned frame.
6. Require the configured minimum number of alignment atoms.
7. Record the transform, RMSD, coverage, matched-atom count, and competing candidates.
8. Leave the mapping unresolved when evidence remains insufficient.

Selection must use a versioned rule table rather than an undocumented aggregate score. Exact identifier candidates are selected only when unique and residue/entity evidence is consistent. Exact/aligned atom candidates must meet configured minimum matched atoms, residue coverage, and maximum RMSD. A `probable_chain_match` requires an explicit lower threshold and must never be reported as exact. Candidates within the configured tie margin remain `ambiguous`; do not choose one by row order. Record all threshold values and rule versions in configuration and candidate rows.

### 18.4 Coordinate provenance rule

RCSB coordinates may be used only for:

* mapping;
* alignment;
* validation;
* ligand-instance identification;
* structural annotations.

They must not be used to replace coordinates in:

```text
pocket_atoms
pocket_coordinates
pocket_geometry_content_hash
pocket_derivation_hash
```

### 18.5 UniProt mappings

Store one row per:

```text
pocket chain
  × RCSB polymer entity
  × UniProt accession
  × contiguous mapping segment
```

Record:

* accession;
* isoform identifier when supplied;
* mapping source;
* mapped sequence segment;
* mapping coverage;
* mapping status;
* whether multiple UniProt accessions map to the same entity.

Do not collapse isoforms or discontinuous segments. Compute accession-level coverage from the union of its segments and store it in a separate accession summary row/table; never infer coverage by spanning gaps.

### 18.6 Mapping statuses

Use:

```text
exact_identifier_match
exact_atom_match
aligned_atom_match
probable_chain_match
ambiguous
unresolved
rcsb_unavailable
not_a_protein
```

---

## 19. Ligand-instance mapping to RCSB

The PDBBind ligand SDF or MOL2 may not preserve a complete PDB ligand-instance identifier.

Attempt to associate it with an RCSB non-polymer or BIRD instance using:

1. element composition;
2. heavy-atom count;
3. CCD or BIRD identity when inferable;
4. graph isomorphism;
5. author chain and residue information when present;
6. coordinates after any protein-based alignment;
7. coordinate RMSD.

Store every candidate in `rcsb_ligand_mapping_candidates.parquet`, including:

* `rcsb_nonpolymer_entity_id`;
* `label_asym_id`;
* `auth_asym_id`;
* `auth_seq_id`;
* CCD ID;
* BIRD ID where applicable;
* match RMSD;
* match method;
* match status.

Allowed statuses:

```text
exact
probable
ambiguous
unresolved
not_attempted
```

Include a stable candidate ID, evidence components, rank, selected flag, and ambiguity group. Support non-polymer, branched/BIRD, and polymer-peptide candidates rather than assuming every ligand is non-polymer. Do not require successful mapping for geometry acceptance or LMDB entry.

---

## 20. Bibliographic enrichment

### 20.1 Optional post-build scope, sources, and priority

Bibliographic enrichment is an optional post-geometry stage. The default build extracts mmCIF/RCSB citation metadata already present in cache but does not query Crossref, PubMed, or authenticated PDBbind pages. Geometry-to-LMDB acceptance is complete when every measurement has `not_attempted` or another valid final reference status.

Enrich references from:

1. mmCIF `_citation` and `_citation_author`;
2. RCSB Data API or GraphQL;
3. PDBbind/PDBbind+ entry page only when access and automated-use authorization are configured;
4. Crossref using DOI or bibliographic matching;
5. PubMed using PMID, DOI, title, author, journal, year, or ECitMatch;
6. manual overrides.

The mmCIF `_citation.id` value `primary` identifies the citation considered most pertinent by the depositor, but that does not prove that the article contains the affinity measurement.

The `_citation_author` category connects authors to citation records using `citation_id`.

Crossref exposes deposited scholarly metadata through a public REST API.

PubMed E-utilities requests must identify the tool and a contact email and use an API key at higher rates. Crossref requests should supply a contact address. Read credentials/contact values only from the configured environment variables; never persist their values. Store request/response provenance through the content-addressed cache contract.

### 20.2 Legacy PDF token

Do not store:

```text
pdb_bind_reference_token
```

Do not store any PDF-looking token removed from post-`//` index text anywhere. Apply the same scrubber to exception and external-page diagnostic text before logging.

### 20.3 Citation identity

Prefer identifiers in this order:

1. DOI;
2. PMID;
3. normalized bibliographic fingerprint;
4. local source-citation ID scoped to the PDB entry.

A global citation record must not contain a scalar `pdb_id`. Connect global citations to PDB entries through `pdb_citation_links.parquet`. A bibliographic fingerprint may contain:

```text
normalized title
publication year
journal
volume
first page
first author
```

### 20.4 Reference-status model

Every binding measurement must have exactly one final affinity-reference adjudication status:

```text
exact_affinity_reference
probable_affinity_reference
probable_structural_reference
structural_reference_only
conflicting_references
reference_unresolved
no_reference_available
not_attempted
```

Definitions:

* `exact_affinity_reference`: The source has been directly verified to report the relevant measurement or is explicitly identified by a sufficiently authoritative structured source as the measurement source.
* `probable_affinity_reference`: Strong bibliographic and contextual evidence links the source to the measurement, but the measurement itself has not been directly verified.
* `probable_structural_reference`: The source is probably the structural article for the PDB entry, with no evidence that it reports the affinity value.
* `structural_reference_only`: The source is confirmed as a structural citation and is known not to be established as the affinity source.
* `conflicting_references`: Different sources imply incompatible affinity-reference assignments.
* `reference_unresolved`: Citations exist, but none can be assigned confidently.
* `no_reference_available`: No usable citation metadata is available.
* `not_attempted`: Bibliographic enrichment was disabled or not run.

Candidate links and their evidence are many-to-many; the final adjudication is a separate one-row-per-measurement table. An RCSB primary citation alone must not produce `exact_affinity_reference`.

Use a versioned deterministic adjudication precedence:

1. a valid manual override wins and records the overridden automatic result;
2. conflicting valid manual overrides produce `conflicting_references` and a fatal override-validation issue;
3. direct verified measurement evidence produces `exact_affinity_reference`;
4. incompatible high-confidence candidates produce `conflicting_references`;
5. the highest supported probable affinity or structural rule applies;
6. existing but unassignable citations produce `reference_unresolved`;
7. absence of usable citation metadata produces `no_reference_available`;
8. disabled enrichment produces `not_attempted` unless cached/source-only evidence already establishes a stronger status.

Define every automatic evidence rule, confidence contribution, tie threshold, and source-priority value in a versioned rule table. Confidence is float64 in `[0,1]` and is never interpreted without the rule version.

### 20.5 Manual overrides

Support:

```text
data/config/affinity_reference_overrides.parquet
```

or CSV with:

```text
complex_id
measurement_id
citation_id
reference_status
evidence_note
verified_by
verified_at_utc
```

Manual overrides must be auditable and must not modify source-derived citation records.

---

## 21. Sidecar data model

Parquet sidecars are the authoritative metadata store. `schemas.py` must define, for every table, exact Arrow types, nullability, primary key, foreign keys, allowed enums, canonical sort order, and any volatile audit columns. Writers must instantiate these schemas even for empty tables and reject implicit type inference. Use UTF-8 strings, signed fixed-width integers, float64 analytical values, UTC timestamps, and `list<string>` where specified. JSON columns are allowed only for genuinely open-ended diagnostics; normalized scientific relationships require tables. Every table must include schema name and semantic version in Parquet metadata.

Volatile audit values such as filesystem modification time, download time, issue creation time, and run timestamps must never enter identifiers, derivation hashes, scientific sort keys, or logical table digests. Each table's logical digest is computed over its schema-defined nonvolatile columns and canonical row order; reports must identify excluded columns. Physical Parquet checksums still cover every byte.

### 21.1 `complexes.parquet`

One row per unique PDBbind complex. Physical duplicate occurrences belong in `index_record_occurrences.parquet`.

Required columns:

```text
complex_id
pdb_id
distribution_id
nominal_complex_set_version
structure_processing_version
index_revision_date
primary_index_line_number
index_line_redacted
source_line_sha256
release_year
resolution_raw
resolution_angstrom
experimental_method_hint
ligand_label
index_comment
complex_directory
protein_file_id
ligand_sdf_file_id
ligand_mol2_file_id
pdbbind_pocket_file_id
rcsb_entry_status
processing_status
geometry_quality_tier
pocket_comparison_quality
structure_mapping_quality
bibliography_quality
warning_count
error_count
```

### 21.2 `ligand_instances.parquet`

One row per PDBBind-selected ligand geometry.

```text
ligand_instance_id
complex_id
pdb_id
selected_source_format
selected_source_file_id
ligand_geometry_content_hash
ligand_derivation_hash
rdkit_parse_status
rdkit_sanitization_status
canonical_smiles
isomeric_smiles
inchi
inchikey
molecular_formula
formal_charge
molecular_weight
atom_count
heavy_atom_count
component_count
element_counts_json
stereochemistry_status
sdf_mol2_comparison_status
sdf_mol2_coordinate_rmsd
rcsb_ligand_match_overall_status
warnings
```

`ligand_components.parquet` contains one row per selected ligand component with `ligand_instance_id`, `component_index`, atom-index list, atom/heavy-atom counts, element counts, formal charge, centroid coordinates, minimum separation to another component, and `is_pocket_defining`.

### 21.3 `pockets.parquet`

One row per re-extracted pocket.

```text
pocket_instance_id
complex_id
ligand_instance_id
pdb_id
pocket_geometry_content_hash
pocket_derivation_hash
extraction_schema_version
distance_cutoff_angstrom
selected_model_id
model_count
altloc_policy
hydrogen_policy
contact_atom_count
residue_expanded_atom_count
exported_atom_count
contact_residue_count
drugclip_export_view
contributing_chain_count
contributing_auth_chain_ids
minimum_ligand_distance_min
minimum_ligand_distance_mean
minimum_ligand_distance_median
minimum_ligand_distance_max
crop_applied
crop_max_atoms
maximum_retained_ligand_distance
minimum_discarded_ligand_distance
all_elements_supported
processing_status
geometry_quality_tier
pocket_comparison_quality
structure_mapping_quality
bibliography_quality
warning_codes
error_codes
lmdb_profile_memberships
```

### 21.4 `pocket_residues.parquet`

One row per `contact_residue` residue. Scalar RCSB fields are populated only for a unique selected mapping; all candidates remain in mapping-candidate tables.

```text
pocket_instance_id
pdb_id
model_id
auth_chain_id
auth_residue_number
insertion_code
residue_name
minimum_ligand_distance
selected_atom_count
total_heavy_atom_count
rcsb_mapping_status
rcsb_label_asym_id
rcsb_label_seq_id
rcsb_polymer_entity_id
```

### 21.5 `pocket_atoms.parquet`

One row per locally classified protein heavy atom in the `residue_expanded_atom` view.

```text
pocket_instance_id
pdbbind_atom_key
source_order
model_id
record_type
auth_chain_id
auth_residue_number
insertion_code
residue_name
atom_name
altloc
element
occupancy
b_factor
x
y
z
minimum_ligand_distance
in_contact_atom_view
in_residue_expanded_atom_view
retained_after_crop
export_order
element_supported_by_drugclip
rcsb_atom_mapping_status
rcsb_label_asym_id
rcsb_label_seq_id
rcsb_atom_id
rcsb_polymer_entity_id
```

This table is potentially large but remains manageable for approximately tens of thousands of pockets.

### 21.6 `pocket_comparisons.parquet`

One row per re-extracted/PDBbind-pocket comparison view.

```text
pocket_instance_id
comparison_view
pdbbind_pocket_file_id
comparison_status
reextracted_atom_count
pdbbind_atom_count
reextracted_heavy_atom_count
pdbbind_heavy_atom_count
common_atom_exact_count
common_atom_fallback_count
only_reextracted_atom_count
only_pdbbind_atom_count
atom_jaccard
reextracted_residue_count
pdbbind_residue_count
common_residue_count
only_reextracted_residue_count
only_pdbbind_residue_count
residue_jaccard
reextracted_chain_ids
pdbbind_chain_ids
chain_sets_equal
reextracted_subset_of_pdbbind
pdbbind_subset_of_reextracted
common_atom_coordinate_rmsd
common_atom_max_coordinate_difference
reextracted_maximum_ligand_distance
pdbbind_maximum_ligand_distance
reextracted_mean_ligand_distance
pdbbind_mean_ligand_distance
reextracted_median_ligand_distance
pdbbind_median_ligand_distance
reextracted_p95_ligand_distance
pdbbind_p95_ligand_distance
warning_codes
```

### 21.7 `pocket_atom_differences.parquet`

One row per compared atom.

```text
pocket_instance_id
comparison_view
comparison_class
reextracted_atom_key
pdbbind_atom_key
match_method
coordinate_distance
auth_chain_id
auth_residue_number
insertion_code
residue_name
atom_name
element
```

### 21.8 `protein_chains.parquet`

One row per contributing PDBBind chain.

```text
pocket_instance_id
pdb_id
pdbbind_auth_chain_id
selected_atom_count
selected_residue_count
rcsb_mapping_status
warnings
```

### 21.9 `chain_mapping_candidates.parquet`

One row per PDBbind-chain/RCSB-instance candidate.

```text
chain_mapping_candidate_id
pocket_instance_id
pdb_id
pdbbind_auth_chain_id
rcsb_label_asym_id
rcsb_auth_asym_id
rcsb_polymer_entity_id
entity_type
entity_description
organism_name
taxonomy_id
sequence_length
mapping_method
mapping_status
candidate_rank
selected
identifier_match_count
atom_match_count
residue_coverage
alignment_atom_count
alignment_rmsd
transform_json
evidence_codes
warning_codes
```

### 21.10 `chain_uniprot_mappings.parquet`

One row per chain/entity/UniProt accession summary. Discontinuous details belong in the segment table.

```text
pocket_instance_id
pdb_id
pdbbind_auth_chain_id
rcsb_label_asym_id
rcsb_polymer_entity_id
uniprot_accession
uniprot_isoform
mapping_source
mapping_status
mapping_coverage
```

`chain_uniprot_mapping_segments.parquet` contains one row per contiguous segment with a stable segment ID and PDB/UniProt begin and end positions.

### 21.11 `rcsb_ligand_mapping_candidates.parquet`

One row per possible RCSB non-polymer, branched/BIRD, or polymer-peptide ligand association, with candidate ID, ligand instance ID, entity/instance identifiers, CCD/BIRD identifiers, match method, composition and graph evidence, RMSD, rank, selected flag, ambiguity group, and status.

### 21.12 `citations.parquet`

One row per distinct citation.

```text
citation_id
title
journal
year
volume
issue
first_page
last_page
doi
pmid
crossref_id
publication_status
source_priority
bibliographic_fingerprint
metadata_sources
conflict_status
```

`pdb_citation_links.parquet` connects `citation_id` to `pdb_id`, source-specific citation ID, role, source, source priority, and evidence. It preserves multiple PDB associations without duplicating a global DOI/PMID citation.

### 21.13 `citation_authors.parquet`

One row per citation author.

```text
citation_id
ordinal
author_name
orcid
source
```

Author rows from conflicting sources must remain source-qualified. Their key is `(citation_id, source, ordinal)`; do not merge ordinals silently.

### 21.14 `affinity_reference_links.parquet`

One row per binding-measurement/reference candidate.

```text
measurement_id
complex_id
citation_id
candidate_status
confidence
evidence_sources
evidence_note
automatic_or_manual
verified_by
verified_at_utc
```

`affinity_reference_adjudications.parquet` has exactly one row per `measurement_id`, containing final `reference_status`, selected citation ID when any, rule/manual-override version, confidence, evidence summary, adjudicator fields, and timestamp. Candidate links never substitute for this final row.

### 21.15 `nearby_nonprotein_components.parquet`

Record excluded but nearby components such as metals, cofactors, waters, and non-polymer residues.

“Nearby” means at least one component atom lies within `nearby_nonprotein_cutoff_angstrom` of a pocket-defining ligand heavy atom. Use one row per residue/component instance, not one row per element; element composition belongs in `element_counts_json`.

```text
pocket_instance_id
component_type
component_id
auth_chain_id
auth_seq_id
residue_name
insertion_code
atom_count
element_counts_json
minimum_ligand_distance
included_in_drugclip_tensor
exclusion_reason
```

### 21.16 `processing_issues.parquet`

One row per warning or error.

```text
issue_id
stage
complex_id
pocket_instance_id
severity
issue_code
message
exception_type
source_file_id
created_at_utc
details_json
```

`issue_id` must be deterministic from stage, affected identifiers, issue code, and a canonical details fingerprint; timestamps are not ID inputs. Exception text and details pass through the PDF-token and secret scrubbers.

Severity:

```text
info
warning
error
fatal
```

### 21.17 `lmdb_records.parquet`

One row per exported LMDB record.

```text
library_profile
lmdb_path
record_index
lmdb_key
pocket_instance_id
pocket_geometry_content_hash
pocket_derivation_hash
atom_count
serialized_record_sha256
logical_record_sha256
```

---

## 22. Quality-control system

### 22.1 Geometry processing status

Use:

```text
accepted
accepted_with_warnings
rejected
not_processed
```

### 22.2 Independent quality dimensions

Store four independent dimensions:

```text
geometry_quality_tier = A | B | C | rejected | not_processed
pocket_comparison_quality = concordant | moderate_difference | severe_difference | unavailable | not_processed
structure_mapping_quality = exact | aligned | ambiguous | unresolved | unavailable | not_processed
bibliography_quality = exact | probable | unresolved | unavailable | not_attempted
```

Only `geometry_quality_tier` controls the default LMDB profile. RCSB, UniProt, ligand-instance, pocket-comparison, or bibliography outcomes never downgrade geometry. Multiple contributing chains are informational and never change any tier by themselves.

Implement a versioned `quality_rules.toml` or equivalent code table mapping every issue code to quality dimension, severity, and tier effect. Unknown issue codes are validation errors. The resolved rules and SHA-256 belong in the run manifest.

For geometry schema version 1:

* Tier A is accepted geometry with none of the Tier B or Tier C geometry codes below.
* Tier B is accepted geometry with `DETERMINISTIC_CROP_APPLIED`, `LIGAND_MOL2_FALLBACK`, `SDF_MOL2_DISAGREEMENT`, `MODIFIED_RESIDUE_INCLUDED`, or another rule explicitly mapped to B.
* Tier C is accepted but excluded from the default profile because at least one explicitly computed C rule applies.

Version 1 Tier C rules are computable as follows:

| Issue code | Rule |
| --- | --- |
| `PROBABLE_COVALENT_CONTACT` | Minimum ligand/protein heavy-atom distance is at most the sum of configured covalent radii plus `0.40 Å`; label as probable, not confirmed. |
| `EXCLUDED_COMPONENT_BRIDGES_CONTACT` | A non-water excluded heavy atom is within `3.0 Å` of both a pocket-defining ligand heavy atom and a contact protein heavy atom. |
| `VERY_SMALL_POCKET` | Exported atom count is below `minimum_pocket_atoms_warning`. |
| `UNSUPPORTED_ATOM_EXCLUDED` | A locally classified pocket atom was removed by an explicit element policy. |
| `SPATIALLY_SEPARATED_LIGAND_COMPONENTS` | Two pocket-defining ligand components have minimum intercomponent heavy-atom distance above a configured threshold, default `8.0 Å`. |
| `LOCAL_MISSING_ATOM_RECORD` | The processed PDB file explicitly reports a missing atom or residue affecting a contact residue; absence inferred only from RCSB does not trigger this geometry rule. |
| `LIGAND_CHEMISTRY_UNUSUAL` | A configured subrule fires, such as unresolved element, impossible formal charge parsing, or failed valence perception; each subrule must be named in details. |

Do not infer “strong dependence,” “substantial missing atoms,” covalency, or unusual chemistry without one of these explicit rules. Tier C pockets remain in sidecars and require an explicit LMDB profile.

### 22.3 Hard rejection conditions

Reject a pocket from LMDB export when:

* both ligand formats are unusable;
* ligand coordinates are absent or nonfinite;
* protein parsing fails;
* no pocket protein atoms are selected;
* atom and coordinate array lengths differ;
* exported coordinates are not shaped (N\times3);
* an exported element is unsupported and no explicit policy handles it;
* deterministic cropping fails;
* either required geometry/derivation hash cannot be generated;
* the pocket identifier is duplicated with a different geometry;
* a fatal source-integrity issue exists.

Bibliographic, comparison, RCSB, or UniProt failure is never by itself a geometry rejection condition.

---

## 23. DrugCLIP dictionary compatibility

Read the active pocket dictionary from:

```text
data/DrugCLIP/data/dict_pkt.txt
```

or a configured equivalent.

Before export:

1. Remove hydrogen atoms.
2. Apply the verified `AffinityPocketDataset.pocket_atom` transformation to a copy of every proposed input token and assert that it is unchanged.
3. Verify that every post-transformation token exists in the dictionary.
4. Record unsupported or lossy elements in sidecars.
5. Do not silently map elements such as selenium to sulfur.
6. Permit explicit versioned element mappings only through configuration and assign the appropriate geometry issue/tier.
7. Include the dictionary, loader, task, helper, and library-contract SHA-256 values in the manifest and LMDB profile metadata. Do not hash or record a checkpoint as though it affected LMDB serialization.

---

## 24. Final LMDB export

### 24.1 Source of truth

LMDB files must be generated only from completed sidecars. Do not generate LMDB records directly during raw-file parsing.

This permits:

* re-exporting different library profiles;
* changing quality-tier filters without repeating extraction;
* validating every LMDB record against sidecars;
* auditing inclusion and exclusion decisions.

### 24.2 Required default LMDB

Create:

```text
lmdb/candidate_pockets.lmdb
```

This is a single-file LMDB opened with:

```python
lmdb.open(path, subdir=False, ...)
```

Records must be stored under dense numeric ASCII keys:

```text
b"0"
b"1"
b"2"
...
```

Write to a sibling temporary LMDB with `lock=False`, close and validate it, then atomically replace the destination according to `--overwrite`. No lock file is a final artifact. Estimate `map_size` from record sizes plus configured headroom and page alignment; a numeric override remains available. Do not allocate a fixed 1 TiB map by default.

### 24.3 Record ordering

Order records by:

```text
pocket_instance_id ascending
```

The ordering must be independent of:

* filesystem directory order;
* worker count;
* completion order;
* operating system;
* network response order.

### 24.4 Required record schema

Each value is a pickled Python dictionary:

```python
{
    "pocket": str,
    "pocket_atoms": list[str],
    "pocket_coordinates": np.ndarray,
}
```

Constraints:

```text
record["pocket"] == pocket_instance_id
record["pocket_coordinates"].dtype == np.float32
record["pocket_coordinates"].dtype.str == "<f4"
record["pocket_coordinates"].flags.c_contiguous
record["pocket_coordinates"].shape == (N, 3)
len(record["pocket_atoms"]) == N
1 <= N <= max_pocket_atoms
no hydrogen atoms
all coordinates finite
all tokens supported by dict_pkt.txt
pickle protocol == 4
```

Do not include binding measurements, UniProt accessions, citations, quality annotations, residue lists, or other large metadata in the DrugCLIP LMDB record.

The value in `pocket` must be the unique `pocket_instance_id`, not merely the four-character PDB ID. DrugCLIP only requires a pocket-name field and pocket atom/coordinate arrays.

### 24.5 Optional LMDB profiles

Support deterministic filtered exports such as:

```text
lmdb/candidate_pockets_tier_a.lmdb
lmdb/candidate_pockets_tiers_ab.lmdb
lmdb/candidate_pockets_all_usable.lmdb
```

Every profile must have:

* an explicit filter expression;
* a profile name;
* a record count;
* an LMDB checksum;
* an accompanying entry in `lmdb_records.parquet`;
* an entry in the manifest.

The required default profile is geometry Tier A plus Tier B unless configured otherwise. Filters must refer to named, versioned sidecar columns and be parsed by a restricted expression grammar; do not execute arbitrary Python expressions.

Every profile has adjacent `*.profile.json` metadata containing its filter AST, schema versions, record count, physical and logical checksums, the versioned library-contract and dictionary/task/loader/helper hashes, serialization protocol, and source sidecar hashes. Encoder checkpoint identity belongs in a separate embedding or encoder-integration artifact keyed to this profile's logical checksum.

BioSensIA-DC pocket embedding caches are not keyed by LMDB content. Overwriting any candidate LMDB must invalidate or namespace the corresponding embedding cache using the LMDB logical checksum. The build/report output must state the required cache action.

### 24.6 Optional lookup LMDB

An optional lookup LMDB may map:

```text
pocket_instance_id → numeric LMDB key
```

It must be stored separately:

```text
lmdb/candidate_pockets_lookup.lmdb
```

The DrugCLIP-compatible LMDB must remain a dense numeric-key database.

---

## 25. Retrieval-result enrichment contract

A target-fishing result will contain:

```text
pocket_instance_id
DrugCLIP similarity score
rank
```

The application must join `pocket_instance_id` to sidecars to display:

* PDB ID;
* bound reference ligand;
* binding measurement;
* explicit statement about which ligand the measurement applies to;
* contributing chains;
* protein names;
* UniProt accessions;
* organism;
* all four quality dimensions;
* PDBBind-pocket comparison summary;
* citation and affinity-reference status;
* warnings.

The binding measurement must not be presented as a predicted affinity for the target-fishing query unless the standardized query molecule matches the deposited PDBBind ligand.

---

## 26. CLI specification

Provide:

```bash
python -m biosensia_pocket_library.cli build \
    --config config/pdbbind-pocket-library.toml
```

Required subcommands:

```text
check-drugclip-contract
inventory
parse-index
download-rcsb
build-sidecars
export-lmdb
validate
report
build
```

Examples:

```bash
python -m biosensia_pocket_library.cli inventory \
    --index-dir data/raw/index \
    --complex-root data/raw/P-L

python -m biosensia_pocket_library.cli download-rcsb \
    --config config/pdbbind-pocket-library.toml \
    --workers 8

python -m biosensia_pocket_library.cli build-sidecars \
    --config config/pdbbind-pocket-library.toml \
    --resume

python -m biosensia_pocket_library.cli export-lmdb \
    --run-dir data/processed/pdbbind_2020_v2024p_20250804/<run_id> \
    --profile tiers-ab

python -m biosensia_pocket_library.cli validate \
    --run-dir data/processed/pdbbind_2020_v2024p_20250804/<run_id>
```

Useful filters for development:

```text
--pdb-id 2l3r
--pdb-ids-file selected_pdb_ids.txt
--limit 100
--year-from
--year-to
--offline
--refresh-cache
--fail-fast
--resume
--overwrite-run
```

`--limit` must be applied after deterministic sorting by PDB ID so development runs are reproducible.

All PDB/year/limit filters form a canonical selection specification whose hash participates in `run_id`; inventory counts must distinguish the complete source from the selected run subset. `--resume` requires an exact compatible run identity. `--overwrite-run` may replace only the explicitly resolved run directory after validating that it is beneath `output_root`; it must never target `output_root` itself. Prefer moving the previous run to a timestamped backup unless `--discard-existing-run` is separately supplied.

Define exit codes: `0` for a completed command even when individual complexes are rejected as recorded data; `1` for validation failure or an incomplete required acceptance threshold; `2` for configuration/usage errors; and `3` for fatal infrastructure or source-integrity failure. `--fail-fast` stops on the first per-complex error and exits nonzero. Reject incompatible combinations such as `--offline --refresh-cache`.

---

## 27. Logging

Use structured logging.

Create:

```text
logs/pipeline.log
logs/events.jsonl
```

Each JSONL event should include:

```text
timestamp_utc
level
stage
complex_id
pocket_instance_id
event_code
message
worker_id
details
```

Do not include entire molecular structures or large coordinate arrays in logs.

Pass all messages/details through centralized secret and PDF-token scrubbers before either log sink. Unit tests must show that parser exceptions containing an unredacted source line cannot leak a token. Do not log environment secret values, complete external API payloads, molecular structures, or coordinate arrays.

---

## 28. Parallelism

Parallelize complex-level CPU work, but maintain deterministic output.

Requirements:

* workers must not write directly to the same Parquet file;
* workers return structured records to a coordinator or write isolated partitions;
* final tables are sorted before writing;
* LMDB export is single-writer;
* network requests use bounded concurrency and retries;
* cached responses are written atomically;
* a failed worker must produce a processing issue rather than silently dropping a complex.

The final output must be logically reproducible. Canonical table contents, row ordering, content/derivation hashes, LMDB record bytes, and LMDB logical digests must be deterministic across worker counts and completion order. Physical Parquet or LMDB byte identity is desirable but not required across library/filesystem versions; timestamps, compression metadata, and physical page layout are excluded from logical equality. Record all library versions needed to explain physical differences.

---

## 29. Validation

### 29.1 Input validation

Check:

* all expected index files exist;
* the declared PL count in the index header is parsed;
* parsed data-line count matches the declared count;
* parsed unique-complex count and equivalent-duplicate occurrence count are reported separately;
* each PL index record has no more than one discovered complex directory;
* source-file checksums are recorded.

### 29.2 Sidecar validation

Check:

* primary-key uniqueness;
* foreign-key integrity;
* no duplicate `pocket_instance_id`;
* duplicate content hashes are allowed and reported, but never with inconsistent atom tokens/coordinates;
* every derivation hash maps to exactly one derivation description and content hash;
* every accepted pocket has a ligand and complex;
* every exported pocket has atom rows;
* every chain mapping points to a known pocket;
* every affinity-reference link points to known measurement and citation rows;
* every measurement has exactly one final affinity-reference adjudication;
* mapping candidates preserve ambiguity and selected flags are unique only when status permits;
* UniProt mapping segments do not overlap inconsistently and accession coverage is regenerated from their union;
* no forbidden legacy PDF field exists;
* no column name or persisted redacted string matches the configured PDF-token detector;
* schemas, nullability, enums, primary keys, foreign keys, and canonical sort order exactly match `schemas.py`.

### 29.3 LMDB validation

For every record:

```python
assert set(record) == {
    "pocket",
    "pocket_atoms",
    "pocket_coordinates",
}
```

Unless the current BioSensIA-DC retrieval loader explicitly requires another field, no extra fields should be added.

Check:

```python
assert isinstance(record["pocket"], str)
assert isinstance(record["pocket_atoms"], list)
assert isinstance(record["pocket_coordinates"], np.ndarray)
assert record["pocket_coordinates"].dtype == np.float32
assert record["pocket_coordinates"].dtype.str == "<f4"
assert record["pocket_coordinates"].flags.c_contiguous
assert record["pocket_coordinates"].ndim == 2
assert record["pocket_coordinates"].shape[1] == 3
assert len(record["pocket_atoms"]) == len(record["pocket_coordinates"])
assert np.isfinite(record["pocket_coordinates"]).all()
```

Verify:

* keys are dense numeric ASCII strings;
* key count equals sidecar profile count;
* record order agrees with `lmdb_records.parquet`;
* each content hash, derivation hash, serialized-record hash, and logical profile digest can be regenerated;
* each token is unchanged by the verified loader atom-token transformation and exists in the configured dictionary;
* no record exceeds `max_pocket_atoms`.

### 29.4 DrugCLIP integration test

Load the produced LMDB through the actual BioSensIA-DC target-fishing pocket loader.

The complete integration coverage has a loader phase (items 1 and 3–7) and a separately runnable encoder phase (items 2 and 8). The loader phase must not require or hash a checkpoint. The encoder phase must hash the checkpoint because it actually uses it. Together the phases must:

1. Load at least one record.
2. Encode at least one pocket with the configured checkpoint.
3. Verify that no dictionary, shape, dtype, hydrogen, or cropping error occurs.
4. Compare the record immediately before and after the loader's hydrogen-removal/cropping wrappers and verify identical atom count, token order, and coordinates.
5. Verify that loader normalization centers coordinates while preserving all pairwise distances.
6. Verify that the pocket name returned by retrieval is the complete `pocket_instance_id`.
7. Join that ID to `pockets.parquet`.
8. Record the linked BioSensIA-DC commit, library-contract file hashes, exact checkpoint hash, and input LMDB logical checksum in the encoder-test artifact.

---

## 30. Reports

Generate:

```text
reports/build_summary.json
reports/build_summary.md
reports/quality_counts.parquet
reports/failure_counts.parquet
reports/pocket_size_distribution.parquet
reports/pocket_comparison_distribution.parquet
reports/chain_mapping_status_counts.parquet
reports/uniprot_mapping_status_counts.parquet
reports/reference_status_counts.parquet
```

The Markdown summary should include:

* declared and parsed complex counts;
* discovered directory count;
* number successfully processed;
* rejected count by reason;
* geometry Tier A, B, C, rejected, and not-processed counts;
* comparison, mapping, and bibliography quality counts separately;
* ligand SDF success count;
* ligand MOL2 fallback count;
* both-formats-failed count;
* cropped-pocket count;
* multi-chain-pocket count;
* unresolved chain mapping count;
* unresolved UniProt mapping count;
* exact/probable/unresolved reference counts;
* PDBBind-pocket comparison distributions;
* final LMDB profile counts;
* checksums and run ID.

---

## 31. Tests

### 31.1 Unit tests

Implement tests for:

* index-header parsing;
* data-line parsing;
* removal and non-persistence of multiple/path-like/mixed-case PDF tokens;
* binding-measurement parsing;
* unit normalization;
* inequality reversal on the (pX) scale;
* inverse-molar `Ka` normalization without relation reversal;
* ligand SDF selection;
* MOL2 fallback;
* SDF/MOL2 disagreement;
* multiple SDF records versus disconnected components;
* deterministic symmetry-aware ligand atom mapping;
* alternate-location selection;
* PDB element inference;
* heavy-atom filtering;
* distance-based pocket extraction;
* all three pocket representations and their membership relationships;
* multi-chain pockets;
* deterministic cropping;
* deterministic atom ordering;
* independent content and derivation hashing, length framing, endian handling, and negative zero;
* pocket comparison in both views, empty-set Jaccard, and one-to-one fallback matching;
* chain identifier mapping;
* ambiguous chain and ligand mapping candidates;
* discontinuous UniProt mapping segments and coverage union;
* RCSB cache validation;
* citation-author joins;
* reference-status rules;
* sidecar schema validation;
* LMDB record serialization;
* verified loader token transformation, hydrogen no-op, crop no-op, and centroid normalization;
* every quality-rule trigger and proof that enrichment/comparison cannot change geometry tier.

### 31.2 Synthetic fixtures

Create small synthetic PDB, SDF, MOL2, and mmCIF fixtures that cover:

* one chain and one ligand;
* two chains contributing to one pocket;
* alternate locations;
* insertion codes;
* modified residues;
* nearby metal;
* multiple ligand components;
* multiple SDF records;
* spatially separated ligand components;
* more than 256 pocket atoms;
* unresolved RCSB mapping;
* NMR-style multiple models;
* equivalent and conflicting duplicate index records;
* an empty comparison view and ambiguous coordinate fallback;

Do not commit actual PDBBind structures unless the license clearly permits redistribution.

### 31.3 Local-data integration tests

Provide tests marked:

```python
@pytest.mark.pdbbind_data
```

They should skip automatically when the local PDBBind files are unavailable.

Provide `pytest_addoption` for `--pdb-id`; the marked test must use the option when present and otherwise choose a documented deterministic fixture ID.

Include a development test for a user-specified PDB ID such as:

```bash
pytest -m pdbbind_data --pdb-id 2l3r
```

### 31.4 Regression test

For a fixed local input and configuration:

* store expected counts and hashes in a local, noncommitted regression artifact;
* verify that rerunning produces identical pocket IDs, content hashes, derivation hashes, canonical table contents, serialized LMDB records, and logical LMDB digest;
* verify that changing the distance cutoff changes the semantic configuration and every derivation hash, while content hashes change only for pockets whose exported tokens or float32 coordinates change;
* verify identical logical outputs with worker counts 1 and greater than 1;
* verify that refreshed enrichment data cannot change pocket identifiers, geometry-content hashes, or derivation hashes.

---

## 32. Acceptance criteria

The implementation is complete when:

1. The linked BioSensIA-DC library contract, dictionary, and helper pass compatibility validation and their relevant hashes are recorded; checkpoint validation is required only for encoder integration or embedding artifacts.
2. A development run for one PDB ID completes from raw files to LMDB.
3. A full run parses every physical record and reports every unique complex in `INDEX_general_PL.2020R1.lst`.
4. Every artifact inside the defined checksum boundary has a checksum, and every LMDB has a logical digest.
5. The distance cutoff and every geometry-affecting policy are configurable and derivation-hashed.
6. Pocket extraction, ordering, cropping, hashing, and logical outputs are deterministic.
7. Ligand SDF is preferred, MOL2 fallback works, and multiple records/components are handled explicitly.
8. Contact-atom, contact-residue, and residue-expanded views are retained and compared correctly with PDBbind-provided pockets.
9. Contributing protein chains are identified from PDBbind geometry alone.
10. RCSB responses are content-addressed and snapshotted for the run.
11. Chain, entity, ligand, and UniProt mappings preserve candidates, multiplicity, ambiguity, isoforms, and discontinuous segments.
12. No persisted field or redacted text contains a removed PDF-looking token or secret.
13. Binding measurements, including dimensionally distinct `Ka`, are parsed, normalized, and retained.
14. Every measurement has exactly one final affinity-reference status; `not_attempted` is valid when optional enrichment is disabled.
15. RCSB availability or content never changes PDBbind-derived atoms, identifiers, geometry-content hashes, or derivation hashes.
16. Every quality outcome follows the versioned computable rule table and independent quality dimensions.
17. Parquet sidecars pass exact schema, primary-key, foreign-key, enum, and sort-order validation.
18. LMDB records use dense numeric ASCII keys, protocol-4 pickled dictionaries, and lossless supported tokens.
19. Every LMDB record contains a stable unique pocket identifier and regenerable content/derivation hashes.
20. The produced LMDB passes the actual BioSensIA-DC loader and encoder integration test without loader-side removal or cropping.
21. Final manifest and summary reports are generated, embedding-cache invalidation is stated, and no restricted artifacts are Git-tracked.

---

## 33. Suggested implementation sequence

Implement in the following pull-request-sized stages.

### Stage 0: DrugCLIP contract and schema foundations

Implement:

* linked-checkout discovery and revision capture;
* task/loader/helper/dictionary characterization and versioned library-contract identity;
* separate checkpoint characterization only in encoder integration;
* exact Arrow schema registry;
* canonical serialization and hash primitives;
* versioned quality-rule registry.

### Stage 1: Inventory and index parsing

Implement:

* configuration;
* source inventory;
* manifest;
* PL index parser;
* binding-measurement parser;
* basic `complexes.parquet`;
* `binding_measurements.parquet`.

### Stage 2: Ligand and protein parsing

Implement:

* RDKit SDF/MOL2 policy;
* protein PDB parsing;
* model and alternate-location policies;
* ligand and source sidecars.

### Stage 3: Pocket extraction and comparison

Implement:

* configurable distance extraction;
* contact-atom, contact-residue, and residue-expanded views;
* versioned DrugCLIP export view;
* deterministic cropping;
* content and derivation hashes;
* PDBBind-pocket comparison;
* pocket atom and residue sidecars.

### Stage 4: RCSB enrichment

Implement:

* mmCIF download/cache;
* content-addressed RCSB/API cache and run snapshots;
* auth-to-label chain mapping;
* polymer entity mapping;
* UniProt mapping;
* ligand-instance mapping;
* coordinate-alignment fallback.

### Stage 5: Bibliographic enrichment

Implement:

* mmCIF citation extraction;
* citation-author extraction;
* optional Crossref/PubMed enrichment;
* optional authorized PDBbind-page enrichment;
* affinity-reference statuses;
* manual overrides.

### Stage 6: Quality and LMDB export

Implement:

* independent quality dimensions and computable tiers;
* profile filters;
* DrugCLIP dictionary validation;
* required LMDB schema;
* lookup and record mapping;
* integration test.

### Stage 7: Full-run validation and reporting

Implement:

* complete validation suite;
* build reports;
* resume behavior;
* performance improvements;
* full special-distribution run.
