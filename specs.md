# Specification: Build a BioSensIA-DC Candidate Pocket Library from PDBBind 2020R1

## 1. Objective

Implement a reproducible pipeline that builds a candidate protein-pocket library from the protein–ligand complexes distributed with PDBBind 2020R1.

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
PDBBind release
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

The PDBBind processed files are the source of the geometry supplied to DrugCLIP. RCSB data is used to provide structural identity, entity relationships, sequence mappings, UniProt annotations, experimental metadata, and bibliographic enrichment. RCSB coordinates must not silently replace PDBBind coordinates.

---

## 2. Scope

### 2.1 Included input set

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
* use RCSB coordinates instead of the PDBBind processed geometry;
* assume that the primary PDB citation is the source of the PDBBind affinity measurement;
* store the PDBBind legacy PDF filename token;
* commit PDBBind data to Git because its distribution may be license-restricted.

---

## 3. Grounding constraints

The DrugCLIP training-data schema uses the keys `pocket`, `pocket_atoms`, and `pocket_coordinates` for pocket data.

The original DrugCLIP task removes pocket hydrogens, normalizes pocket coordinates, and applies a maximum-pocket-atom setting whose default is 256.

The original `CroppingPocketDataset` performs a seeded but sampled crop when a pocket contains more than the maximum number of atoms. The new library builder must pre-crop deterministically so a candidate’s geometry does not depend on retrieval seed, epoch, or record order.

The current BioSensIA-DC LMDB helper writes pickled dictionaries under dense numeric ASCII keys. The new exporter should reuse that helper or preserve the same storage convention.

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
├── rcsb_download.py
├── rcsb_metadata.py
├── structure_mapping.py
├── citation_enrichment.py
├── quality.py
├── sidecars.py
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

Recommended optional additions:

```text
tenacity
rapidfuzz
```

Use `gemmi` for legacy PDB and PDBx/mmCIF parsing. Use RDKit for ligand chemistry. Use `httpx` for cached HTTP requests. Use Polars and PyArrow for sidecar generation.

All new dependencies must be pinned or bounded in `pyproject.toml`.

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
│   └── rcsb/
│       ├── mmcif/
│       │   └── <pdb_id>.cif.gz
│       ├── data_api/
│       │   └── <pdb_id>.json
│       ├── graphql/
│       │   └── <request_hash>.json
│       ├── crossref/
│       │   └── <request_hash>.json
│       └── pubmed/
│           └── <request_hash>.json
└── processed/
    └── pdbbind_2020r1/
        └── <run_id>/
            ├── manifest.json
            ├── config.resolved.toml
            ├── sidecars/
            ├── lmdb/
            ├── reports/
            ├── logs/
            └── checkpoints/
```

The `<run_id>` must include a short configuration hash:

```text
pdbbind2020r1-r6.0-v1-<config_hash_8>
```

Do not use only a timestamp as the run identifier. A timestamp may be appended for readability, but the configuration hash is required.

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
rcsb_cache_dir = "data/cache/rcsb"
output_root = "data/processed/pdbbind_2020r1"

[pocket]
distance_cutoff_angstrom = 6.0
distance_uses_ligand_heavy_atoms = true
include_protein_hydrogens = false
include_polymer_hetatm = true
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

[ligand]
primary_format = "sdf"
fallback_format = "mol2"
sanitize = true
allow_multiple_components = true
require_3d_coordinates = true

[rcsb]
download_mmcif = true
download_compressed = true
timeout_seconds = 60
maximum_retries = 5
requests_per_second = 2.0
use_data_api = true
use_graphql = true

[bibliography]
enabled = true
query_crossref = true
query_pubmed = true
allow_pdbbind_page_lookup = true

[lmdb]
map_size = 1099511627776
overwrite = false
include_quality_tiers = ["A", "B"]
```

Every resolved configuration value must be copied into:

```text
config.resolved.toml
manifest.json
```

---

## 8. Pipeline execution model

The pipeline must consist of resumable stages.

```text
inventory
parse-index
download-rcsb
parse-structures
extract-pockets
compare-pockets
map-structures
enrich-citations
quality-control
write-sidecars
export-lmdb
validate
report
```

Each stage must:

1. Read immutable inputs or outputs from preceding stages.
2. Write its result atomically.
3. Produce a completion marker containing the stage configuration hash.
4. Be safely restartable.
5. Refuse to reuse a checkpoint whose configuration hash is incompatible.
6. Record failures as data instead of terminating the entire run, unless `--fail-fast` is specified.

Intermediate checkpoints may use Parquet files under:

```text
checkpoints/<stage-name>/
```

Temporary files must be written under the run directory and renamed atomically only after successful completion.

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
  "started_at_utc": "...",
  "completed_at_utc": null,
  "status": "running",
  "python_version": "...",
  "platform": "...",
  "dependency_versions": {},
  "dataset": {
    "name": "PDBBind",
    "nominal_version": "2020R1",
    "distribution_revision": "Aug 2025",
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

The distribution revision should be read from the index header when possible, not hardcoded.

### 9.2 Source file inventory

Create:

```text
sidecars/source_files.parquet
```

One row per relevant source or downloaded file.

Required columns:

| Column               | Type               | Description                                                                                                                             |
| -------------------- | ------------------ | --------------------------------------------------------------------------------------------------------------------------------------- |
| `source_file_id`     | string             | Stable hash-based identifier                                                                                                            |
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

The manifest must include the checksums of all index files and all final outputs. It may refer to `source_files.parquet` rather than embedding all complex-file checksums directly.

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
* original line text.

The `*.pdf` token after `//` must be discarded and must not be persisted in any table, manifest field, log message, or LMDB record.

The parser may recognize it transiently only to remove it from the remaining text.

### 10.3 Complex identifier

Use:

```text
complex_id = "pdbbind-2020r1:<lowercase_pdb_id>"
```

Example:

```text
pdbbind-2020r1:2l3r
```

The parser must reject duplicate PDB IDs unless the duplicates are byte-for-byte equivalent. Conflicting duplicate records must receive a fatal index-validation issue.

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
Ka=...
```

Support:

* measurement types such as `Kd`, `Ki`, `IC50`, `EC50`, `Ka`, and unknown types;
* relations `=`, `<`, `>`, `<=`, `>=`, and approximate values;
* units `M`, `mM`, `uM`, `µM`, `nM`, `pM`, and `fM`;
* whitespace variants;
* scientific notation;
* unparseable raw values.

Normalize concentration measurements to molar units:

$$
v_{\mathrm{M}} = v \times f_{\mathrm{unit}}
$$

When a positive normalized molar value exists, calculate:

$$
pX = -\log_{10}(v_{\mathrm{M}})
$$

The relation must be reversed when transformed to the (pX) scale:

```text
Kd < 10 µM  →  pKd > 5
Kd > 10 µM  →  pKd < 5
```

Never replace the raw measurement string with the normalized value.

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
| `p_measurement_name`          | nullable string  |
| `p_relation`                  | nullable string  |
| `p_value`                     | nullable float64 |
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

A usable ligand requires:

* at least one atom;
* at least one heavy atom;
* finite three-dimensional coordinates;
* matching atom and coordinate counts.

Sanitization failure is a warning if the geometry and element identities remain usable. It becomes fatal only if neither format supplies a usable ligand.

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
* coordinate RMSD after atom mapping when feasible.

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

Do not remove salts or disconnected components from the bound ligand geometry. Record them and optionally derive a separate standardized parent identifier for future ligand matching.

### 13.5 Ligand geometry hash

Create a deterministic ligand geometry hash from:

* selected source-file SHA-256;
* ordered atomic numbers;
* ordered coordinates encoded as little-endian float64 or float32 according to a documented schema;
* component and conformer index;
* ligand-parser schema version.

The coordinate serialization and floating-point type must be fixed and versioned.

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

Resolve alternate locations deterministically per atom identity.

Default precedence:

1. highest occupancy;
2. blank alternate-location identifier;
3. alternate location `A`;
4. lexicographically smallest alternate-location identifier.

Record the number of discarded alternate-location records and any zero-occupancy selections.

### 14.4 Element determination

Use the PDB element field when present.

When absent, infer the element conservatively from the atom name and residue context. Record every inferred element.

Do not confuse atom names such as `CA` for an alpha carbon with calcium.

### 14.5 Protein versus non-protein atoms

By default:

* retain polymer protein atoms;
* exclude water;
* exclude non-polymer ligands, ions, and cofactors from the DrugCLIP pocket tensor;
* record nearby non-protein components separately;
* permit modified amino-acid residues represented with `HETATM` only when RCSB/mmCIF mapping identifies the residue as part of a polymer entity.

If RCSB metadata is not yet available, initially parse all records but defer final polymer classification until the mapping stage.

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

Define a canonical `pdbbind_atom_key` from the stable identity fields. Do not use only atom serial numbers.

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

Do not discard insertion codes.

---

## 15. Pocket extraction

### 15.1 Biological residue view

A protein residue belongs to the biological pocket when at least one of its heavy atoms lies within the configured distance cutoff of at least one ligand heavy atom.

For residue (r), include it when:

$$
\min_{a \in r,\ \ell \in L}
\lVert x_a-x_\ell\rVert
\leq d
$$

where:

* (L) is the selected ligand heavy-atom set;
* (d) is `distance_cutoff_angstrom`.

### 15.2 DrugCLIP atom view

A protein heavy atom belongs to the DrugCLIP pocket when:

$$
\min_{\ell \in L}
\lVert x_a-x_\ell\rVert
\leq d
$$

This atom-level set is the geometry exported to the LMDB.

The biological residue view and DrugCLIP atom view must both be retained in sidecars.

### 15.3 Distance computation

Requirements:

* use ligand heavy atoms only by default;
* use protein heavy atoms only for the LMDB view;
* use vectorized NumPy computation or a spatial index;
* avoid constructing an unnecessarily large full distance matrix for extreme structures;
* calculate each selected protein atom’s minimum ligand distance;
* preserve unrounded distances internally;
* store distances in Ångström;
* reject NaN or infinite coordinates.

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

The exact ordering fields and null ordering must be fixed in the schema version.

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

This prevents the original DrugCLIP loader from performing its sampled crop because the exported record will already be at or below the configured maximum. The original loader uses a default maximum of 256 and samples when the pocket is larger.

### 15.7 Pocket geometry hash

Calculate a SHA-256 hash from a canonical binary serialization containing:

* extraction schema version;
* PDBBind protein-file SHA-256;
* selected ligand geometry hash;
* model policy and selected model;
* alternate-location policy;
* distance cutoff;
* hydrogen policy;
* protein/polymer classification policy;
* deterministic crop settings;
* ordered atom identities;
* ordered element symbols;
* ordered coordinates;
* ordered minimum-ligand distances.

The hash must change if any input or extraction parameter capable of changing the exported geometry changes.

### 15.8 Pocket instance identifier

Use:

```text
pocket_instance_id =
    "pb20r1:<pdb_id>:<first_16_hex_chars_of_geometry_hash>"
```

Example:

```text
pb20r1:2l3r:8f3a25b19c173e62
```

Do not include an RCSB ligand-instance mapping in this identifier because enrichment may be unresolved or may improve later. The identifier must depend only on source geometry and versioned extraction rules.

---

## 16. PDBBind-provided pocket comparison

### 16.1 Source

Parse:

```text
<pdb_id>_pocket.pdb
```

using the same model, alternate-location, element, and hydrogen policies used for the re-extracted pocket.

Do not crop either pocket before performing the primary comparison.

### 16.2 Comparison levels

Perform comparisons at:

1. chain level;
2. residue level;
3. atom-identity level;
4. coordinate level;
5. ligand-distance level.

### 16.3 Required summary metrics

Record:

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
* (95%) quantile of minimum-ligand distance;
* whether one atom set is a subset of the other;
* whether one residue set is a subset of the other.

The atom-set Jaccard similarity is:

$$
J(A,B)=\frac{|A\cap B|}{|A\cup B|}
$$

### 16.4 Atom matching

Primary matching must use the canonical atom identity.

When atom identities cannot be matched because PDBBind processing altered identifiers, permit a fallback match using:

* chain;
* residue number;
* insertion code;
* residue name;
* atom name;
* coordinate tolerance.

Fallback matches must be marked separately from exact identity matches.

### 16.5 Difference details

Store every atom present in only one representation in a normalized sidecar, not only as counts.

Required classifications:

```text
only_reextracted
only_pdbbind
matched_exact_identity
matched_fallback
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
PDBBIND_POCKET_EXTENDS_BEYOND_CUTOFF
REEXTRACTED_POCKET_EMPTY
```

Thresholds must be configurable and must affect quality flags, not geometry itself.

---

## 17. RCSB mmCIF caching

### 17.1 Download

For each distinct PDB ID, cache:

```text
data/cache/rcsb/mmcif/<lowercase_pdb_id>.cif.gz
```

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
* never overwrite a valid cached file unless `--refresh-cache` is specified.

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
* use replacement metadata only in explicitly marked supplemental fields.

### 17.4 Offline mode

With `--offline`:

* make no network requests;
* use cached files and API responses only;
* mark missing enrichment as unresolved;
* allow geometry processing to continue.

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

For every PDBBind contributing `auth_chain_id`:

1. Find matching polymer atom records in mmCIF using `auth_asym_id`.
2. Determine corresponding `label_asym_id` values.
3. Determine the polymer entity ID for each `label_asym_id`.
4. Determine whether the entity is a protein.
5. Retrieve UniProt mappings for the polymer entity.
6. Record all mappings; do not select one accession arbitrarily.

### 18.3 Ambiguity resolution

If a direct `auth_asym_id` mapping is ambiguous:

1. Match PDBBind residues to mmCIF residues using author residue number, insertion code, residue name, and atom name.
2. Compare coordinates under the configured tolerance.
3. If necessary, align the PDBBind processed protein to the RCSB asymmetric unit using shared alpha-carbon atoms and the Kabsch algorithm.
4. Reattempt residue and atom matching in the aligned frame.
5. Require the configured minimum number of alignment atoms.
6. Record RMSD and matched-atom count.
7. Leave the mapping unresolved when evidence remains insufficient.

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
pocket_geometry_hash
```

### 18.5 UniProt mappings

Store one row per:

```text
pocket chain
  × RCSB polymer entity
  × UniProt accession
```

Record:

* accession;
* isoform identifier when supplied;
* mapping source;
* mapped sequence range;
* mapping coverage;
* mapping status;
* whether multiple UniProt accessions map to the same entity.

Do not collapse isoforms unless a separate policy explicitly requests it.

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

Store:

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

Do not require successful ligand-instance mapping for a geometrically valid candidate to enter the default library. It affects metadata quality, not the PDBBind-derived geometry.

---

## 20. Bibliographic enrichment

### 20.1 Sources and priority

Enrich references from:

1. mmCIF `_citation` and `_citation_author`;
2. RCSB Data API or GraphQL;
3. PDBBind/PDBBind+ entry page, when accessible;
4. Crossref using DOI or bibliographic matching;
5. PubMed using PMID, DOI, title, author, journal, year, or ECitMatch;
6. manual overrides.

The mmCIF `_citation.id` value `primary` identifies the citation considered most pertinent by the depositor, but that does not prove that the article contains the affinity measurement.

The `_citation_author` category connects authors to citation records using `citation_id`.

Crossref exposes deposited scholarly metadata through a public REST API.

PubMed E-utilities support searching, fetching records, and bibliographic citation matching. Requests should identify the tool and an email address, and higher request rates should use an API key.

### 20.2 Legacy PDF token

Do not store:

```text
pdb_bind_reference_token
```

Do not store the `*.pdf` filename anywhere.

### 20.3 Citation identity

Prefer identifiers in this order:

1. DOI;
2. PMID;
3. normalized bibliographic fingerprint;
4. local citation ID scoped to the PDB entry.

A bibliographic fingerprint may contain:

```text
normalized title
publication year
journal
volume
first page
first author
```

### 20.4 Reference-status model

Every binding measurement must have one affinity-reference status:

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

An RCSB primary citation alone must not produce `exact_affinity_reference`.

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

Parquet sidecars are the authoritative metadata store. Every table must include a schema version in Parquet metadata.

### 21.1 `complexes.parquet`

One row per PDBBind protein–ligand index record.

Required columns:

```text
complex_id
pdb_id
pdbbind_nominal_version
pdbbind_distribution_revision
index_line_number
index_line_raw
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
quality_tier
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
ligand_geometry_hash
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
rcsb_ligand_match_status
rcsb_nonpolymer_entity_id
rcsb_label_asym_id
rcsb_auth_asym_id
rcsb_auth_seq_id
ccd_id
bird_id
warnings
```

### 21.3 `pockets.parquet`

One row per re-extracted pocket.

```text
pocket_instance_id
complex_id
ligand_instance_id
pdb_id
pocket_geometry_hash
extraction_schema_version
distance_cutoff_angstrom
selected_model_id
model_count
altloc_policy
hydrogen_policy
raw_atom_count
raw_heavy_atom_count
exported_atom_count
biological_residue_count
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
quality_tier
warning_codes
error_codes
lmdb_profile_memberships
```

### 21.4 `pocket_residues.parquet`

One row per biological-pocket residue.

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

One row per uncropped re-extracted pocket atom.

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

One row per re-extracted/PDBBind-pocket comparison.

```text
pocket_instance_id
pdbbind_pocket_file_id
comparison_status
reextracted_atom_count
pdbbind_atom_count
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
warning_codes
```

### 21.7 `pocket_atom_differences.parquet`

One row per compared atom.

```text
pocket_instance_id
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
rcsb_label_asym_id
rcsb_polymer_entity_id
entity_type
entity_description
organism_name
taxonomy_id
sequence_length
alignment_atom_count
alignment_rmsd
warnings
```

### 21.9 `chain_uniprot_mappings.parquet`

One row per chain/entity/UniProt association.

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
pdb_sequence_begin
pdb_sequence_end
uniprot_sequence_begin
uniprot_sequence_end
mapping_coverage
```

### 21.10 `citations.parquet`

One row per distinct citation.

```text
citation_id
pdb_id
mmcif_citation_id
citation_role
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

### 21.11 `citation_authors.parquet`

One row per citation author.

```text
citation_id
ordinal
author_name
orcid
source
```

### 21.12 `affinity_reference_links.parquet`

One row per binding-measurement/reference candidate.

```text
measurement_id
complex_id
citation_id
reference_status
confidence
evidence_sources
evidence_note
automatic_or_manual
verified_by
verified_at_utc
```

### 21.13 `nearby_nonprotein_components.parquet`

Record excluded but nearby components such as metals, cofactors, waters, and non-polymer residues.

```text
pocket_instance_id
component_type
component_id
auth_chain_id
auth_seq_id
element
atom_count
minimum_ligand_distance
included_in_drugclip_tensor
exclusion_reason
```

### 21.14 `processing_issues.parquet`

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

Severity:

```text
info
warning
error
fatal
```

### 21.15 `lmdb_records.parquet`

One row per exported LMDB record.

```text
library_profile
lmdb_path
record_index
lmdb_key
pocket_instance_id
pocket_geometry_hash
atom_count
record_sha256
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

### 22.2 Suggested quality tiers

#### Tier A

A pocket qualifies for Tier A when:

* a ligand format was parsed successfully;
* a bound three-dimensional ligand geometry is available;
* the protein file parsed successfully;
* the pocket contains at least one supported protein heavy atom;
* all exported elements are supported by the configured DrugCLIP pocket dictionary;
* chain-to-RCSB mapping is unambiguous;
* at least one protein entity is identified;
* no fatal or severe geometry warning is present;
* deterministic cropping is either unnecessary or completed successfully.

#### Tier B

A pocket qualifies for Tier B when the geometry is usable but one or more nonfatal issues exist, such as:

* unresolved RCSB or UniProt mapping;
* SDF/MOL2 disagreement;
* deterministic cropping;
* moderate PDBBind-pocket disagreement;
* ambiguous ligand-instance mapping;
* citation/reference unresolved;
* modified protein residues;
* multiple contributing chains.

Multiple chains alone must not be treated as an error.

#### Tier C

Tier C includes geometrically usable but unusual or potentially misleading cases, such as:

* covalent ligand;
* strong dependence on an excluded metal or cofactor;
* substantial missing atoms;
* severe PDBBind-pocket disagreement;
* very small pocket;
* unusual ligand chemistry;
* unsupported protein atom removed by an explicit policy;
* low-confidence mapping after structural alignment.

Tier C pockets must remain in sidecars. Inclusion in an LMDB profile must be explicit.

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
* a geometry hash cannot be generated;
* the pocket identifier is duplicated with a different geometry;
* a fatal source-integrity issue exists.

Bibliographic or UniProt enrichment failure is not, by itself, a geometry rejection condition.

---

## 23. DrugCLIP dictionary compatibility

Read the active pocket dictionary from:

```text
external/DrugCLIP/data/dict_pkt.txt
```

or a configured equivalent.

Before export:

1. Remove hydrogen atoms.
2. Verify that every exported element token exists in the dictionary.
3. Record unsupported elements in sidecars.
4. Do not silently map elements such as selenium to sulfur.
5. Permit explicit versioned element mappings only through configuration.
6. Include the dictionary SHA-256 in the manifest and LMDB profile metadata.

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

The current BioSensIA-DC helper already writes this layout using pickled dictionaries.

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
record["pocket_coordinates"].shape == (N, 3)
len(record["pocket_atoms"]) == N
1 <= N <= max_pocket_atoms
no hydrogen atoms
all coordinates finite
all tokens supported by dict_pkt.txt
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

The required default profile is Tier A plus Tier B unless configured otherwise.

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
* quality tier;
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
    --run-dir data/processed/pdbbind_2020r1/<run_id> \
    --profile tiers-ab

python -m biosensia_pocket_library.cli validate \
    --run-dir data/processed/pdbbind_2020r1/<run_id>
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

Do not log the discarded PDBBind PDF filename token.

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

The final output must be byte-for-byte reproducible where library versions permit it, except for manifest timestamps and nondeterministic compression metadata. If Parquet byte identity cannot be guaranteed, table contents and row ordering must still be deterministic.

---

## 29. Validation

### 29.1 Input validation

Check:

* all expected index files exist;
* the declared PL count in the index header is parsed;
* parsed data-line count matches the declared count;
* each PL index record has no more than one discovered complex directory;
* source-file checksums are recorded.

### 29.2 Sidecar validation

Check:

* primary-key uniqueness;
* foreign-key integrity;
* no duplicate `pocket_instance_id`;
* no duplicate geometry hash with conflicting metadata;
* every accepted pocket has a ligand and complex;
* every exported pocket has atom rows;
* every chain mapping points to a known pocket;
* every affinity-reference link points to known measurement and citation rows;
* no legacy PDF token field exists;
* no column name contains `pdb_bind_reference_token`.

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
assert record["pocket_coordinates"].ndim == 2
assert record["pocket_coordinates"].shape[1] == 3
assert len(record["pocket_atoms"]) == len(record["pocket_coordinates"])
assert np.isfinite(record["pocket_coordinates"]).all()
```

Verify:

* keys are dense numeric ASCII strings;
* key count equals sidecar profile count;
* record order agrees with `lmdb_records.parquet`;
* each record geometry hash can be regenerated from sidecars;
* every token exists in the configured pocket dictionary;
* no record exceeds `max_pocket_atoms`.

### 29.4 DrugCLIP integration test

Load the produced LMDB through the actual BioSensIA-DC target-fishing pocket loader.

The integration test must:

1. Load at least one record.
2. Encode at least one pocket with the configured checkpoint.
3. Verify that no dictionary, shape, dtype, hydrogen, or cropping error occurs.
4. Verify that the pocket name returned by retrieval is the complete `pocket_instance_id`.
5. Join that ID to `pockets.parquet`.

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
* Tier A, B, and C counts;
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
* binding-measurement parsing;
* unit normalization;
* inequality reversal on the (pX) scale;
* ligand SDF selection;
* MOL2 fallback;
* SDF/MOL2 disagreement;
* alternate-location selection;
* PDB element inference;
* heavy-atom filtering;
* distance-based pocket extraction;
* biological residue extraction;
* multi-chain pockets;
* deterministic cropping;
* deterministic atom ordering;
* geometry hashing;
* pocket comparison;
* chain identifier mapping;
* RCSB cache validation;
* citation-author joins;
* reference-status rules;
* sidecar schema validation;
* LMDB record serialization.

### 31.2 Synthetic fixtures

Create small synthetic PDB, SDF, MOL2, and mmCIF fixtures that cover:

* one chain and one ligand;
* two chains contributing to one pocket;
* alternate locations;
* insertion codes;
* modified residues;
* nearby metal;
* multiple ligand components;
* more than 256 pocket atoms;
* unresolved RCSB mapping;
* NMR-style multiple models.

Do not commit actual PDBBind structures unless the license clearly permits redistribution.

### 31.3 Local-data integration tests

Provide tests marked:

```python
@pytest.mark.pdbbind_data
```

They should skip automatically when the local PDBBind files are unavailable.

Include a development test for a user-specified PDB ID such as:

```bash
pytest -m pdbbind_data --pdb-id 2l3r
```

### 31.4 Regression test

For a fixed local input and configuration:

* store expected counts and hashes in a local, noncommitted regression artifact;
* verify that rerunning produces identical pocket IDs and geometry hashes;
* verify that changing the distance cutoff changes the configuration hash and affected pocket hashes.

---

## 32. Acceptance criteria

The implementation is complete when:

1. A development run for one PDB ID completes from raw files to LMDB.
2. A full run parses every record in `INDEX_general_PL.2020R1.lst`.
3. Every raw input and output has a checksum.
4. The distance cutoff is configurable.
5. Pocket extraction is deterministic.
6. Ligand SDF is preferred and MOL2 fallback works.
7. PDBBind-provided pockets are compared and detailed atom differences are stored.
8. Contributing protein chains are identified.
9. RCSB mmCIF files are downloaded and cached.
10. PDBBind chains are mapped to RCSB chain instances and polymer entities when possible.
11. UniProt mappings preserve multiplicity and ambiguity.
12. Bibliographic metadata is stored without the PDBBind PDF filename token.
13. Binding measurements are parsed, normalized, and retained.
14. Affinity-reference statuses are implemented.
15. RCSB geometry never replaces PDBBind model geometry.
16. Parquet sidecars pass primary-key and foreign-key validation.
17. LMDB records use dense numeric ASCII keys and pickled dictionaries.
18. Every LMDB record contains a stable unique pocket identifier.
19. The produced LMDB loads successfully in the BioSensIA-DC target-fishing pipeline.
20. A final manifest and summary report are generated.

---

## 33. Suggested implementation sequence

Implement in the following pull-request-sized stages.

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
* biological residue view;
* DrugCLIP atom view;
* deterministic cropping;
* geometry hashes;
* PDBBind-pocket comparison;
* pocket atom and residue sidecars.

### Stage 4: RCSB enrichment

Implement:

* mmCIF download/cache;
* RCSB API cache;
* auth-to-label chain mapping;
* polymer entity mapping;
* UniProt mapping;
* ligand-instance mapping;
* coordinate-alignment fallback.

### Stage 5: Bibliographic enrichment

Implement:

* mmCIF citation extraction;
* citation-author extraction;
* Crossref/PubMed enrichment;
* PDBBind-page enrichment where available;
* affinity-reference statuses;
* manual overrides.

### Stage 6: Quality and LMDB export

Implement:

* quality tiers;
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
* full PDBBind 2020R1 run.
