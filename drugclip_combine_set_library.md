# DrugCLIP combine_set candidate-pocket library

## Revised recommendation

The implementation uses a separate source adapter and build orchestrator, not a
fork of the complete PDBbind builder and not source-condition branches inside
the existing geometry extractor. Both builders converge on the same sidecar,
LMDB, hashing, validation, reporting, and DrugCLIP-contract layers.

The real source audit changes one detail from the initial design proposal. The
pickles contain DrugCLIP **input-stage** pocket geometry: they commonly include
hydrogen and commonly exceed 256 atoms. DrugCLIP target fishing subsequently
applies `AffinityPocketDataset`, `RemoveHydrogenPocketDataset`,
`CroppingPocketDataset`, and normalization. Therefore the candidate LMDB must
preserve the raw pickle atom sequence and coordinates. Pre-filtering or
pre-cropping these records would change the established loader behavior.

The two representations are explicit:

- `post_transform`: existing PDBbind re-extractions, already hydrogen-free and
  bounded by `pocket.max_pocket_atoms`.
- `source_pickle`: raw combine_set geometry; hydrogen and more than 256 atoms
  are valid because the linked DrugCLIP loader performs those transforms.

`lmdb_records.record_representation` records this distinction. The profile
metadata fingerprints the loader files that define the downstream transform.

## Source and identity policy

`data.pkl` is authoritative for the LMDB atom order and coordinates. The
neighboring files have narrower roles:

- ligand SDF/MOL2: chemical identity, components, and bound ligand geometry;
- pocket and pocket6A PDB: residue and chain metadata after element-coordinate
  mapping;
- protein PDB: retained source provenance and future enrichment input.

Mapping cannot alter, reorder, add, or remove pickle atoms. An exact ordered
match is preferred. Otherwise the mapper accepts unique element-and-coordinate
matches within the configured tolerance and records unresolved or ambiguous
atoms. Mapping quality is separate from geometry acceptance.

Complex, ligand, and pocket identifiers include source or derivation hashes.
Repeated geometry is not collapsed merely because records share a PDB ID or a
content hash.

## Trust boundary

Python pickle is code-executing serialization. A build fails before discovery
and deserialization unless `combine_set.trusted_pickles = true`, or the CLI is
given `--trust-pickles`. Inventory does not require trust because it only
enumerates and hashes files. The manifest records that trusted deserialization
was authorized and records every selected source checksum.

## Schema compatibility

Sidecar schema v2 adds source-neutral provenance without removing v1 fields:

- `geometry_origin` and `geometry_source_file_id`;
- `derivation_method` and source geometry counts;
- `source_atom_key`, `source_mapping_status`, and
  `included_in_lmdb_source`;
- source-neutral left and right geometry roles on comparison rows;
- `lmdb_records.record_representation`.

The reader and validator accept completed v1 sidecars when their schema differs
only by these additive v2 columns. New PDBbind and combine_set builds both write
v2.

## Commands

Audit selection and file availability without executing pickle payloads:

```bash
uv run biosensia-pocket-library inventory-combine-set \
  --config config/drugclip-combine-set-library.toml \
  --pdb-id 2h2h
```

Build sidecars only:

```bash
uv run biosensia-pocket-library build-combine-set-sidecars \
  --config config/drugclip-combine-set-library.toml \
  --pdb-id 2h2h
```

Build sidecars plus the default candidate LMDB:

```bash
uv run biosensia-pocket-library build-combine-set \
  --config config/drugclip-combine-set-library.toml \
  --pdb-id 2h2h
```

Selection is deterministic: directories are sorted by identifier before
`--limit` is applied. File contents, selected identifiers, semantic
configuration, and the DrugCLIP library contract all participate in run
identity.

## Validation

A completed run validates:

- required pickle fields, one finite `(n, 3)` pocket array, and atom alignment;
- transformed token coverage by the linked DrugCLIP pocket dictionary;
- source file checksums and relational joins;
- raw source geometry hashes regenerated from atom sidecars;
- dense numeric LMDB keys and record-to-sidecar order;
- little-endian contiguous float32 coordinates in serialized records;
- physical, logical, profile, and per-record checksums;
- explicit structure-mapping ambiguity and processing issues.
