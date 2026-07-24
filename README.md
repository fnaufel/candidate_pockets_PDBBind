# BioSensIA candidate-pocket libraries

This package builds reproducible, auditable DrugCLIP-compatible candidate-pocket
LMDBs from the licensed PDBbind 2020/v2024-reprocessed special distribution and
DrugCLIP's vendored `pdb/combine_set` bundles. It retains ligand, protein,
geometry, comparison, quality, and enrichment provenance in explicit Parquet
sidecars.

## PDBbind re-extracted library

```bash
uv sync
uv run biosensia-pocket-library check-drugclip-contract \
  --config config/pdbbind-pocket-library.toml

uv run biosensia-pocket-library build \
  --config config/pdbbind-pocket-library.toml \
  --pdb-id 2tpi --offline
```

For the full design, algorithms, commands, output contract, and troubleshooting guidance, see [`pdbbind_pocket_library.md`](pdbbind_pocket_library.md).

## DrugCLIP combine_set library

The second builder uses trusted `data.pkl` files under DrugCLIP's
`data/pdb/combine_set` as authoritative pocket geometry while using adjacent
SDF, MOL2, and PDB files for ligand identity and structure mapping.

```bash
uv run biosensia-pocket-library inventory-combine-set \
  --config config/drugclip-combine-set-library.toml --limit 10

uv run biosensia-pocket-library build-combine-set \
  --config config/drugclip-combine-set-library.toml --pdb-id 2h2h
```

Pickle deserialization can execute code. The dedicated configuration explicitly
authorizes these vendored pickles; do not enable that option for untrusted data.
See the detailed Quarto guide
[`drugclip_combine_set_library.qmd`](drugclip_combine_set_library.qmd), or its
concise Markdown companion
[`drugclip_combine_set_library.md`](drugclip_combine_set_library.md), for the
revised design and data contract.

The raw and derived coordinate-bearing data are intentionally ignored by Git. They inherit access and redistribution constraints from PDBbind and must not be redistributed without authorization.
