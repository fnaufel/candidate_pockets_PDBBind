# BioSensIA PDBbind candidate-pocket library

This package builds a reproducible, auditable DrugCLIP-compatible candidate-pocket LMDB from the licensed PDBbind 2020/v2024-reprocessed special distribution. It retains the ligand, protein, extraction, comparison, quality, and enrichment provenance in explicit Parquet sidecars.

```bash
uv sync
uv run biosensia-pocket-library check-drugclip-contract \
  --config config/pdbbind-pocket-library.toml

uv run biosensia-pocket-library build \
  --config config/pdbbind-pocket-library.toml \
  --pdb-id 2tpi --offline
```

For the full design, algorithms, commands, output contract, and troubleshooting guidance, see [`pdbbind_pocket_library.md`](pdbbind_pocket_library.md).

The raw and derived coordinate-bearing data are intentionally ignored by Git. They inherit access and redistribution constraints from PDBbind and must not be redistributed without authorization.
