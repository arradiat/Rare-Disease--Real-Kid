# File Migration

The source was a ZIP-style snapshot without `.git` history. Before reorganization, a complete sibling backup was created as `Rare-Disease--Real-Kid-main.pre-reorganization-backup-20260919`.

[`original_file_manifest.json`](../artifacts/manifests/original_file_manifest.json) records all 131 original files with their source path, SHA-256, size, category, destination, and migration action.

## Policies

- Source code moved into the installable `src/rare_disease_real_kid` package.
- Tests moved into unit, integration, and smoke groups.
- Generated results, figures, checkpoints, reports, and methods forms moved under `artifacts/`.
- Original documentation replaced for readers was retained under `docs/archive/`.
- The September 8 33-page PDF is canonical; the distinct September 7 26-page PDF is archived.
- The second September 7 PDF and duplicate submission CSV were removed only after identical SHA-256 hashes were confirmed.
- The legacy temporal checkpoint is preserved, and a converted checkpoint with new module paths is provided.

Run `rdkid verify` to check destination coverage, backup hashes when the sibling backup is available, root cleanliness, documentation links, model imports, checkpoint loading, and submission invariants.

