# External Data

No patient data or large public benchmark downloads are committed here.

Recommended layout when `RDKID_DATA_DIR` points to this directory:

```text
data/
├── private/                       # gated child bundle; never commit
└── open/
    ├── clinvar_temporal/
    ├── phenopacket_store/
    ├── varpb/
    └── opencravat_modules/
```

The gated SageBio challenge dataset is available under the challenge access terms from <https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data>. Public sources, versions, and hashes are recorded in [`docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md).

Use `rdkid doctor --profile wgs` before downloading or processing the full bundle. The downloader recommends 150 GB of free space for the approximately 85 GB compressed input plus working room.

