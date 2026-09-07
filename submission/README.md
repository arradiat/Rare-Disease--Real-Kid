# Track 1 submission


**[`model1_bonamame_wide-rank50-sbnn.csv`](model1_bonamame_wide-rank50-sbnn.csv)**

It is the canonical output of the selected wide rank-50 SBNN, not a renamed deterministic result. The first row is the accepted compound-heterozygous `BUB1B` pair. The file has six rows and the required columns:

```text
proband_id, chrom_1, pos_1, ref_1, alt_1,
chrom_2, pos_2, ref_2, alt_2, epcr, finding_type, notes
```

## Live Hugging Face result

The platform accepted this file as Submission 2 for `PROBAND01` under team `bonamame` and account `mametoure`.

| Metric | Result |
|---|---:|
| Full-match rank | 1 |
| Rank Points | 100.0 / 100 |
| F-max | 1.000 |
| F-max EPCR threshold | 0.3724070924857989 |

This threshold is exactly the EPCR of the first-ranked `BUB1B` pair in the committed CSV and matches the independent local evaluator.


The companion [`model1_bonamame_wide-rank50-sbnn.manifest.json`](model1_bonamame_wide-rank50-sbnn.manifest.json) records the model configuration, training provenance, input prediction hash, output hash, local score audit, and live platform result. Regenerate both files with:

```bash
python generate_track1_submission.py
```

The exporter accepts only the audited `strict_wide_rank50_sbnn.csv` prediction produced by `verify_temporal_models_track1.py`. It fails instead of publishing a CSV when the schema or exact score check does not pass.

The complete generation chain is:

```text
generate_track1_submission.py
    -> models/wide_rank50_sbnn.py
    -> verify_temporal_models_track1.py
    -> results/track1_temporal_transfer_20260906/strict_wide_rank50_sbnn.csv
    -> export_track1_sbnn_submission.py
    -> submission/model1_bonamame_wide-rank50-sbnn.csv
```

`generate_track1_submission.py` is the public entry point. `models/wide_rank50_sbnn.py` contains the selected model configuration. `verify_temporal_models_track1.py` performs posterior inference and pair scoring. `export_track1_sbnn_submission.py` performs final validation and writes the upload file.
