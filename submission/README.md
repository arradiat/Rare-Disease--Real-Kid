# Canonical Track 1 Submission

[`model1_bonamame_wide-rank50-sbnn.csv`](model1_bonamame_wide-rank50-sbnn.csv) is the upload-ready output of the selected wide rank-50 SBNN.

The live evaluator accepted it as Submission 2 for team `bonamame`, account `mametoure`, and proband `PROBAND01`: full-match rank 1, 100 Rank Points, and F-max 1.000. The adjacent [manifest](model1_bonamame_wide-rank50-sbnn.manifest.json) records the hashes, model configuration, provenance, local score, and live receipt.

Regenerate or revalidate it with:

```bash
rdkid submission generate
```

The guarded chain is:

```text
rdkid submission generate
  -> rare_disease_real_kid.models.wide_rank50_sbnn
  -> rare_disease_real_kid.evaluation.track1_transfer
  -> artifacts/results/track1_temporal_transfer_20260906/strict_wide_rank50_sbnn.csv
  -> rare_disease_real_kid.evaluation.submission_export
  -> submission/model1_bonamame_wide-rank50-sbnn.csv
```

The exporter refuses to publish when the source audit, schema, or locally reproduced score check fails.

