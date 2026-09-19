# Reproducibility

## Install and verify

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[report,test]"
rdkid doctor --profile core
rdkid verify
python -m pytest -q
```

The reference lightweight environment used Python 3.12.7. Exact recorded package versions are in [`requirements/report.txt`](../requirements/report.txt).

## Regenerate the canonical CSV

Without external inputs, the exporter revalidates the committed audited model output:

```bash
rdkid submission generate
```

With the compact candidate export:

```bash
rdkid submission generate \
  --ranking-zip /path/to/ranking.zip \
  --checkpoint artifacts/checkpoints/temporal/trained_models.pkl
```

The full 85 GB child bundle is not required for these commands.

## Retrain the temporal models

```bash
rdkid submission generate \
  --historical data/open/clinvar_temporal/clinvar_20231226.vcf.gz \
  --ranking-zip /path/to/ranking.zip \
  --checkpoint artifacts/checkpoints/temporal/retrained_models.pkl
```

The historical VCF is read only when the requested checkpoint does not already exist.

## Full WGS preflight and execution

```bash
rdkid doctor --profile wgs
rdkid pipeline all \
  --data-dir /path/to/gated_bundle \
  --reference /path/to/matching_GRCh38.fasta \
  --threads 32
```

The original WGS environment used `pysam==0.24.0`, CPU execution, six alignment threads, 2,000 read pairs per Python batch, 2,000,000 read pairs per worker shard, and a 10,000,000-base BWA chunk. Exact versions of pybwa, WhatsHap, NumPy, and SciPy from that cluster run were not retained and are not invented here.

## Public data versions and hashes

| Dataset | Version | SHA-256 |
|---|---|---|
| ClinVar GRCh38 | 2023-12-26 | `510c54e3fd04b0786f0a110bab63853e9e8fc3e9dde9e07ce7bc83f10cb63786` |
| ClinVar GRCh38 | 2026-08-29 | `e72f5619f003c80cb3051afcb0c1f11c4409741054b48bbfd0ced03046ec91f9` |
| Phenopacket Store | 0.1.27 | `d0c70005bb09b87035087516252b6b539fbc0415be7c664b73c9cde7eedf205a` |
| HPO OBO | 2023-10-09 | `d3c1280f2780628cc2adc41d43ad92eab8fa4c2bab0f421e6334c7b1e6cd4e52` |
| VarPB background | v3 | `4606c256319d42819934d225ff30a78b40256eeda31c57ccc960375c12bd42ba` |
| VarPB dominant causals | 2023 | `84f1b31adc45e7fda118d9b2b670a0ce5a6cb96ea2d12888a78165be161b6158` |
| VarPB recessive causals | 2023 | `914e8fced43c31554701336382e62cef67b296781e74e4b854e5a78b58e22d4d` |

## Build the report

```bash
rdkid doctor --profile report
rdkid report build
```

The build regenerates all report figures from committed result artifacts, runs the tests and submission checks, and writes `artifacts/reports/rebuilt_from_source.pdf`. The canonical 33-page submitted report remains available separately for exact historical reference.

## Reproducibility boundary

- Submission verification is self-contained.
- Checkpoint inference additionally needs `ranking.zip`.
- Public retraining needs versioned public downloads.
- Full WGS execution needs gated data and missing bioinformatics dependencies.
- Phenopacket/HPO and VarPB predictions, ranks, summaries, audits, and the augmented checkpoint are committed and hash-checked by `rdkid verify`.
