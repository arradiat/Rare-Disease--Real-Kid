# Reproducibility

## Reproduction levels

The repository supports three distinct levels of reproduction. Keeping them separate prevents the private 85 GB genome bundle from being mistaken for a requirement of the lightweight ML checks.

1. **Submission verification:** validates the committed CSV, recomputes Rank Points and F-max, and requires no external data.
2. **Checkpoint inference:** recreates all model-specific child rankings from the retained 2 MB checkpoint and Bonaventure's compact `ranking.zip` candidate export.
3. **Full training:** rebuilds temporal models from the archived ClinVar 2023 VCF. It does not require the child's FASTQ files.

The full WGS stage is separate. It requires the gated patient data and the cluster bioinformatics environment described in the report.

## CPU environment

The local analysis, plotting, and verification were run on a MacBook Pro `Mac14,7` with Apple M2, 8 CPU cores, 8 GB RAM, macOS 15.2 arm64, and Python 3.12.7 from Anaconda.

Recorded core packages:

| Package | Version |
|---|---:|
| NumPy | 1.26.4 |
| SciPy | 1.13.1 |
| scikit-learn | 1.5.1 |
| pandas | 2.2.2 |
| Matplotlib | 3.9.2 |
| joblib | 1.4.2 |
| threadpoolctl | 3.5.0 |
| pytest | 7.4.4 |
| Pillow | 10.4.0 |
| huggingface-hub | 0.36.0 |

The generated machine manifest is stored at [`../report/environment_manifest.json`](../report/environment_manifest.json). `requirements-report.txt` covers verification and figures; `requirements-bayesian.txt` and `requirements-sbnn.txt` contain model-specific dependencies.

## Verify the committed submission

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-report.txt

python audit_submission.py submission/model1_bonamame_wide-rank50-sbnn.csv
python verify_track1_exact_scores.py
python -m pytest -q
```

Expected test result: 34 passed. The exact score audit must report `full_match_rank: 1`, `rank_points: 100.0`, `f_max: 1.0`, and `exact_perfect_score: true` for `selected_wide_rank50_sbnn`.

## Recreate the CSV from the fitted checkpoint

The required compact input is `ranking.zip`, produced by the upstream WGS pipeline. It contains candidate TSV files, not FASTQ reads or a full VCF.

```bash
python generate_track1_submission.py \
  --ranking-zip /path/to/ranking.zip \
  --checkpoint results/track1_temporal_transfer_20260906/trained_models.pkl
```

This command regenerates strict and broad Track 1 outputs for every retained temporal model, selects only `strict_wide_rank50_sbnn.csv`, validates it, and publishes it under the upload-ready name `submission/model1_bonamame_wide-rank50-sbnn.csv`.

## Retrain the temporal models

The 2023 ClinVar archive is needed only if the checkpoint is absent:

```bash
python generate_track1_submission.py \
  --historical data/open/clinvar_temporal/clinvar_20231226.vcf.gz \
  --ranking-zip /path/to/ranking.zip \
  --checkpoint results/track1_temporal_transfer_20260906/retrained_models.pkl
```

When the checkpoint path does not exist, the script trains gradient boosting, a dense MLP, a five-network ensemble, the original SBNN, dense BBB, the wide rank-50 SBNN, and a width-matched dense BBB. It then saves the fitted collection and runs the child transfer.

## Public datasets

| Dataset | Version | Source | Purpose |
|---|---|---|---|
| ClinVar GRCh38 | 2023-12-26 | [NCBI archive](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/) | Training labels and training-time gene history |
| ClinVar GRCh38 | 2026-08-29 | [NCBI current archive](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) | Temporal holdout |
| Phenopacket Store | 0.1.27 | [Monarch release](https://github.com/monarch-initiative/phenopacket-store/releases/tag/0.1.27) | Solved patients, variants, genes, and HPO terms |
| Human Phenotype Ontology | 2023-10-09 | [HPO releases](https://github.com/obophenotype/human-phenotype-ontology/releases) | Leakage-controlled phenotype matching |
| Variant Prioritization Benchmark | v3 | [Zenodo](https://doi.org/10.5281/zenodo.10291114) | Executed 107-genome full-background spike-in evaluation |

## Recompute patient-level generalization

The executed case-level outputs are evaluated independently of the model code:

```bash
python patient_generalization_benchmark.py \
  --cases results/patient_generalization_protocol_20260907/patient_cases.csv \
  --predictions results/patient_generalization_protocol_20260907/patient_predictions.csv.gz \
  --output-dir results/patient_generalization_protocol_20260907/evaluation
```

See [`GENERALIZATION_BENCHMARK.md`](GENERALIZATION_BENCHMARK.md) for the complete Phenopacket reconstruction and VarPB scoring commands.

## Reproduce the OpenCRAVAT VarPB annotation

The separate `.venv-opencravat` environment uses OpenCRAVAT 3.1.1. Its module
directory is `data/open/opencravat_modules/`; installed versions are hg38 mapper
1.47.2, gnomAD 2.2.0 with data v2.1, REVEL 2024.07.30 with data v4.3a, and
AlphaMissense 1.1.1. The environment and modules use about 28 GB. CADD is not a
dependency of this run.

The complete preparation, annotation, export, scoring, and evaluation commands
are recorded in [`GENERALIZATION_BENCHMARK.md`](GENERALIZATION_BENCHMARK.md).
The resulting summary is
[`../results/varpb_2023_opencravat_gnomad/evaluation/summary.json`](../results/varpb_2023_opencravat_gnomad/evaluation/summary.json),
and the independent predictor comparison is
[`../results/varpb_2023_external_predictors/evaluation/summary.json`](../results/varpb_2023_external_predictors/evaluation/summary.json).

The predictor-augmented follow-up, including the exact 24-cell factorial grid,
training-only selection rule, pairwise loss, and locked test result, is recorded
step by step in
[`POST_SUBMISSION_GENERALIZATION_2026-09-08.md`](POST_SUBMISSION_GENERALIZATION_2026-09-08.md).
The executable entry point is `augmented_sbnn_varpb_experiment.py`; its `train`
subcommand does not accept a VarPB causal-table argument, preventing those test
labels from entering training or model selection.

Recorded SHA-256 hashes:

```text
510c54e3fd04b0786f0a110bab63853e9e8fc3e9dde9e07ce7bc83f10cb63786  clinvar_20231226.vcf.gz
e72f5619f003c80cb3051afcb0c1f11c4409741054b48bbfd0ced03046ec91f9  clinvar_20260829.vcf.gz
d0c70005bb09b87035087516252b6b539fbc0415be7c664b73c9cde7eedf205a  all_phenopackets_v0.1.27.zip
d3c1280f2780628cc2adc41d43ad92eab8fa4c2bab0f421e6334c7b1e6cd4e52  hp_2023-10-09.obo
4606c256319d42819934d225ff30a78b40256eeda31c57ccc960375c12bd42ba  allSamplesExonicVariants.txt.gz
84f1b31adc45e7fda118d9b2b670a0ce5a6cb96ea2d12888a78165be161b6158  2023_autosomalDominant_notInClinvarAtAllBefore.txt
914e8fced43c31554701336382e62cef67b296781e74e4b854e5a78b58e22d4d  2023_autosomalRecessive_notInClinvarAtAllBefore.txt
```

## Full report build

```bash
./build_report.sh
```

The script captures the environment, regenerates report figures, runs the tests, exports the canonical SBNN submission, recomputes exact challenge scores, and compiles the LaTeX report. If `ranking.zip` is available at `../ranking.zip`, it first recreates the temporal child outputs from the checkpoint. A different location can be supplied as `RANKING_ZIP=/path/to/ranking.zip ./build_report.sh`. Without that compact input, the build uses the committed audited SBNN prediction. A TeX installation with `pdflatex` and `bibtex` is required.

## Private challenge data

The approximately 85 GB SageBio bundle is gated and must not be committed. It was processed on Mila cluster storage. The original WGS environment used CPU execution, six alignment threads, 2,000 read pairs per Python batch, 2,000,000 read pairs per process-isolated shard, a 10,000,000-base BWA chunk, and `pysam==0.24.0`.

The retained local files do not establish the exact versions of pybwa, WhatsHap, OpenCRAVAT, NumPy, or SciPy used in that original cluster run. The report states this limitation rather than inventing versions.
