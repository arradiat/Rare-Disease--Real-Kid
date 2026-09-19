# Change log

## 2026-09-19: judge-ready repository reorganization

- Converted the ZIP snapshot into an installable `src/` package with the unified `rdkid` CLI.
- Moved tests, checkpoints, results, figures, reports, methods forms, and manifests into named subdirectories.
- Replaced machine-specific data defaults with `RDKID_DATA_DIR` and centralized repository paths.
- Preserved all 131 original files in a SHA-256 migration ledger and removed only two hash-identical duplicates.
- Converted the temporal checkpoint to the new module namespace and added exact fixed-seed prediction-parity coverage.
- Rewrote the README and active documentation around judge navigation, model comparisons, generalization evidence, and reproducibility levels.
- Restored a green 30-test suite, guarded submission verification, Markdown-link checks, environment diagnostics, and a successful LaTeX build.

## 2026-09-08: post-submission full-background extension

- Froze the 7 September report and maximum-score Track 1 submission as the baseline before beginning the extension.
- Downloaded only the 96.7 MB VarPB v3 background and two small 2023 causal tables; verified all checksums.
- Installed OpenCRAVAT 3.1.1 with hg38 1.47.2, gnomAD 2.2.0/data v2.1, REVEL 2024.07.30/data v4.3a, and AlphaMissense 1.1.1. CADD was excluded because its module alone is about 86.6 GB.
- Deduplicated 17,872,513 sample rows into 894,327 variants and retained 100% of the 2,577 causal keys, with explicit fallback values for two mapper failures.
- Ran all eight saved models on 275,739 VarPB spike-in trials. The original small SBNN was best of those models but remained weak at median rank 522 and Top-100 0.20%.
- Added coverage-preserving REVEL and AlphaMissense comparators. REVEL reached median rank 8.5 and Top-20 65.90%; AlphaMissense reached median rank 29 and Top-20 44.24%.
- Downloaded and checksum-verified the exact ClinVar 26 December 2023 archive used by the original temporal experiment.
- Deterministically selected 146,810 unambiguous labeled missense SNVs (54,548 pathogenic, 92,262 benign) and annotated 146,786 with the same predictor versions; 24 mapper failures remain explicit missing values.
- Added a Bayesian pairwise RankNet objective that scores each positive-negative pair under one shared posterior draw.
- Screened the full 24-cell factorial grid over ranks, KL weight, posterior initialization scale, and prior scale on training-only data, then retrained the selected rank `(50, 50, 1)` configuration on 117,700 variants.
- Kept VarPB completely out of training and selection. On its locked test, the augmented SBNN reached median rank 14, Top-20 65.70%, and Top-100 86.63%, making it competitive with REVEL for shortlist recall while remaining weaker at Top-1 and MRR.

## 2026-09-07: patient-level generalization contract

- Added an evaluator that separates candidate-generation recall from conditional reranking and counts missing causal hypotheses as end-to-end failures.
- Added Top-1/5/10, MRR, patient-bootstrap intervals, calibration, selective prediction, runtime, inheritance, unseen-gene and unseen-disease metrics.
- Exported all 719 post-2023 Phenopacket cases and every model prediction in a compressed, auditable table.
- Added fixed 20%/40% phenotype deletion and 20% phenotype-noise tests without retuning HPO weights.
- Added a bounded-memory VarPB scorer for all eight model families and documented the public 107-genome follow-up.

## 2026-09-07: public repository release

- Recorded the live Hugging Face result for the selected wide rank-50 SBNN: Submission 2, full match at rank 1, 100 Rank Points, F-max 1.000, and EPCR threshold 0.3724070924857989.
- Added a repository-level overview of the complete genome-to-submission pipeline.
- Exposed seven fitted models and the temporal consensus as separate modules under `models/`.
- Added a model registry and tests tying report names to importable implementations.
- Added one canonical SBNN submission at `submission/model1_bonamame_wide-rank50-sbnn.csv`.
- Added a guarded exporter and machine-readable provenance manifest for the submission.
- Added public model, results, data, and reproducibility documentation.
- Removed private coordination notes and duplicate cluster-package files from the public release.
- Updated checkpoint inference so the ClinVar training archive is required only for retraining.
- Added the selected SBNN output to exact Track 1 score verification.
- Verified 29 tests, the submission schema, rank 1, 100 Rank Points, and F-max 1.000.

## 2026-09-06: wide Bayesian comparison and generalization

- Evaluated a 50-200-100-1 network so low-rank choices were meaningful at the tested scale.
- Screened coupled rank, KL, posterior-scale, and prior-scale configurations.
- Added dense Bayes-by-Backprop and width-matched comparisons.
- Selected the rank `(50, 50, 1)` SBNN for its temporal top-1 ranking and PR-AUC.
- Retained the five-network ensemble as the strongest tested unseen-gene detector.
- Added repeated post-2023 Phenopacket/HPO ranking experiments.

## 2026-09-05: temporal benchmark and Track 1 transfer

- Added a leakage-controlled ClinVar temporal split with exact Variation ID exclusion.
- Added gradient boosting, dense MLP, ensemble, SBNN, BBB, and consensus baselines.
- Added classification, ranking, calibration, selective-prediction, and OOD metrics.
- Transferred the trained allele models to strict and broad child candidate tables.
- Verified that the accepted `BUB1B` pair remains rank 1 across the tested adapters.

## 2026-09-02: Bayesian pipeline foundation

- Added the shared candidate-pair feature contract.
- Added the conditional-logit linear Bayesian baseline and Laplace posterior.
- Added the low-rank neural Bayesian implementation and Monte Carlo inference.
- Added independent submission-schema validation and exact score verification.
- Preserved the original genome-wide pipeline and leaderboard-confirmed deterministic CSV.
