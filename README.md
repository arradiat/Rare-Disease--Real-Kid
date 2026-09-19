# Rare Disease, Real Kid

An end-to-end, uncertainty-aware variant-prioritization system for Track 1 of the 2026 MVA Rare Disease, Real Kid hackathon. The project turns whole-genome evidence into an auditable ranked diagnosis by combining biological filtering, phenotype matching, deterministic baselines, dense Bayesian models, deep ensembles, and low-rank Singular Bayesian Neural Networks (SBNNs).

The submitted **wide rank-50 SBNN** ranked the accepted compound-heterozygous `BUB1B` pair first and earned the maximum live Hugging Face score: **100 Rank Points and F-max 1.000**. The work goes substantially beyond the challenge case, with temporal ClinVar evaluation, 719 post-2023 Phenopacket cases, and 275,739 VarPB genome-background spike-in trials. On the hardest full-background test, the predictor-augmented pairwise SBNN reached **86.63% Top-100 recall**, exceeded REVEL at Top-100, and nearly matched it at Top-20 while retaining a compact Bayesian posterior.

**Team:** Mame Diarra Toure and Bonaventure F. P. Dossou  

## Project highlights

| Contribution | Result |
|---|---|
| Maximum challenge performance | Accepted `BUB1B` pair at rank 1; **100 Rank Points; F-max 1.000** |
| Strong temporal transfer | Wide rank-50 SBNN achieved **84.35% Top-1**, PR-AUC **0.9736**, and F-max **0.9235** |
| Patient-oriented generalization | **719** post-2023 solved Phenopackets evaluated across five candidate-set seeds |
| Full-background generalization | Augmented SBNN achieved **65.70% Top-20** and **86.63% Top-100** over **275,739** VarPB trials |
| Competitive with leading predictors | Within **0.20 percentage points** of REVEL at Top-20 and **1.81 points better** at Top-100 |
| Compact Bayesian integration | **56,204** variational scalars: 9.9% fewer than matched dense BBB and 64.0% fewer than a matched five-network ensemble |
| Reproducible evidence | Nine model implementations, two load-tested checkpoints, 132 hash-verified generalization artifacts, and 31 passing tests |

## Start here

Python 3.11 or 3.12 is supported.

```bash
python -m pip install -e ".[report,test]"
rdkid verify
python -m pytest -q
```

`rdkid verify` checks the submission schema and score, all registered models, both bundled checkpoints, the committed generalization-artifact hashes, and all internal documentation links. The canonical upload is [`submission/model1_bonamame_wide-rank50-sbnn.csv`](submission/model1_bonamame_wide-rank50-sbnn.csv), with SHA-256 `308eb684858fc9bfc97cd8b81f0689914408f78239a1b42563419498157510eb`.

## The scientific task

The challenge asks which pair of genetic variants best explains a child's disease. A whole genome contains millions of variants, so the system separates two responsibilities:

1. **Candidate generation** enforces genotype quality, rarity, molecular consequence, inheritance, phase, and phenotype compatibility.
2. **Probabilistic ranking** compares the remaining hypotheses and quantifies model uncertainty.

```text
gated child FASTQ + supplied VCF + phenotype
                    |
                    v
alignment, QC, normalization, annotation and phasing
                    |
                    v
rare, damaging, inheritance-compatible alleles
                    |
                    v
same-gene compound-heterozygous pairs + HPO similarity
                    |
                    v
deterministic / ensemble / BBB / low-rank SBNN comparison
                    |
                    v
wide rank-50 SBNN posterior draws for both alleles
                    |
                    v
pair score = phenotype_similarity * sqrt(P1 * P2)
                    |
                    v
audited six-row Track 1 submission
```

The Track 1 adapter deliberately excludes the submitted rank, previous pair score, challenge answer, ClinVar significance, CADD, REVEL, and SpliceAI. Candidate construction and phenotype similarity remain separate from learned pathogenicity scoring, reducing the risk of direct answer leakage.

## Final Track 1 result

| Check | Result | Evidence |
|---|---:|---|
| Accepted `BUB1B` pair rank | **1** | [Submission manifest](submission/model1_bonamame_wide-rank50-sbnn.manifest.json) |
| Rank Points | **100** | [Exact score audit](artifacts/results/audits/exact_track1_score_audit.json) |
| F-max | **1.000** | [Exact score audit](artifacts/results/audits/exact_track1_score_audit.json) |
| Submission rows | 6 | [Canonical CSV](submission/model1_bonamame_wide-rank50-sbnn.csv) |
| Schema and value audit | Pass | [Submission audit](artifacts/results/audits/submission_audit.json) |

The same accepted pair was ranked first by all eight submission-era learned models. The wide rank-50 SBNN was selected before submission because it had the strongest temporal variant-ranking result, not because alternatives were hidden after the challenge result was known.

## Why nine models were built

The project treats the challenge result as the beginning of the analysis rather than its endpoint. Nine related methods are implemented under shared data and evaluation contracts so that predictive accuracy, uncertainty, parameter efficiency, and robustness can be compared directly.

| Model | What it tests | Implementation |
|---|---|---|
| Gradient boosting | Strong deterministic tabular baseline | [`gradient_boost.py`](src/rare_disease_real_kid/models/gradient_boost.py) |
| Dense MLP | Single-network neural baseline | [`dense_mlp.py`](src/rare_disease_real_kid/models/dense_mlp.py) |
| Five-network deep ensemble | Model disagreement and OOD detection | [`deep_ensemble.py`](src/rare_disease_real_kid/models/deep_ensemble.py) |
| Small SBNN | Initial low-rank Bayesian baseline | [`small_sbnn.py`](src/rare_disease_real_kid/models/small_sbnn.py) |
| Dense BBB | Independent posterior over every dense parameter | [`full_rank_bbb.py`](src/rare_disease_real_kid/models/full_rank_bbb.py) |
| Wide rank-50 SBNN | Selected Track 1 model and scalable Bayesian comparison | [`wide_rank50_sbnn.py`](src/rare_disease_real_kid/models/wide_rank50_sbnn.py) |
| Wide dense BBB | Width-matched full-weight Bayesian comparison | [`wide_dense_bbb.py`](src/rare_disease_real_kid/models/wide_dense_bbb.py) |
| Temporal consensus | Complementarity of tree, ensemble, and SBNN predictions | [`temporal_consensus.py`](src/rare_disease_real_kid/models/temporal_consensus.py) |
| Augmented pairwise SBNN | Predictor integration and full-background ranking | [`augmented_pairwise_sbnn.py`](src/rare_disease_real_kid/models/augmented_pairwise_sbnn.py) |

The machine-readable model index is [`registry.py`](src/rare_disease_real_kid/models/registry.py). The submitted 50-200-100-1 SBNN uses ranks `(50, 50, 1)`, KL weight `0.01`, posterior-scale multiplier `0.1`, prior scale `0.3`, seed `20260905`, and 128 posterior draws.

## Beyond the child

The strongest part of this work is that it does not stop after solving the challenge case. The modeling ideas are tested at three complementary scales:

- **Challenge ranking:** recovering the accepted compound-heterozygous pair from the child's shortlist.
- **Temporal and patient transfer:** ranking later ClinVar variants and independent post-2023 solved Phenopackets.
- **Genome-background transfer:** finding causal variants among realistic backgrounds from 107 genomes.

### Temporal ClinVar

Models were trained on the verified 26 December 2023 ClinVar release and evaluated on later records excluded by Variation ID. The selected wide rank-50 SBNN delivered the strongest temporal ranking result: **84.35% Top-1 among 100 candidates, PR-AUC 0.9736, and F-max 0.9235**. The deep ensemble complemented this result with the strongest unseen-gene OOD discrimination, giving the system separate high-performing choices for ranking and shift detection.

Evidence: [`wide_rank_results.json`](artifacts/results/wide_rank_temporal_20260906/wide_rank_results.json), [`full_bbb_results.json`](artifacts/results/full_bbb_temporal_20260906/full_bbb_results.json), and [`consensus_holdout.json`](artifacts/results/clinvar_temporal_20260905/consensus_holdout.json).

### Post-2023 Phenopackets and HPO

The patient-oriented study contains **719 temporally new solved cases**, including **64 two-variant cases**, evaluated across five independently resampled ten-candidate cohorts. HPO conditioning improved ranking across every evaluated model. The dense BBB averaged **58.89% ± 0.96 percentage points Top-1**, while the wide dense BBB reached **68.13% mean Top-1 on the challenge-like two-variant subset**. These consistent gains show that molecular evidence and phenotype similarity contribute complementary information.

The committed evidence includes the case table, compressed prediction table, fixed-cohort robustness conditions, bootstrap intervals, calibration metrics, novelty strata, and every repeated-seed summary:

- [Patient cases](artifacts/results/patient_generalization_protocol_20260907/patient_cases.csv)
- [Patient predictions](artifacts/results/patient_generalization_protocol_20260907/patient_predictions.csv.gz)
- [Patient evaluation](artifacts/results/patient_generalization_protocol_20260907/evaluation/patient_generalization_results.json)
- [Five-seed HPO summary](artifacts/results/phenopacket_hpo2023_generalization_repeated_20260906/repeated_seed_summary.json)

### VarPB full-background ranking

The full-background study is the most demanding generalization experiment in the repository: **2,577 causal variants crossed with 107 real genome backgrounds, producing 275,739 locked spike-in trials**. Applying the frozen submission-era models first established that whole-genome prioritization requires richer protein-level evidence than compact shortlist ranking alone.

That finding motivated a new 54-feature pairwise SBNN trained only on pre-2024 ClinVar missense records. The model integrates REVEL and AlphaMissense with the original molecular, frequency, consequence, sequence, and missingness features; uses a gene-disjoint split and Bayesian RankNet loss; and never accesses VarPB causal answers during training or model selection.

| Locked VarPB result | Median rank | Top-1 | Top-20 | Top-100 | MRR |
|---|---:|---:|---:|---:|---:|
| Predictor-augmented pairwise SBNN | 14 | 7.26% | 65.70% | **86.63%** | 0.1596 |
| REVEL | **8.5** | **12.48%** | **65.90%** | 84.82% | **0.2516** |
| AlphaMissense | 29 | 8.84% | 44.24% | 69.72% | 0.1632 |

The result is particularly strong for shortlist recall:

- At **Top-20**, the SBNN reaches 65.70%, only **0.20 percentage points from REVEL's 65.90%** and 21.46 points above AlphaMissense.
- At **Top-100**, the SBNN reaches **86.63%**, outperforming REVEL by **1.81 percentage points** and AlphaMissense by **16.91 points**.
- Its median causal rank is **14**, compared with 29 for AlphaMissense.
- It produces posterior draws, predictive standard deviation, and mutual-information estimates from one fitted model—capabilities not supplied by the standard deterministic predictor outputs.
- Its 56,204 fitted variational scalars are fewer than matched dense BBB and far fewer than the matched five-network ensemble.

REVEL remains strongest at Top-1 and MRR, while the augmented SBNN is strongest at Top-100. Together, these results establish the SBNN as a competitive uncertainty-aware evidence integrator for clinically relevant shortlist construction.

Evidence: [augmented evaluation](artifacts/results/varpb_augmentation_20260908/evaluation/summary.json), [external predictors](artifacts/results/varpb_2023_external_predictors/evaluation/summary.json), [submission-era models](artifacts/results/varpb_2023_opencravat_gnomad/evaluation/summary.json), and the complete [post-submission protocol](docs/POST_SUBMISSION_GENERALIZATION.md).

### Frozen transfer of the augmented model

The augmented checkpoint also transferred back to Track 1 without retraining on the challenge child. It preserved the accepted pair at rank 1 in both strict and broad candidate sets and achieved the local official-equivalent maximum of **100 Rank Points and F-max 1.000**. This provides a valuable consistency result: the model's large-background gains did not sacrifice the original child-level solution. The selected wide rank-50 SBNN is the live-submitted model; the augmented result is locally evaluated.

Evidence: [Track 1 augmented audit](artifacts/results/track1_augmented_sbnn_20260908/track1_augmented_sbnn_audit.json) and [bundled checkpoint](artifacts/checkpoints/augmented_pairwise_sbnn/augmented_sbnn.pkl).

## Command-line interface

The project exposes one interface through both `rdkid` and `python -m rare_disease_real_kid`:

```bash
rdkid doctor --profile core
rdkid verify
rdkid submission generate
rdkid pipeline --help
rdkid benchmark clinvar --help
rdkid benchmark phenopacket --help
rdkid benchmark varpb --help
rdkid report build
```

Use `rdkid doctor --profile wgs` before attempting genome processing. It reports missing packages, executables, data, and disk capacity instead of treating unavailable stages as executed.

## Reproduction levels

| Level | Command | Required inputs |
|---|---|---|
| Verify committed results | `rdkid verify` | Core Python dependencies only |
| Run automated tests | `python -m pytest -q` | Test dependencies |
| Rebuild figures and report | `rdkid report build` | Reporting dependencies and TeX |
| Recreate checkpoint predictions | `rdkid submission generate --ranking-zip /path/to/ranking.zip` | Compact private candidate export |
| Retrain public benchmarks | `rdkid benchmark <clinvar\|phenopacket\|varpb> ...` | Versioned public datasets |
| Run the child WGS workflow | `rdkid pipeline all --data-dir /secure/path ...` | Gated child data, reference, bioinformatics tools, and substantial storage |

Exact public-data versions, hashes, commands, and environment requirements are in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md). The lightweight verification path is self-contained; retraining and full WGS execution intentionally require external inputs.

## Repository map

```text
.
├── README.md                         # Project overview and judge reading path
├── pyproject.toml                    # Package metadata, dependencies and rdkid command
├── CITATION.cff                      # Citation and authorship metadata
├── src/rare_disease_real_kid/
│   ├── cli.py                        # Unified command-line interface
│   ├── paths.py                      # Repository, artifact and data path resolution
│   ├── checkpoint.py                 # Trusted bundled-checkpoint loading
│   ├── pipeline/                     # WGS inspection, alignment, annotation and ranking
│   ├── features/                     # Stable allele and candidate-pair features
│   ├── models/                       # Nine model definitions and central registry
│   ├── experiments/                  # ClinVar, Phenopacket, VarPB and tuning workflows
│   ├── evaluation/                   # Submission, patient, Track 1 and VarPB evaluation
│   └── reporting/                    # Reproducible figures and environment capture
├── tests/
│   ├── unit/                         # Metrics, layers, registry and schema tests
│   ├── integration/                  # Cross-module export and ranking tests
│   └── smoke/                        # Checkpoint loading and prediction parity
├── artifacts/
│   ├── checkpoints/                  # Temporal and augmented fitted models
│   ├── results/                      # Predictions, audits, metrics and benchmark tables
│   ├── figures/                      # Figures regenerated from committed results
│   ├── reports/                      # Canonical report, archive and rebuilt report
│   ├── methods-forms/                # Submitted methods workbooks
│   └── manifests/                    # Integrity, environment and provenance records
├── submission/                       # Canonical upload CSV and provenance manifest
├── report/                           # LaTeX source and bibliography
├── docs/                             # Architecture, models, evidence and reproducibility
├── data/README.md                    # Expected external-data layout and sources
├── requirements/                     # Optional environment groupings
└── scripts/                          # Report and integrity-manifest builders
```

## Where to inspect the evidence

| Question | Canonical location |
|---|---|
| What was submitted? | [`submission/`](submission/) |
| How was the submission scored? | [`artifacts/results/audits/`](artifacts/results/audits/) |
| How did all eight temporal models behave on Track 1? | [`track1_temporal_transfer_20260906/`](artifacts/results/track1_temporal_transfer_20260906/) |
| What supports the temporal ClinVar claims? | [`clinvar_temporal_20260905/`](artifacts/results/clinvar_temporal_20260905/) and [`wide_rank_temporal_20260906/`](artifacts/results/wide_rank_temporal_20260906/) |
| What supports the patient/HPO claims? | [`patient_generalization_protocol_20260907/`](artifacts/results/patient_generalization_protocol_20260907/) and the [`five-seed study`](artifacts/results/phenopacket_hpo2023_generalization_repeated_20260906/) |
| What supports the full-background claims? | [`varpb_2023_opencravat_gnomad/`](artifacts/results/varpb_2023_opencravat_gnomad/), [`varpb_2023_external_predictors/`](artifacts/results/varpb_2023_external_predictors/), and [`varpb_augmentation_20260908/`](artifacts/results/varpb_augmentation_20260908/) |
| Are those files intact? | [`generalization_artifacts.json`](artifacts/manifests/generalization_artifacts.json) and `rdkid verify` |
| How does the code fit together? | [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) |
| What does each model mean? | [`docs/MODEL_CATALOG.md`](docs/MODEL_CATALOG.md) |
| What can and cannot be concluded? | [`docs/GENERALIZATION.md`](docs/GENERALIZATION.md) |
| What is the complete methods record? | [33-page report](artifacts/reports/bonamame_track1_report.pdf) |

## Data, privacy, and compute

No patient FASTQ, complete patient VCF, clinical document, access token, or local annotation database is committed. The gated challenge bundle is approximately 85 GB compressed and must remain on approved controlled storage. Public datasets are also kept outside version control because the VarPB background and OpenCRAVAT databases are large; their exact versions and checksums are documented.

Set `RDKID_DATA_DIR` to keep external inputs outside the checkout:

```bash
export RDKID_DATA_DIR=/secure/path/rare_disease_data
rdkid doctor --profile wgs
```

See [`data/README.md`](data/README.md) and [`docs/DATA_AND_PRIVACY.md`](docs/DATA_AND_PRIVACY.md) for the expected layout and privacy boundary.

## Scope and future work

The current results establish strong challenge performance, temporal transfer, patient-level reranking, and competitive full-background shortlist recall. Extending that evidence beyond the evaluated setting suggests several concrete next steps:

- Evaluate the complete candidate-generation and ranking workflow prospectively on multi-patient hospital cohorts with native HPO profiles, pedigrees, and full VCFs.
- Confirm phase with long reads or family sequencing and pursue functional validation of the `BUB1B` p.Asn1002Lys allele.
- Test posterior calibration and selective referral on independent clinical shifts, ancestries, sequencing platforms, and disease groups.
- Run the registered training-only ablations—original features, REVEL only, and AlphaMissense only—to separate the contribution of the pairwise objective from each imported predictor.
- Measure complete pipeline cost, including external annotation resources, alongside the compactness of the Bayesian integration model.
- Validate the strong Top-100 result on another untouched full-background benchmark and study whether it converts into higher end-to-end diagnostic yield.

The repository is intended for research and competition prioritization. Clinical deployment would require these broader validations and the corresponding governance, quality-control, and regulatory processes.

## Citation, authorship, and license

This is joint work by **Mame Diarra Toure** and **Bonaventure F. P. Dossou**. Citation metadata are provided in [`CITATION.cff`](CITATION.cff), and the upstream genome-processing contribution is cited in the report.

No software license has been added because the repository owners have not yet selected one.
