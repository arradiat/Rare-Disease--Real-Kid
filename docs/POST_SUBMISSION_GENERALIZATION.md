# Post-submission generalization extension: 8 September 2026

## Boundary and purpose

The Track 1 report and the live Hugging Face submission were complete on 7
September 2026. At that point the child result was stable, every retained model
ranked the accepted `BUB1B` pair first, and the selected wide rank-50 SBNN had
received 100 Rank Points and F-max 1.000. The report explicitly named a
full-background VarPB evaluation as remaining work. Everything below was added
after that stopping point and must be described as a post-submission extension.

The extension tests a harder claim: can a model trained without the challenge
child or VarPB causal answers prioritize later pathogenic variants inside real
107-genome exonic backgrounds? The failure criterion was fixed in advance: a
perfect child score or strong ten-candidate result would not compensate for poor
full-background rank.

## Step 1: acquire the smallest sufficient public benchmark

Only three files from VarPB v3 were downloaded from
<https://doi.org/10.5281/zenodo.10291114>: the 96.7 MB background table and the
2023 dominant and recessive causal tables. The 624 MB convenience archive and
precomputed method matrices were unnecessary. Checksums, exact filenames, and
commands are in `data/README.md` and `docs/REPRODUCIBILITY.md`.

The background contains 17,872,513 sample-variant rows from 107 Personal Genome
Project genomes. `prepare_varpb_opencravat.py` used a disk-backed SQLite key to
deduplicate these to 891,753 background variants, then added 2,574 unique causal
keys, producing a union of 894,327 hg19 variants. This avoided loading the table
into memory and did not require the 85 GB challenge bundle.

## Step 2: run local annotation and audit key preservation

OpenCRAVAT 3.1.1 was installed in `.venv-opencravat`. The local modules were hg38
mapper 1.47.2, gnomAD 2.2.0/data v2.1, REVEL 2024.07.30/data v4.3a, and
AlphaMissense 1.1.1. Their installed footprint is about 28 GB. CADD was not
installed because its OpenCRAVAT module alone is about 86.6 GB.

OpenCRAVAT mapped 893,268 union variants. The export recovered original hg19 VCF
keys rather than OpenCRAVAT's internally trimmed indel representation. It kept
1,059 mapper failures as explicit missing values, including two causal variants,
giving 894,327 output rows and 100% causal-key coverage. Dropping these rows
would have biased the benchmark.

## Step 3: stress-test the eight submission-era models

All eight saved models were scored with ClinVar-review and gene-history features
neutralized. Each method was evaluated on 2,577 causal variants crossed with 107
backgrounds: 275,739 trials. The best submission-era model was the original
small SBNN, but its median causal rank was 522, MRR 0.00213, and Top-100 0.20%.
The wide rank-50 SBNN had median rank 1,562. This was a negative result: models
trained for ClinVar classification and compact candidate reranking did not
transfer to whole-background prioritization.

## Step 4: establish external predictor targets

REVEL and AlphaMissense were exported from the same OpenCRAVAT database. Missing
scores were set to zero and every union key stayed in the output. They are
external pretrained comparators, not models trained by this project.

| Method | Median rank | Top-1 | Top-20 | Top-100 | MRR |
|---|---:|---:|---:|---:|---:|
| REVEL | **8.5** | **12.48%** | **65.90%** | 84.82% | **0.2516** |
| AlphaMissense | 29 | 8.84% | 44.24% | 69.72% | 0.1632 |

This comparison identified the missing information: the original 50-feature
contract did not contain modern protein-level predictor signals.

## Step 5: build a leakage-controlled augmented training set

The exact NCBI ClinVar GRCh38 archive from 26 December 2023 was downloaded and
matched the previously recorded SHA-256. A deterministic filter retained only
unambiguous pathogenic/likely-pathogenic or benign/likely-benign, biallelic SNVs
whose `MC` field contains `missense_variant`:

- 2,340,600 input ClinVar records;
- 1,084,930 unambiguous labeled records;
- 146,810 labeled missense SNVs;
- 54,548 pathogenic and 92,262 benign examples.

The same REVEL and AlphaMissense versions annotated 146,786 records. Twenty-four
mapper failures were encoded by a zero value plus a separate missingness flag.
The augmented input has the original 50 features plus REVEL value/missingness and
AlphaMissense value/missingness, for 54 total features.

The first implementation attempt sent the complete 2,340,600-record archive to
OpenCRAVAT. It was stopped during mapping because unlabeled records and long
indels cannot contribute REVEL or AlphaMissense scores and produced extensive
mapper-error output. This was an engineering failure, not an evaluated model.
The replacement filter was written as a tested, deterministic script before the
annotation was rerun. Its retained-row counts are stored in
`clinvar_training_filter.json`; the incomplete OpenCRAVAT directory is ignored
and is not used by any experiment.

VarPB causal tables were not read by the training program. Genes, rather than
individual rows, were split 80/20, yielding 117,700 training and 29,110
validation variants. This prevents variants from the same gene appearing on
both sides of model selection.

The exact training-data commands were:

```bash
curl -L -k \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2023/clinvar_20231226.vcf.gz \
  -o data/open/clinvar_temporal/clinvar_20231226.vcf.gz

python -m rare_disease_real_kid.experiments.prepare_clinvar_predictors \
  --source data/open/clinvar_temporal/clinvar_20231226.vcf.gz \
  --output data/open/clinvar_temporal/clinvar_20231226_labeled_missense.vcf.gz \
  --manifest artifacts/results/varpb_augmentation_20260908/clinvar_training_filter.json

.venv-opencravat/bin/oc run \
  data/open/clinvar_temporal/clinvar_20231226_labeled_missense.vcf.gz \
  -i vcf -l hg38 -a revel alphamissense --mp 6 --endat aggregator \
  -n clinvar_20231226_labeled_missense_predictors \
  -d data/open/clinvar_temporal/opencravat_filtered --user-email-opt-out
```

## Step 6: align the SBNN objective with ranking

The earlier SBNN minimized binary cross-entropy. The extension adds a Bayesian
RankNet loss:

```text
loss = log(1 + exp(-(score(pathogenic) - score(benign)))) + scaled KL
```

Both variants in a comparison use the same posterior weight draw. Sampling two
different networks would confound model uncertainty with the desired ordering.
The KL is scaled once by the number of training comparisons and warmed up over
the first eight epochs.

Model selection used a complete 3 x 2 x 2 x 2 factorial screen:

- ranks: `(16,16,1)`, `(32,32,1)`, `(50,50,1)`;
- KL weights: `0.001`, `0.01`;
- posterior initialization scales: `0.05`, `0.1`;
- prior scales: `0.3`, `1.0`.

All 24 configurations used the same split, architecture widths 200 and 100,
screening budget, and metrics. The selected configuration was rank `(50,50,1)`,
KL `0.001`, posterior scale `0.1`, and prior scale `1.0`. After retraining on the
full training partition, gene-disjoint validation reached Top-1 62.15%, Top-10
96.15%, MRR 0.7449, and mean rank 2.684 among 100 candidates.

The executed training command was:

```bash
python -m rare_disease_real_kid.experiments.augmented_varpb train \
  --historical data/open/clinvar_temporal/clinvar_20231226.vcf.gz \
  --predictor-database data/open/clinvar_temporal/opencravat_filtered/clinvar_20231226_labeled_missense_predictors.sqlite \
  --output-dir artifacts/results/varpb_augmentation_20260908 \
  --screen-per-class 3000 --screen-epochs 5 --screen-pairs 4000 \
  --final-epochs 20 --final-pairs 30000 --seed 20260908
```

## Step 7: unlock VarPB once

The selected model produced 32 posterior draws for all 894,327 union variants.
The locked VarPB result is:

```bash
python -m rare_disease_real_kid.experiments.augmented_varpb score-varpb \
  --checkpoint artifacts/checkpoints/augmented_pairwise_sbnn/augmented_sbnn.pkl \
  --annotated-variants data/open/varpb/varpb_2023_opencravat_gnomad.tsv.gz \
  --predictor-database data/open/varpb/opencravat/varpb_2023_gnomad.sqlite \
  --output artifacts/results/varpb_augmentation_20260908/varpb_scores/bonamame_augmented_pairwise_sbnn_varpb.tsv \
  --draws 32 --chunk-size 25000 --seed 20260908

python -m rare_disease_real_kid.evaluation.varpb_suite \
  --background data/open/varpb/allSamplesExonicVariants.txt.gz \
  --dominant data/open/varpb/2023_autosomalDominant_notInClinvarAtAllBefore.txt \
  --recessive data/open/varpb/2023_autosomalRecessive_notInClinvarAtAllBefore.txt \
  --scores-dir artifacts/results/varpb_augmentation_20260908/varpb_scores \
  --output-dir artifacts/results/varpb_augmentation_20260908/evaluation
```

| Method | Median rank | Top-1 | Top-20 | Top-100 | MRR |
|---|---:|---:|---:|---:|---:|
| Augmented pairwise SBNN | 14 | 7.26% | 65.70% | **86.63%** | 0.1596 |
| REVEL | **8.5** | **12.48%** | **65.90%** | 84.82% | **0.2516** |
| AlphaMissense | 29 | 8.84% | 44.24% | 69.72% | 0.1632 |

The model is now competitive for shortlist recall: its Top-20 is within 0.20
percentage points of REVEL and its Top-100 is 1.81 points higher. It substantially
improves over AlphaMissense at Top-20 and Top-100. It does not beat REVEL at the
very top of the ranking, and it is slightly below AlphaMissense in Top-1 and MRR.
Those limitations must remain visible.

The result does not show that the SBNN independently discovered REVEL's biology:
REVEL and AlphaMissense are input features. The added value tested here is
whether a scalable Bayesian ranker can combine those predictors with sequence,
frequency, consequence, and missingness evidence while retaining competitive
full-background recall and a posterior predictive distribution.

## Frozen transfer to Track 1

The selected augmented checkpoint was also applied to the challenge shortlist
without retraining or model selection. For each posterior draw, both alleles of
a compound-heterozygous hypothesis were evaluated with the same sampled SBNN,
then combined as `phenotype_score * sqrt(p1 * p2)`. Across 128 draws, the
accepted BUB1B pair ranked first in both the strict and broad candidate sets.
The strict six-row CSV achieved 100 Rank Points and F-max 1.000 under the local
reproduction of the official evaluator. Its EPCR threshold was 0.4332283540 and
the posterior standard deviation of the top pair score was 0.00788. This result
is locally official-equivalent and has not been submitted to the live platform.

## Model size and built-in uncertainty

The parameter comparison must use the same 54-200-100-1 architecture. The
rank-50 SBNN stores a mean and scale for each low-rank factor and bias, for
56,204 fitted variational scalars. A width-matched dense BBB stores a mean and
scale for every dense weight and bias, for 62,402 scalars. Five width-matched
deterministic networks store 156,005 weights in total.

| Architecture-matched method | Stored fitted scalars | Difference from SBNN |
|---|---:|---:|
| Rank-50 SBNN | 56,204 | reference |
| Dense BBB | 62,402 | SBNN uses 9.9% fewer |
| Five-network ensemble | 156,005 | SBNN uses 64.0% fewer |

With Adam, model parameters plus the two optimizer-state arrays require 168,612
scalars for the SBNN and 187,206 for dense BBB, excluding activations. Ensemble
members can be trained sequentially to reduce peak optimizer memory, but all
five point estimates are retained for ensemble inference. The SBNN checkpoint,
including its scaler and gene statistics, is 670,999 bytes.

The asymptotic distinction is more important than the modest 9.9% saving in
this small network. For a square width-`d` layer, dense BBB grows as `2*d*d`,
whereas an SBNN of fixed rank `r` grows as approximately `4*r*d`. The relative
advantage therefore increases with width when `r` remains well below `d`.

REVEL and AlphaMissense provide one deterministic pathogenicity score per
variant. Their normal outputs do not provide posterior draws, predictive
standard deviation, or epistemic mutual information. The SBNN produces a
posterior predictive distribution from one fitted model and can report all
three. This is a capability advantage, not yet an empirical uncertainty win:
the augmented model's gene-disjoint validation mean mutual information is
0.000118 and mean predictive standard deviation is 0.00337, and VarPB does not
define an OOD label for validating those quantities. A new independent shift
benchmark is required before claiming calibrated or superior OOD detection.

The 56,204-scalar comparison covers only the Bayesian integration model. Since
REVEL and AlphaMissense are its inputs, the complete pipeline still inherits
their annotation cost. The accurate claim is that the added uncertainty-aware
ranker is compact and rank-scalable, not that the entire system is smaller than
running its upstream predictors.

## Reproducible artifacts

- `prepare_clinvar_predictor_training.py`: deterministic pre-2023 filter.
- `augmented_sbnn_varpb_experiment.py`: joint sweep, pairwise training, and VarPB inference.
- `evaluate_augmented_sbnn_track1.py`: frozen augmented-model transfer to the
  challenge candidate pairs using shared posterior draws for both alleles.
- `artifacts/results/varpb_augmentation_20260908/clinvar_training_filter.json`: filter counts.
- `artifacts/results/varpb_augmentation_20260908/screen.json`: all 24 validation results.
- `artifacts/results/varpb_augmentation_20260908/training_report.json`: selected configuration and final validation metrics.
- `artifacts/results/varpb_augmentation_20260908/evaluation/summary.json`: locked VarPB metrics.
- `artifacts/checkpoints/augmented_pairwise_sbnn/augmented_sbnn.pkl`: fitted predictor-augmented SBNN checkpoint.
- `artifacts/results/track1_augmented_sbnn_20260908/track1_augmented_sbnn_audit.json`:
  exact Track 1 transfer configuration and locally reproduced official score.
- `artifacts/figures/report/varpb_full_background.png`: comparison rebuilt from JSON.

## Next registered experiment

The next step is not to tune against VarPB. It is to freeze this result and test
three training-only ablations: original 50 features, REVEL only, and
AlphaMissense only, all with the same pairwise objective and budget. This will
separate the contribution of the ranking objective from the two imported
predictors. A later independent temporal release or another untouched benchmark
must be used for any further model choice.
