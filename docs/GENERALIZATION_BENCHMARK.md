# Patient-level generalization benchmark

The generalization claim is split into two questions because no unrestricted
public dataset supplies real clinical phenotypes, confirmed answers and complete
patient genomes in one package.

1. **Phenopacket Store:** tests whether molecular scores and HPO matching recover
   causal hypotheses in 719 post-2023 solved cases. It includes 64 two-variant
   cases, 382 cases whose causal gene is absent from the older phenotype
   reference, and 390 cases whose disease identifier is absent.
2. **VarPB:** tests whether a causal missense variant can be found inside the
   exonic background of 107 openly redistributable Personal Genome Project
   samples. This is a full-background variant-ranking test, but it does not
   include patient-specific HPO phenotypes and is not a complete clinical cohort.

These experiments are complementary. Neither should be described as a direct
estimate of diagnostic sensitivity in a hospital.

## Executed Phenopacket protocol

The original HPO mixing coefficient is selected only on older cases. It is then
frozen for the untouched post-2023 cases and for three phenotype-shift tests:
20% term removal, 40% term removal, and addition of 20% unrelated terms.

```bash
python phenopacket_hpo_generalization.py \
  --historical data/open/clinvar_temporal/clinvar_20231226.vcf.gz \
  --current data/open/clinvar_temporal/clinvar_20260829.vcf.gz \
  --phenopackets data/open/phenopacket_store/all_phenopackets_v0.1.27.zip \
  --hpo data/open/phenopacket_store/hp_2023-10-09.obo \
  --model-checkpoint results/track1_temporal_transfer_20260906/trained_models.pkl \
  --output-dir results/patient_generalization_protocol_20260907

python patient_generalization_benchmark.py \
  --cases results/patient_generalization_protocol_20260907/patient_cases.csv \
  --predictions results/patient_generalization_protocol_20260907/patient_predictions.csv.gz \
  --output-dir results/patient_generalization_protocol_20260907/evaluation
```

The second command recomputes candidate recall, end-to-end Top-1/5/10, MRR,
patient-bootstrap confidence intervals, score calibration, unseen-gene OOD,
selective accuracy and inheritance/disease/gene strata from case-level rows.
If candidate generation drops the causal hypothesis, its reciprocal rank is zero
and every end-to-end Top-k metric counts it as a failure.

The current compact benchmark deliberately inserts the truth into ten candidates,
so its candidate-generation recall is 1 by construction. That number is a schema
check, not evidence that the genome-wide filter has 100% sensitivity.

## Executed single-seed robustness check

For candidate seed 20260905, the temporal consensus reaches 59.81% Top-1 with
the original HPO terms, 58.83% after dropping 20%, 58.55% after dropping 40%,
and 59.39% after adding 20% unrelated terms. The wide rank-50 SBNN reaches
57.44%, 57.16%, 56.33%, and 56.88%, respectively. These small changes support
robustness to moderate phenotype noise on this compact candidate benchmark.

The harder result is gene novelty. The temporal consensus reaches 33.77% Top-1
on 382 unseen-gene cases, compared with 89.32% on 337 seen-gene cases. On the
390 disease-unseen cases it reaches 34.36%. This gap is the main unsolved
generalization problem and must remain visible in the report.

## Executed VarPB full-background protocol

VarPB version 3 is archived at <https://doi.org/10.5281/zenodo.10291114>.
The complete archive is 624 MB; the only benchmark input needed for this run is
`allSamplesExonicVariants.txt.gz` (96.7 MB), plus the small causal-variant tables
for 2023. The run used all 17,872,513 rows from 107 backgrounds and all 2,577
2023 causal variants. It did not use the precomputed result matrices.

The exact executed preparation and annotation commands were:

```bash
python prepare_varpb_opencravat.py prepare-vcf \
  --background data/open/varpb/allSamplesExonicVariants.txt.gz \
  --causal data/open/varpb/2023_autosomalDominant_notInClinvarAtAllBefore.txt \
           data/open/varpb/2023_autosomalRecessive_notInClinvarAtAllBefore.txt \
  --output data/open/varpb/varpb_2023_union_hg19.vcf.gz \
  --index data/open/varpb/varpb_2023_union.sqlite

.venv-opencravat/bin/oc run data/open/varpb/varpb_2023_union_hg19.vcf.gz \
  -i vcf -l hg19 -a gnomad revel alphamissense \
  --mp 6 --endat aggregator \
  -n varpb_2023_gnomad -d data/open/varpb/opencravat

python prepare_varpb_opencravat.py export-annotations \
  --database data/open/varpb/opencravat/varpb_2023_gnomad.sqlite \
  --fallback-vcf data/open/varpb/varpb_2023_union_hg19.vcf.gz \
  --output data/open/varpb/varpb_2023_opencravat_gnomad.tsv.gz
```

The exported annotation table contains exactly 894,327 input keys. OpenCRAVAT
provided mapped annotations for 893,268; 1,059 mapper failures were retained as
explicit missing values rather than silently removed. Only two of the 2,577
causal keys needed this fallback. This gives 100% causal-key coverage.

All eight saved models were then scored in bounded-memory chunks:

```text
Chromosome  Position  Ref  Alt  Gene  Consequence  PopulationAF
```

Then score every model in bounded-memory chunks:

```bash
python score_varpb_models.py \
  --annotated-variants data/open/varpb/varpb_2023_opencravat_gnomad.tsv.gz \
  --checkpoint results/track1_temporal_transfer_20260906/trained_models.pkl \
  --output-dir results/varpb_2023_opencravat_gnomad

python evaluate_varpb_suite.py \
  --background data/open/varpb/allSamplesExonicVariants.txt.gz \
  --dominant data/open/varpb/2023_autosomalDominant_notInClinvarAtAllBefore.txt \
  --recessive data/open/varpb/2023_autosomalRecessive_notInClinvarAtAllBefore.txt \
  --scores-dir results/varpb_2023_opencravat_gnomad \
  --output-dir results/varpb_2023_opencravat_gnomad/evaluation
```

The evaluator performed 2,577 causal variants by 107 backgrounds, or 275,739
spike-in trials per method. Scores cover every background row. Gene-history and
ClinVar-review inputs were neutralized because they are unavailable for generic
background variants, preventing label leakage.

## Full-background results

| Method | Median rank | MRR | Top-20 | Top-100 |
|---|---:|---:|---:|---:|
| Original small SBNN | 522 | 0.00213 | 0.000% | 0.200% |
| Wide rank-50 SBNN | 1,562 | 0.00113 | 0.000% | 0.078% |
| Temporal consensus | 2,071 | 0.00088 | 0.000% | 0.039% |
| Dense BBB | 2,966 | 0.00053 | 0.000% | 0.039% |
| Deep ensemble | 4,416 | 0.00028 | 0.000% | 0.000% |
| Gradient boosting | 4,252 | 0.00029 | 0.000% | 0.000% |
| Predictor-augmented pairwise SBNN | 14 | 0.15958 | 65.699% | **86.632%** |
| REVEL, external comparator | **8.5** | **0.25157** | **65.897%** | **84.823%** |
| AlphaMissense, external comparator | 29 | 0.16317 | 44.235% | 69.723% |

REVEL and AlphaMissense were evaluated as external pretrained comparators, not
as models trained by this project. Missing predictor values were assigned zero
and all 894,327 union variants remained in the score files, so their results are
not inflated by dropping unannotated variants.

The first run was a useful negative result: real annotation improved the saved
models substantially over the declared-missense ablation, but their coarse 2023
ClinVar feature contract did not transfer to full-background exome ranking. The
post-submission extension then retrained a rank-50 SBNN on pre-2023 ClinVar only,
using frozen REVEL and AlphaMissense inputs and a pairwise ranking loss. It
essentially matches REVEL at Top-20 and exceeds it at Top-100, but remains weaker
at median rank, Top-1, and MRR. See
[`POST_SUBMISSION_GENERALIZATION_2026-09-08.md`](POST_SUBMISSION_GENERALIZATION_2026-09-08.md)
for the exact chronology, feature additions, 24-cell joint screen, and locked
evaluation. Because the external predictor values are inputs, the result supports
a Bayesian evidence-integration claim rather than independent replacement of
REVEL or AlphaMissense.

## Failure criteria fixed before evaluation

- A method does not generalize end to end if causal-candidate recall is poor,
  even when conditional reranking looks strong.
- Seen-gene performance cannot substitute for unseen-gene performance.
- A Bayesian model does not provide useful OOD detection merely because it has a
  posterior; OOD AUROC and selective accuracy must improve empirically.
- Moderate phenotype corruption should not cause a large Top-k collapse.
- VarPB results apply to missense spike-ins, not structural, deep-intronic,
  repeat-expansion or complete compound-heterozygous diagnosis.

## Sources

- [PhEval framework and workflow](https://monarch-initiative.github.io/pheval/)
- [Phenopacket Store](https://github.com/monarch-initiative/phenopacket-store)
- [VarPB repository](https://github.com/ZuchnerLab/VariantPrioritizationBenchmark)
- [VarPB version 3 data](https://doi.org/10.5281/zenodo.10291114)
