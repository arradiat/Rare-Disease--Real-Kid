# Generalization Evidence

The project deliberately separates success on the solved challenge child from transfer to later variants, independent solved cases, and realistic genome-scale backgrounds.

| Question | Dataset and protocol | Committed evidence | Interpretation |
|---|---|---|---|
| Does the model rank the accepted challenge pair? | Strict and broad child candidate tables | [Eight-model transfer audit](../artifacts/results/track1_temporal_transfer_20260906/track1_temporal_model_audit.json) and [augmented-model audit](../artifacts/results/track1_augmented_sbnn_20260908/track1_augmented_sbnn_audit.json) | Directly verified shortlist ranking |
| Does training transfer forward in time? | ClinVar 2023 training versus later unseen Variation IDs | [Temporal comparison](../artifacts/results/clinvar_temporal_20260905/consensus_holdout.json), [dense BBB](../artifacts/results/full_bbb_temporal_20260906/full_bbb_results.json), and [wide-rank study](../artifacts/results/wide_rank_temporal_20260906/wide_rank_results.json) | Temporal variant ranking, calibration, and OOD comparison |
| Does phenotype conditioning help across solved cases? | 719 post-2023 Phenopackets with resampled candidate sets | [Five-seed summary](../artifacts/results/phenopacket_hpo2023_generalization_repeated_20260906/repeated_seed_summary.json), [cases](../artifacts/results/patient_generalization_protocol_20260907/patient_cases.csv), [predictions](../artifacts/results/patient_generalization_protocol_20260907/patient_predictions.csv.gz), and [patient evaluation](../artifacts/results/patient_generalization_protocol_20260907/evaluation/patient_generalization_results.json) | Patient-oriented conditional reranking, robustness, calibration, and novelty strata |
| Does ranking transfer to large genome backgrounds? | 2,577 VarPB causal variants crossed with 107 backgrounds | [Submission-era models](../artifacts/results/varpb_2023_opencravat_gnomad/evaluation/summary.json), [external predictors](../artifacts/results/varpb_2023_external_predictors/evaluation/summary.json), and [augmented SBNN](../artifacts/results/varpb_augmentation_20260908/evaluation/summary.json) | Full-background shortlist recall on 275,739 spike-in trials |

## Evidence included in the repository

- Exact strict and broad Track 1 predictions for all eight submission-era models.
- Temporal ClinVar metrics, calibration outputs, hyperparameter studies, and model-comparison figures.
- All 719 patient cases, compressed per-candidate predictions, bootstrap evaluations, phenotype-corruption conditions, and five repeated candidate-set seeds.
- Full VarPB score exports for every submission-era model, REVEL, AlphaMissense, and the predictor-augmented SBNN.
- Per-trial compressed rank tables, summary CSV/JSON files, comparison figures, training filter counts, the 24-cell model-selection screen, and the augmented training report.
- The fitted augmented SBNN checkpoint and its frozen Track 1 transfer audit.
- SHA-256 and size records for the complete evidence set in [`generalization_artifacts.json`](../artifacts/manifests/generalization_artifacts.json).

## Main findings

The wide rank-50 SBNN was the strongest temporal ranker, reaching 84.35% Top-1 among 100 candidates. The deep ensemble was the strongest evaluated unseen-gene detector, so the best ranker and the best OOD discriminator were not the same model.

Across five resampled Phenopacket candidate cohorts, HPO conditioning improved Top-1 ranking for all compared models. Dense BBB averaged 58.89% Top-1 across all 719 cases; wide dense BBB averaged 68.13% on the 64 two-variant cases. These values measure ranking within constructed ten-candidate sets.

The frozen submission-era models performed poorly on VarPB's full backgrounds: the best, small SBNN, reached only 0.20% Top-100 with median rank 522. The predictor-augmented pairwise SBNN raised Top-100 to 86.63% and median rank to 14, while REVEL remained stronger at Top-1 and MRR. This negative-to-positive sequence demonstrates why compact shortlist performance cannot be treated as genome-wide generalization.

## Leakage controls

- Temporal evaluation excludes Variation IDs found in the training release.
- Patient candidate generation is repeated under fixed seeds, and gene/disease novelty is evaluated explicitly.
- The augmented SBNN uses a gene-disjoint train/validation split.
- VarPB causal tables are not accepted by the augmented model's training command.
- The Track 1 transfer excludes the answer, prior rank/score, disease name, ClinVar significance, and direct challenge predictor columns.
- Both alleles in a pair use the same posterior weight draw before their probabilities are combined.

## Limits of the evidence

Constructed candidate sets measure conditional ranking, not end-to-end diagnostic sensitivity. VarPB is a spike-in benchmark rather than a prospective hospital cohort. REVEL and AlphaMissense are inputs to the augmented SBNN, so its result supports learned evidence integration rather than independent replacement of those predictors. The two child alleles are not directly proven to be in trans, and the software is not a clinically validated diagnostic device.
