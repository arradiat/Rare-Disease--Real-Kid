# Model Catalog

The registry in [`registry.py`](../src/rare_disease_real_kid/models/registry.py) is the authoritative list. Every registered module is imported in the automated tests.

| Registered name | Family | Architecture or rule | Role |
|---|---|---|---|
| `gradient_boost` | Deterministic tree | 220 boosting iterations, at most 31 leaves | Tabular baseline |
| `dense_mlp` | Deterministic neural network | 50-32-16-1 | Single-network baseline |
| `deep_ensemble` | Ensemble | Five 50-32-16-1 MLPs | Disagreement/OOD baseline |
| `small_sbnn` | Low-rank Bayesian network | 50-32-16-1, ranks 8-8-1 | Original SBNN baseline |
| `dense_bbb` | Dense Bayesian network | 50-32-16-1 | Full-weight BBB comparison |
| `wide_rank50_sbnn` | Low-rank Bayesian network | 50-200-100-1, ranks 50-50-1 | Selected submission model |
| `wide_dense_bbb` | Dense Bayesian network | 50-200-100-1 | Width-matched BBB comparison |
| `temporal_consensus` | Probability mixture | 0.3 tree + 0.2 ensemble + 0.5 small SBNN | Complementary ranking model |
| `augmented_pairwise_sbnn` | Predictor-augmented SBNN | 54-200-100-1, ranks 50-50-1 | Full-background evidence integration |

## Selected Track 1 model

The selected wide rank-50 SBNN fixes:

- input width: 50;
- hidden widths: 200 and 100;
- ranks: 50, 50, and 1;
- KL weight: 0.01;
- posterior-scale multiplier: 0.1;
- prior scale: 0.3;
- seed: 20260905;
- Track 1 posterior draws: 128.

The committed transfer audit is [`track1_temporal_model_audit.json`](../artifacts/results/track1_temporal_transfer_20260906/track1_temporal_model_audit.json). The canonical submission is independently checked by [`track1_scores.py`](../src/rare_disease_real_kid/evaluation/track1_scores.py).

## Interpreting the comparisons

The best model depends on the question. The wide rank-50 SBNN was selected for temporal variant-ranking performance; the deep ensemble showed stronger unseen-gene OOD discrimination. The augmented pairwise SBNN adds external predictor evidence and a ranking loss, so it tests evidence integration rather than whether an SBNN independently replaces REVEL or AlphaMissense.

For the matched 54-200-100-1 comparison, the rank-50 SBNN has 56,204 variational scalars, dense BBB has 62,402, and five deterministic networks retain 156,005 weights. The counts, validation protocol, winning hyperparameters, and training history are committed in [`training_report.json`](../artifacts/results/varpb_augmentation_20260908/training_report.json).

## Native uncertainty outputs

- Deterministic tree/MLP models produce no epistemic statistic.
- The ensemble exposes member disagreement.
- SBNN and BBB models expose posterior predictive samples, predictive deviation, and mutual information.
- The existence of these outputs is not evidence that they are calibrated; calibration and OOD performance are evaluated separately.
