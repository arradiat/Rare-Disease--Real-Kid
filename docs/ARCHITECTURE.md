# Architecture

## System boundary

The project has two connected but independently auditable layers:

1. The WGS pipeline discovers biologically plausible, inheritance-compatible candidate pairs from the gated child data.
2. The modeling layer ranks those candidate pairs and measures transfer on public temporal and patient-oriented benchmarks.

The separation prevents a learned pathogenicity score from bypassing genotype quality, frequency, consequence, inheritance, phase, or phenotype checks.

## Data flow

```text
private inputs -> pipeline/ -> candidate tables -> features/
                                              -> models/
                                              -> evaluation/
public datasets -> experiments/ -------------------^
                                                      |
                                     artifacts/results + submission/
```

## Package ownership

- `pipeline/` owns data discovery, HPO extraction, reference preparation, alignment, normalization, phasing, annotation, pair construction, and Track 1 table creation.
- `features/` owns the stable allele and pair feature contracts.
- `models/` contains only model definitions, fitting/prediction adapters, and the registry.
- `experiments/` owns temporal splits, hyperparameter studies, public-dataset preparation, and training protocols.
- `evaluation/` owns challenge scoring, submission validation, patient metrics, and VarPB scoring.
- `reporting/` reads committed artifacts and produces figures or environment manifests.

All repository paths are defined in [`paths.py`](../src/rare_disease_real_kid/paths.py). External data default to `data/` and can be redirected with `RDKID_DATA_DIR`. An installed command can target another checkout with `RDKID_PROJECT_ROOT`.

## Scientific invariants

Feature order, splits, random seeds, model widths, ranks, KL weights, posterior scales, prior scales, and pair-score formulas are explicit and covered by tests. The smoke tests in [`test_checkpoint_parity.py`](../tests/smoke/test_checkpoint_parity.py) compare the two temporal checkpoint representations under identical inputs and seeds and confirm that the augmented pairwise checkpoint loads through the trusted artifact boundary.

The canonical checkpoint uses the new package module paths. The original pickle is retained as `trained_models.legacy.pkl`, and [`checkpoint.py`](../src/rare_disease_real_kid/checkpoint.py) remaps only the two historical project module names needed to load it.

## Trust boundary

Pickle can execute code during deserialization. The default loader accepts only checkpoints inside `artifacts/`; explicit callers must opt in to an external trusted checkpoint. Patient data and authentication tokens are outside the repository boundary.
