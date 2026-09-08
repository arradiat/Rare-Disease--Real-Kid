#!/usr/bin/env python3
"""Train a predictor-augmented SBNN without reading VarPB labels.

Training and model selection use only ClinVar 2023. The final command can score
the locked VarPB union after training, but no VarPB causal table is accepted by
this program. This separation makes answer leakage harder by construction.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import pickle
import sqlite3
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import StandardScaler

from clinvar_temporal_experiment import Adam, BinarySBNN, Variant, featurize, gene_statistics, parse_release, ranking_metrics
from score_varpb_models import as_variant, chunks, official_chromosome, open_text


PREDICTOR_COLUMNS = ("revel__score", "alphamissense__am_pathogenicity")
AUGMENTED_FEATURES = ("revel_score", "revel_missing", "alphamissense_score", "alphamissense_missing")


def clean_key(chrom, pos, ref, alt):
    return str(chrom).removeprefix("chr"), int(pos), str(ref).upper(), str(alt).upper()


def predictor_lookup(database: Path) -> dict[tuple[str, int, str, str], tuple[float | None, float | None]]:
    query = '''
        SELECT original_input__chrom, extra_vcf_info__pos, extra_vcf_info__ref,
               extra_vcf_info__alt, revel__score, alphamissense__am_pathogenicity
        FROM variant
    '''
    result = {}
    with sqlite3.connect(database) as connection:
        for chrom, pos, ref, alt, revel, alpha in connection.execute(query):
            if pos is None or ref is None or alt is None:
                continue
            key = clean_key(chrom, pos, ref, alt)
            old = result.get(key, (None, None))
            values = []
            for previous, current in zip(old, (revel, alpha)):
                available = [float(value) for value in (previous, current) if value is not None]
                values.append(max(available) if available else None)
            result[key] = tuple(values)
    return result


def augmented_features(records: list[Variant], stats, predictors) -> np.ndarray:
    base = featurize(records, stats)
    extra = np.empty((len(records), 4), dtype=np.float32)
    for index, record in enumerate(records):
        revel, alpha = predictors.get(
            clean_key(record.chrom, record.pos, record.ref, record.alt), (None, None)
        )
        extra[index] = (
            0.0 if revel is None else revel,
            float(revel is None),
            0.0 if alpha is None else alpha,
            float(alpha is None),
        )
    return np.concatenate([base, extra], axis=1)


@dataclass(frozen=True)
class Config:
    ranks: tuple[int, int, int]
    kl_weight: float
    posterior_scale: float
    prior_scale: float

    @property
    def name(self):
        return f"rank{'-'.join(map(str, self.ranks))}_kl{self.kl_weight:g}_post{self.posterior_scale:g}_prior{self.prior_scale:g}"


def screen_configs():
    return [
        Config(rank, kl, posterior, prior)
        for rank in ((16, 16, 1), (32, 32, 1), (50, 50, 1))
        for kl in (0.001, 0.01)
        for posterior in (0.05, 0.1)
        for prior in (0.3, 1.0)
    ]


def matched_parameter_counts(input_dim: int = 54, hidden=(200, 100), ranks=(50, 50, 1)) -> dict:
    """Count fitted scalars for architecture-matched uncertainty methods."""
    sbnn = BinarySBNN(input_dim, hidden=hidden, ranks=ranks, seed=0)
    sbnn_variational = int(sum(value.size for _, value in sbnn.parameters()))
    widths = (input_dim, *hidden, 1)
    deterministic_single = int(sum(widths[i] * widths[i + 1] + widths[i + 1] for i in range(len(widths) - 1)))
    dense_bbb_variational = 2 * deterministic_single
    deep_ensemble_five = 5 * deterministic_single
    return {
        "augmented_sbnn_variational_scalars": sbnn_variational,
        "matched_dense_bbb_variational_scalars": dense_bbb_variational,
        "matched_single_deterministic_network_scalars": deterministic_single,
        "matched_five_network_ensemble_scalars": deep_ensemble_five,
        "sbnn_reduction_vs_dense_bbb_percent": 100.0 * (1.0 - sbnn_variational / dense_bbb_variational),
        "sbnn_reduction_vs_five_network_ensemble_percent": 100.0 * (1.0 - sbnn_variational / deep_ensemble_five),
        "adam_training_scalars_sbnn": 3 * sbnn_variational,
        "adam_training_scalars_dense_bbb": 3 * dense_bbb_variational,
        "counting_note": "Model plus two Adam state arrays; excludes activations and fixed external annotation resources.",
    }


def fit_pairwise(
    x, y, x_val, y_val, config: Config, *, seed: int, epochs: int,
    pairs_per_epoch: int, batch_size: int = 512,
):
    model = BinarySBNN(
        x.shape[1], hidden=(200, 100), ranks=config.ranks, seed=seed,
        posterior_scale_multiplier=config.posterior_scale,
        prior_sigma1=config.prior_scale,
    )
    optimizer = Adam(model.parameters(), learning_rate=8e-4)
    rng = np.random.default_rng(seed)
    positives = np.flatnonzero(y == 1)
    negatives = np.flatnonzero(y == 0)
    validation_positive = np.flatnonzero(y_val == 1)
    validation_negative = np.flatnonzero(y_val == 0)
    best = None
    best_loss = math.inf
    history = []
    stale = 0
    for epoch in range(epochs):
        p = rng.choice(positives, pairs_per_epoch, replace=pairs_per_epoch > len(positives))
        n = rng.choice(negatives, pairs_per_epoch, replace=pairs_per_epoch > len(negatives))
        order = rng.permutation(pairs_per_epoch)
        for start in range(0, pairs_per_epoch, batch_size):
            batch = order[start:start + batch_size]
            warmup = min(1.0, (epoch + 1) / 8.0)
            _, _, _, gradients = model.pairwise_loss_and_gradients(
                x[p[batch]], x[n[batch]], rng,
                warmup * config.kl_weight / pairs_per_epoch,
            )
            optimizer.step(gradients)
        pair_count = min(5000, len(validation_positive), len(validation_negative))
        vp = rng.choice(validation_positive, pair_count, replace=False)
        vn = rng.choice(validation_negative, pair_count, replace=False)
        joined = np.concatenate([x_val[vp], x_val[vn]])
        probability = model.predict_statistics(joined, draws=12, seed=seed + epoch + 1)[0]
        margin = np.log(np.clip(probability[:pair_count], 1e-7, 1 - 1e-7) / np.clip(1 - probability[:pair_count], 1e-7, 1))
        margin -= np.log(np.clip(probability[pair_count:], 1e-7, 1 - 1e-7) / np.clip(1 - probability[pair_count:], 1e-7, 1))
        value = float(np.mean(np.logaddexp(0.0, -margin)))
        history.append({"epoch": epoch + 1, "validation_pairwise_nll": value})
        if value < best_loss - 1e-4:
            best_loss = value
            best = {name: parameter.copy() for name, parameter in model.parameters()}
            stale = 0
        else:
            stale += 1
            if stale >= 6:
                break
    if best is not None:
        for name, parameter in model.parameters():
            parameter[...] = best[name]
    return model, history


def validation_metrics(model, x, y, seed):
    probability, uncertainty, standard_deviation = model.predict_statistics(x, draws=32, seed=seed)
    ranking = ranking_metrics(y, probability, seed=seed, negatives_per_case=99, max_cases=2000)
    return {
        "ranking": ranking,
        "mean_mutual_information": float(np.mean(uncertainty)),
        "mean_predictive_std": float(np.mean(standard_deviation)),
    }


def train(args):
    records, _ = parse_release(args.historical, collect_all_ids=False)
    records = [
        record for record in records
        if len(record.ref) == len(record.alt) == 1 and "missense_variant" in record.consequence
    ]
    labels = np.asarray([record.label for record in records])
    groups = np.asarray([record.gene for record in records])
    train_index, validation_index = next(GroupShuffleSplit(
        n_splits=1, test_size=0.2, random_state=args.seed
    ).split(np.zeros(len(records)), labels, groups))
    stats = gene_statistics(records, train_index)
    train_records = [records[int(index)] for index in train_index]
    validation_records = [records[int(index)] for index in validation_index]
    predictors = predictor_lookup(args.predictor_database)
    train_raw = augmented_features(train_records, stats, predictors)
    validation_raw = augmented_features(validation_records, stats, predictors)
    scaler = StandardScaler().fit(train_raw)
    x_train = scaler.transform(train_raw).astype(np.float32)
    x_validation = scaler.transform(validation_raw).astype(np.float32)
    y_train = labels[train_index]
    y_validation = labels[validation_index]
    rng = np.random.default_rng(args.seed)
    screen_index = np.concatenate([
        rng.choice(np.flatnonzero(y_train == label), min(args.screen_per_class, int(np.sum(y_train == label))), replace=False)
        for label in (0, 1)
    ])
    validation_screen = np.concatenate([
        rng.choice(np.flatnonzero(y_validation == label), min(7500, int(np.sum(y_validation == label))), replace=False)
        for label in (0, 1)
    ])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    screened = []
    for number, config in enumerate(screen_configs(), 1):
        model, history = fit_pairwise(
            x_train[screen_index], y_train[screen_index],
            x_validation[validation_screen], y_validation[validation_screen],
            config, seed=args.seed + number, epochs=args.screen_epochs,
            pairs_per_epoch=args.screen_pairs,
        )
        metrics = validation_metrics(model, x_validation[validation_screen], y_validation[validation_screen], args.seed + 1000 + number)
        item = {"config": {**asdict(config), "name": config.name}, "metrics": metrics, "history": history}
        screened.append(item)
        (args.output_dir / "screen.json").write_text(json.dumps(screened, indent=2) + "\n", encoding="utf-8")
        print(f"{number:02d}/{len(screen_configs())} {config.name}: MRR={metrics['ranking']['mrr']:.4f}", flush=True)
    winner_row = max(screened, key=lambda item: (item["metrics"]["ranking"]["mrr"], item["metrics"]["ranking"]["top10"]))
    values = winner_row["config"]
    winner = Config(tuple(values["ranks"]), values["kl_weight"], values["posterior_scale"], values["prior_scale"])
    model, history = fit_pairwise(
        x_train, y_train, x_validation, y_validation, winner,
        seed=args.seed, epochs=args.final_epochs, pairs_per_epoch=args.final_pairs,
    )
    metrics = validation_metrics(model, x_validation, y_validation, args.seed + 5000)
    artifact = {
        "model": model,
        "scaler": scaler,
        "gene_statistics": stats,
        "predictor_columns": PREDICTOR_COLUMNS,
        "augmented_features": AUGMENTED_FEATURES,
        "config": asdict(winner),
        "seed": args.seed,
        "parameter_counts": matched_parameter_counts(train_raw.shape[1], (200, 100), winner.ranks),
    }
    with (args.output_dir / "augmented_sbnn.pkl").open("wb") as handle:
        pickle.dump(artifact, handle)
    report = {
        "experiment": "post-submission predictor-augmented pairwise SBNN",
        "historical_source": str(args.historical),
        "varpb_used_for_training_or_selection": False,
        "split": "gene-disjoint 80/20 GroupShuffleSplit",
        "training_records": len(train_records),
        "validation_records": len(validation_records),
        "training_positive": int(y_train.sum()),
        "validation_positive": int(y_validation.sum()),
        "input_features": train_raw.shape[1],
        "parameter_counts": matched_parameter_counts(train_raw.shape[1], (200, 100), winner.ranks),
        "winner": asdict(winner),
        "validation": metrics,
        "history": history,
    }
    (args.output_dir / "training_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


def score(args):
    with args.checkpoint.open("rb") as handle:
        artifact = pickle.load(handle)
    predictors = predictor_lookup(args.predictor_database)
    output = args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with open_text(args.annotated_variants, "rt") as incoming, output.open("w", newline="", encoding="utf-8") as outgoing:
        reader = csv.DictReader(incoming, delimiter="\t")
        writer = csv.writer(outgoing, delimiter="\t")
        writer.writerow(["Chromosome", "Position", "Ref", "Alt", "bonamame_augmented_pairwise_sbnn_score"])
        for block in chunks(reader, args.chunk_size):
            records = [as_variant(row, total + index + 1) for index, row in enumerate(block)]
            raw = augmented_features(records, {}, predictors)
            x = artifact["scaler"].transform(raw).astype(np.float32)
            probability = artifact["model"].predict_statistics(x, draws=args.draws, seed=args.seed + total)[0]
            for row, value in zip(block, probability):
                writer.writerow([official_chromosome(row["Chromosome"]), row["Position"], row["Ref"], row["Alt"], float(value)])
            total += len(block)
            print(f"Scored {total:,} variants", flush=True)
    print({"variants": total, "output": str(output), "draws": args.draws})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    training = subparsers.add_parser("train")
    training.add_argument("--historical", type=Path, required=True)
    training.add_argument("--predictor-database", type=Path, required=True)
    training.add_argument("--output-dir", type=Path, required=True)
    training.add_argument("--seed", type=int, default=20260908)
    training.add_argument("--screen-per-class", type=int, default=5000)
    training.add_argument("--screen-epochs", type=int, default=8)
    training.add_argument("--screen-pairs", type=int, default=8000)
    training.add_argument("--final-epochs", type=int, default=30)
    training.add_argument("--final-pairs", type=int, default=40000)
    scoring = subparsers.add_parser("score-varpb")
    scoring.add_argument("--checkpoint", type=Path, required=True)
    scoring.add_argument("--annotated-variants", type=Path, required=True)
    scoring.add_argument("--predictor-database", type=Path, required=True)
    scoring.add_argument("--output", type=Path, required=True)
    scoring.add_argument("--draws", type=int, default=32)
    scoring.add_argument("--chunk-size", type=int, default=25000)
    scoring.add_argument("--seed", type=int, default=20260908)
    args = parser.parse_args()
    train(args) if args.command == "train" else score(args)


if __name__ == "__main__":
    main()
