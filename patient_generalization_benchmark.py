#!/usr/bin/env python3
"""Evaluate candidate generation, ranking, calibration and OOD by patient.

The evaluator deliberately consumes plain CSV files rather than one particular
pipeline's objects.  This lets the same metrics be used for Phenopacket/PhEval
experiments, VarPB spike-ins and future clinical cohorts without changing the
metric definitions.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_auc_score


CASE_FIELDS = {
    "case_id",
    "condition",
    "causal_hypothesis_id",
    "causal_gene",
    "disease_id",
    "inheritance",
    "gene_seen",
    "disease_seen",
    "raw_variant_count",
}
PREDICTION_FIELDS = {
    "case_id",
    "condition",
    "model",
    "hypothesis_id",
    "score",
    "uncertainty",
    "gene",
    "runtime_seconds",
}


def _read_csv(path: Path, required: set[str]) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return list(reader)


def _boolean(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no"}:
        return False
    raise ValueError(f"Expected a boolean value, got {value!r}")


def _number(value: str, *, nonnegative: bool = False) -> float:
    result = float(value)
    if not math.isfinite(result) or (nonnegative and result < 0):
        raise ValueError(f"Invalid numeric value: {value!r}")
    return result


def _optional_nonnegative(value: str) -> float:
    if value.strip() == "":
        return float("nan")
    return _number(value, nonnegative=True)


def load_cases(path: Path) -> list[dict]:
    cases = []
    seen = set()
    for row in _read_csv(path, CASE_FIELDS):
        key = (row["case_id"], row["condition"])
        if key in seen:
            raise ValueError(f"Duplicate case/condition row: {key}")
        seen.add(key)
        raw_count = int(row["raw_variant_count"])
        if raw_count < 1:
            raise ValueError(f"raw_variant_count must be positive for {key}")
        cases.append(
            {
                **row,
                "key": key,
                "gene_seen": _boolean(row["gene_seen"]),
                "disease_seen": _boolean(row["disease_seen"]),
                "raw_variant_count": raw_count,
            }
        )
    if not cases:
        raise ValueError("The case manifest is empty")
    return cases


def load_predictions(path: Path, valid_case_keys: set[tuple[str, str]]) -> dict:
    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    hypotheses = set()
    for row in _read_csv(path, PREDICTION_FIELDS):
        case_key = (row["case_id"], row["condition"])
        if case_key not in valid_case_keys:
            raise ValueError(f"Prediction refers to an unknown case/condition: {case_key}")
        key = (*case_key, row["model"])
        hypothesis_key = (*key, row["hypothesis_id"])
        if hypothesis_key in hypotheses:
            raise ValueError(f"Duplicate predicted hypothesis: {hypothesis_key}")
        hypotheses.add(hypothesis_key)
        score = _number(row["score"])
        uncertainty = _optional_nonnegative(row["uncertainty"])
        runtime = _optional_nonnegative(row["runtime_seconds"])
        grouped[key].append(
            {**row, "score": score, "uncertainty": uncertainty, "runtime_seconds": runtime}
        )
    if not grouped:
        raise ValueError("The predictions file is empty")
    return grouped


def average_rank(rows: list[dict], causal_hypothesis_id: str) -> float | None:
    truth = next((row for row in rows if row["hypothesis_id"] == causal_hypothesis_id), None)
    if truth is None:
        return None
    greater = sum(row["score"] > truth["score"] for row in rows)
    tied_other = sum(row["score"] == truth["score"] for row in rows) - 1
    return 1.0 + greater + 0.5 * tied_other


def normalized_probabilities(rows: list[dict]) -> np.ndarray:
    """Turn non-negative EPCR-like scores into a per-case ranking distribution."""
    scores = np.asarray([row["score"] for row in rows], dtype=float)
    minimum = float(scores.min())
    if minimum < 0:
        scores = scores - minimum
    total = float(scores.sum())
    if total <= 0:
        return np.full(len(scores), 1.0 / len(scores))
    return scores / total


def bootstrap_interval(values: np.ndarray, seed: int, draws: int) -> list[float]:
    if len(values) == 0:
        return [float("nan"), float("nan")]
    rng = np.random.default_rng(seed)
    samples = rng.integers(0, len(values), size=(draws, len(values)))
    means = values[samples].mean(axis=1)
    return [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]


def selective_metrics(correct: np.ndarray, uncertainty: np.ndarray) -> dict[str, float]:
    finite = np.isfinite(uncertainty)
    if not np.any(finite):
        return {"top1_at_80pct_coverage": float("nan"), "risk_coverage_auc": float("nan")}
    correct = correct[finite]
    uncertainty = uncertainty[finite]
    order = np.argsort(uncertainty, kind="stable")
    retained_correct = correct[order]
    coverage = np.arange(1, len(order) + 1, dtype=float) / len(order)
    risk = 1.0 - np.cumsum(retained_correct) / np.arange(1, len(order) + 1)
    keep = max(1, int(math.ceil(0.8 * len(order))))
    return {
        "top1_at_80pct_coverage": float(retained_correct[:keep].mean()),
        "risk_coverage_auc": float(np.trapz(risk, coverage)),
    }


def summarize(cases: list[dict], grouped: dict, model: str, condition: str, seed: int, draws: int) -> dict:
    selected = [case for case in cases if case["condition"] == condition]
    ranks = []
    generated = []
    top_uncertainty = []
    gene_unseen = []
    nll = []
    brier = []
    candidate_counts = []
    reduction = []
    runtimes = []
    inheritance = defaultdict(list)

    for case in selected:
        rows = grouped.get((*case["key"], model), [])
        rank = average_rank(rows, case["causal_hypothesis_id"]) if rows else None
        present = rank is not None
        generated.append(float(present))
        ranks.append(float(rank) if present else float("inf"))
        candidate_counts.append(len(rows))
        reduction.append(max(0.0, 1.0 - len(rows) / case["raw_variant_count"]))
        inheritance[case["inheritance"]].append(len(ranks) - 1)
        if rows:
            top = max(rows, key=lambda row: row["score"])
            top_uncertainty.append(top["uncertainty"])
            gene_unseen.append(float(not case["gene_seen"]))
            runtimes.append(max(row["runtime_seconds"] for row in rows))
        else:
            top_uncertainty.append(float("nan"))
            gene_unseen.append(float(not case["gene_seen"]))
            runtimes.append(float("nan"))
        if present:
            probabilities = normalized_probabilities(rows)
            truth_index = next(i for i, row in enumerate(rows) if row["hypothesis_id"] == case["causal_hypothesis_id"])
            truth_probability = float(np.clip(probabilities[truth_index], 1e-12, 1.0))
            targets = np.zeros(len(rows), dtype=float)
            targets[truth_index] = 1.0
            nll.append(-math.log(truth_probability))
            brier.append(float(np.mean((probabilities - targets) ** 2)))

    ranks_array = np.asarray(ranks)
    generated_array = np.asarray(generated)
    reciprocal = np.where(np.isfinite(ranks_array), 1.0 / ranks_array, 0.0)
    correct = (ranks_array <= 1).astype(float)
    uncertainty_array = np.asarray(top_uncertainty)
    unseen_array = np.asarray(gene_unseen)
    finite_ood = np.isfinite(uncertainty_array)
    ood = float("nan")
    if len(np.unique(unseen_array[finite_ood])) == 2:
        ood = float(roc_auc_score(unseen_array[finite_ood], uncertainty_array[finite_ood]))

    def stratum(indices: list[int]) -> dict[str, float | int]:
        values = ranks_array[indices]
        return {
            "cases": len(indices),
            "candidate_recall": float(np.mean(np.isfinite(values))),
            "top1_end_to_end": float(np.mean(values <= 1)),
            "top5_end_to_end": float(np.mean(values <= 5)),
            "top10_end_to_end": float(np.mean(values <= 10)),
            "mrr_end_to_end": float(np.mean(np.where(np.isfinite(values), 1.0 / values, 0.0))),
        }

    finite_runtime = np.asarray(runtimes)[np.isfinite(runtimes)]
    report = {
        "cases": len(selected),
        "candidate_generation_recall": float(generated_array.mean()),
        "top1_end_to_end": float(correct.mean()),
        "top5_end_to_end": float(np.mean(ranks_array <= 5)),
        "top10_end_to_end": float(np.mean(ranks_array <= 10)),
        "mrr_end_to_end": float(reciprocal.mean()),
        "top1_conditional_on_truth_retained": float(correct[generated_array.astype(bool)].mean()) if generated_array.any() else float("nan"),
        "top1_bootstrap_95": bootstrap_interval(correct, seed, draws),
        "mrr_bootstrap_95": bootstrap_interval(reciprocal, seed + 1, draws),
        "candidate_count_median": float(np.median(candidate_counts)),
        "candidate_reduction_rate_mean": float(np.mean(reduction)),
        "candidate_distribution_nll": float(np.mean(nll)) if nll else float("nan"),
        "candidate_distribution_brier": float(np.mean(brier)) if brier else float("nan"),
        "unseen_gene_ood_auroc": ood,
        "runtime_seconds_median": float(np.median(finite_runtime)) if len(finite_runtime) else float("nan"),
        "runtime_seconds_p95": float(np.quantile(finite_runtime, 0.95)) if len(finite_runtime) else float("nan"),
        "by_gene_novelty": {
            "seen": stratum([i for i, case in enumerate(selected) if case["gene_seen"]]),
            "unseen": stratum([i for i, case in enumerate(selected) if not case["gene_seen"]]),
        },
        "by_disease_novelty": {
            "seen": stratum([i for i, case in enumerate(selected) if case["disease_seen"]]),
            "unseen": stratum([i for i, case in enumerate(selected) if not case["disease_seen"]]),
        },
        "by_inheritance": {name: stratum(indices) for name, indices in sorted(inheritance.items())},
    }
    report.update(selective_metrics(correct, uncertainty_array))
    return report


def evaluate(case_path: Path, prediction_path: Path, *, seed: int = 20260907, draws: int = 2000) -> dict:
    cases = load_cases(case_path)
    grouped = load_predictions(prediction_path, {case["key"] for case in cases})
    models = sorted({key[2] for key in grouped})
    conditions = sorted({case["condition"] for case in cases})
    results = {
        condition: {
            model: summarize(cases, grouped, model, condition, seed, draws)
            for model in models
        }
        for condition in conditions
    }
    return {
        "claim": "Patient-level generalization is measured without giving a ranker credit when candidate generation loses the causal hypothesis.",
        "metric_contract": {
            "candidate_generation_recall": "Fraction of cases whose exact causal hypothesis reaches the ranker.",
            "end_to_end_top_k": "Missing causal hypotheses count as failures.",
            "conditional_top1": "Ranking accuracy only among cases whose causal hypothesis survived candidate generation.",
            "calibration": "Per-case non-negative scores are normalized across retained hypotheses before NLL and Brier calculation.",
            "ood": "AUROC of top-ranked-hypothesis uncertainty for distinguishing unseen from seen causal genes.",
        },
        "bootstrap": {"unit": "patient case", "draws": draws, "seed": seed},
        "results": results,
    }


def write_summary_csv(path: Path, report: dict) -> None:
    fields = [
        "condition", "model", "cases", "candidate_generation_recall",
        "top1_end_to_end", "top5_end_to_end", "top10_end_to_end",
        "mrr_end_to_end", "top1_conditional_on_truth_retained",
        "unseen_gene_ood_auroc", "candidate_distribution_nll",
        "candidate_distribution_brier", "top1_at_80pct_coverage",
        "risk_coverage_auc", "runtime_seconds_median", "runtime_seconds_p95",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for condition, models in report["results"].items():
            for model, metrics in models.items():
                writer.writerow({"condition": condition, "model": model, **{key: metrics[key] for key in fields[2:]}})


def plot_summary(path: Path, report: dict) -> None:
    conditions = list(report["results"])
    models = sorted({model for condition in conditions for model in report["results"][condition]})
    figure, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    metrics = (
        ("candidate_generation_recall", "Causal candidate retained"),
        ("top1_end_to_end", "End-to-end Top-1"),
        ("mrr_end_to_end", "End-to-end MRR"),
    )
    x = np.arange(len(models))
    width = 0.8 / max(1, len(conditions))
    for axis, (metric, title) in zip(axes, metrics):
        for index, condition in enumerate(conditions):
            offset = (index - (len(conditions) - 1) / 2) * width
            axis.bar(x + offset, [report["results"][condition][model][metric] for model in models], width, label=condition)
        axis.set_title(title)
        axis.set_xticks(x, [model.replace("_", "\n") for model in models], rotation=20, ha="right")
        axis.grid(axis="y", alpha=0.2)
        axis.set_ylim(0, 1)
    axes[0].set_ylabel("Fraction")
    axes[-1].legend()
    figure.suptitle("Patient-level generalization: generation and ranking are scored separately")
    figure.tight_layout()
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def json_safe(value):
    """Replace non-finite metric placeholders with JSON null."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--bootstrap-draws", type=int, default=2000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = evaluate(args.cases, args.predictions, seed=args.seed, draws=args.bootstrap_draws)
    safe_report = json_safe(report)
    (args.output_dir / "patient_generalization_results.json").write_text(
        json.dumps(safe_report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    write_summary_csv(args.output_dir / "patient_generalization_summary.csv", report)
    plot_summary(args.output_dir / "patient_generalization_summary.png", report)
    print(json.dumps(safe_report["results"], indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
