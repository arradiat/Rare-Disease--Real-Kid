#!/usr/bin/env python3
"""Re-score every saved Track 1 prediction with the public challenge rules."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from rare_disease_real_kid.paths import PROJECT_ROOT, RESULTS, SUBMISSION


ROOT = PROJECT_ROOT
TRUE_VARIANTS = frozenset(
    {
        ("chr15", 40209701, "T", "G"),
        ("chr15", 40220612, "T", "G"),
    }
)
PREDICTIONS = {
    "selected_wide_rank50_sbnn": SUBMISSION / "model1_bonamame_wide-rank50-sbnn.csv",
    "linear_bayesian_strict": RESULTS / "bayesian_run_20260905/strict_bayesian_candidates.csv",
    "linear_bayesian_broad": RESULTS / "bayesian_run_20260905/broad_bayesian_candidates.csv",
    "temporal_gradient_boost": RESULTS / "track1_temporal_transfer_20260906/strict_gradient_boost.csv",
    "temporal_dense_mlp": RESULTS / "track1_temporal_transfer_20260906/strict_dense_mlp.csv",
    "temporal_deep_ensemble": RESULTS / "track1_temporal_transfer_20260906/strict_deep_ensemble.csv",
    "temporal_small_sbnn": RESULTS / "track1_temporal_transfer_20260906/strict_small_sbnn.csv",
    "temporal_dense_bbb": RESULTS / "track1_temporal_transfer_20260906/strict_dense_bbb.csv",
    "temporal_wide_rank50_sbnn": RESULTS / "track1_temporal_transfer_20260906/strict_wide_rank50_sbnn.csv",
    "temporal_wide_dense_bbb": RESULTS / "track1_temporal_transfer_20260906/strict_wide_dense_bbb.csv",
    "temporal_consensus": RESULTS / "track1_temporal_transfer_20260906/strict_temporal_consensus.csv",
    "augmented_pairwise_sbnn": RESULTS / "track1_augmented_sbnn_20260908/strict_augmented_pairwise_sbnn.csv",
}


def variant(row: dict[str, str], suffix: str) -> tuple[str, int, str, str]:
    return (
        row[f"chrom_{suffix}"].strip(),
        int(row[f"pos_{suffix}"]),
        row[f"ref_{suffix}"].strip().upper(),
        row[f"alt_{suffix}"].strip().upper(),
    )


def load_rows(path: Path) -> list[tuple[frozenset[tuple[str, int, str, str]], float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        parsed = []
        for row in csv.DictReader(handle):
            variants = {variant(row, "1")}
            if row.get("pos_2", "").strip():
                variants.add(variant(row, "2"))
            parsed.append((frozenset(variants), float(row["epcr"])))
    if len(parsed) > 10:
        raise ValueError(f"{path} has {len(parsed)} rows; the challenge maximum is 10")
    return sorted(parsed, key=lambda item: item[1], reverse=True)


def rank_points(rank: int | None) -> float:
    if rank is None:
        return 0.0
    for last_rank, points in ((1, 100), (3, 50), (5, 25), (10, 10)):
        if rank <= last_rank:
            return float(points)
    return 0.0


def score(path: Path) -> dict[str, object]:
    rows = load_rows(path)
    full_match_rank = next(
        (index for index, (variants, _) in enumerate(rows, start=1) if variants == TRUE_VARIANTS),
        None,
    )

    best_f = 0.0
    best_threshold = None
    best_rows = 0
    for threshold in sorted({epcr for _, epcr in rows}, reverse=True):
        selected = [variants for variants, epcr in rows if epcr >= threshold]
        predicted = set().union(*selected)
        tp = len(predicted & TRUE_VARIANTS)
        fp = len(predicted - TRUE_VARIANTS)
        fn = len(TRUE_VARIANTS - predicted)
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f_score = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        if f_score > best_f:
            best_f = f_score
            best_threshold = threshold
            best_rows = len(selected)

    return {
        "path": str(path.relative_to(ROOT)),
        "rows": len(rows),
        "full_match_rank": full_match_rank,
        "rank_points": rank_points(full_match_rank),
        "f_max": best_f,
        "f_max_threshold": best_threshold,
        "rows_at_f_max": best_rows,
        "exact_perfect_score": rank_points(full_match_rank) == 100.0 and best_f == 1.0,
    }


def main() -> None:
    results = {name: score(path) for name, path in PREDICTIONS.items() if path.exists()}
    missing = {name: str(path.relative_to(ROOT)) for name, path in PREDICTIONS.items() if not path.exists()}
    payload = {
        "scoring_rule": "Public SageBio Track 1 evaluation.py logic, checked 2026-09-06",
        "scoring_source": (
            "https://huggingface.co/spaces/SageBio/"
            "rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py"
        ),
        "scope_note": (
            "Local exact re-scoring against the accepted pair. The selected wide rank-50 "
            "SBNN CSV is confirmed by the live Hugging Face evaluator; the remaining model "
            "adapters are local official-equivalent scores."
        ),
        "leaderboard_confirmation": {
            "source": "Hugging Face submission receipt supplied by the team",
            "team": "bonamame",
            "account": "mametoure",
            "submission_number": 2,
            "proband_id": "PROBAND01",
            "full_match_rank": 1,
            "rank_points": 100.0,
            "f_max": 1.0,
            "f_max_threshold": 0.3724070924857989,
        },
        "accepted_pair": [list(item) for item in sorted(TRUE_VARIANTS)],
        "results": results,
        "missing_predictions": missing,
    }
    output = RESULTS / "audits/exact_track1_score_audit.json"
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
