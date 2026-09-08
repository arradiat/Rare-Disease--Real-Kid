#!/usr/bin/env python3
"""Evaluate VarPB scores by exact spike-in rank without loading 17.9M rows.

This reproduces the rank definition in VarPB's ``benchmarkNewTool.py`` for the
2023 dominant and recessive causal tables.  For each of 107 background genomes,
the rank of every causal variant is computed against that person's non-reference
variants.  Sorted score arrays make this equivalent calculation tractable.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(mode.removesuffix("t"), newline="", encoding="utf-8")


def key(chrom: str, pos: str | int, ref: str, alt: str) -> str:
    return f"{str(chrom).removeprefix('chr')}:{int(pos)}:{ref.upper()}:{alt.upper()}"


def load_scores(path: Path, score_column: str) -> dict[str, float]:
    result = {}
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Chromosome", "Position", "Ref", "Alt", score_column}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")
        for row in reader:
            result[key(row["Chromosome"], row["Position"], row["Ref"], row["Alt"])] = float(row[score_column])
    return result


def load_causal(paths: list[tuple[str, Path]], scores: dict[str, float]) -> list[dict]:
    variants = []
    for inheritance, path in paths:
        with open_text(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                variant_key = key(row["Chromosome"], row["Position"], row["ReferenceAllele"], row["AlternateAllele"])
                if variant_key in scores:
                    variants.append({"inheritance": inheritance, "key": variant_key, "score": scores[variant_key]})
    return variants


def sample_score_arrays(background: Path, scores: dict[str, float]) -> tuple[dict[str, np.ndarray], dict]:
    by_sample: dict[str, list[float]] = {}
    rows = non_reference = scored = 0
    with open_text(background, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            if str(row["genotype"]) == "0":
                continue
            non_reference += 1
            value = scores.get(key(row["hg19_chr"], row["hg19_pos(1-based)"], row["ref"], row["alt"]))
            if value is None:
                continue
            scored += 1
            by_sample.setdefault(row["sampleID"], []).append(value)
    arrays = {sample: np.sort(np.asarray(values, dtype=np.float64)) for sample, values in by_sample.items()}
    return arrays, {
        "background_rows": rows,
        "non_reference_rows": non_reference,
        "scored_non_reference_rows": scored,
        "background_score_coverage": scored / non_reference if non_reference else 0.0,
        "samples": len(arrays),
    }


def rank_against(sorted_scores: np.ndarray, score: float) -> float:
    left = int(np.searchsorted(sorted_scores, score, side="left"))
    right = int(np.searchsorted(sorted_scores, score, side="right"))
    greater = len(sorted_scores) - right
    equal = right - left
    return 1.0 + greater + 0.5 * equal


def auc_to_cutoff(ranks: np.ndarray, cutoff: int) -> float:
    x = np.arange(0, cutoff + 1, dtype=float)
    y = np.asarray([np.mean(ranks <= value) for value in x])
    return float(np.trapz(y, x / cutoff))


def summarize_ranks(ranks: np.ndarray, normalized: np.ndarray) -> dict:
    return {
        "spike_in_trials": int(len(ranks)),
        "top1": float(np.mean(ranks <= 1)),
        "top5": float(np.mean(ranks <= 5)),
        "top10": float(np.mean(ranks <= 10)),
        "top20": float(np.mean(ranks <= 20)),
        "top100": float(np.mean(ranks <= 100)),
        "mrr": float(np.mean(1.0 / ranks)),
        "median_rank": float(np.median(ranks)),
        "mean_rank": float(np.mean(ranks)),
        "auc_top20": auc_to_cutoff(ranks, 20),
        "fraction_within_normalized_0_003": float(np.mean(normalized <= 0.003)),
    }


def evaluate(background: Path, dominant: Path, recessive: Path, scores_path: Path, score_column: str) -> tuple[dict, list[dict]]:
    scores = load_scores(scores_path, score_column)
    causal = load_causal([("dominant", dominant), ("recessive", recessive)], scores)
    samples, coverage = sample_score_arrays(background, scores)
    records = []
    for sample, background_scores in samples.items():
        for variant in causal:
            rank = rank_against(background_scores, variant["score"])
            records.append({
                "sample_id": sample,
                "inheritance": variant["inheritance"],
                "causal_variant": variant["key"],
                "rank": rank,
                "normalized_rank": rank / len(background_scores),
            })
    ranks = np.asarray([row["rank"] for row in records])
    normalized = np.asarray([row["normalized_rank"] for row in records])
    result = {
        "model": score_column,
        "score_rows": len(scores),
        "causal_variants_in_score_file": len(causal),
        "coverage": coverage,
        "overall": summarize_ranks(ranks, normalized),
    }
    for inheritance in ("dominant", "recessive"):
        mask = np.asarray([row["inheritance"] == inheritance for row in records])
        result[inheritance] = summarize_ranks(ranks[mask], normalized[mask])
    return result, records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--dominant", type=Path, required=True)
    parser.add_argument("--recessive", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--score-column", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result, records = evaluate(args.background, args.dominant, args.recessive, args.scores, args.score_column)
    stem = args.score_column.removesuffix("_score")
    (args.output_dir / f"{stem}_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with gzip.open(args.output_dir / f"{stem}_ranks.tsv.gz", "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
