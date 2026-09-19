#!/usr/bin/env python3
"""Evaluate every VarPB model export and build one comparison figure."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt

from rare_disease_real_kid.evaluation.varpb import evaluate, open_text


def score_column(path: Path) -> str:
    with open_text(path, "rt") as handle:
        columns = next(csv.reader(handle, delimiter="\t"))
    if len(columns) != 5 or not columns[-1].endswith("_score"):
        raise ValueError(f"Unexpected VarPB score header in {path}: {columns}")
    return columns[-1]


def run_suite(background: Path, dominant: Path, recessive: Path, scores_dir: Path, output_dir: Path) -> list[dict]:
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for path in sorted(scores_dir.glob("bonamame_*_varpb.tsv")):
        column = score_column(path)
        print(f"Evaluating {column}", flush=True)
        result, records = evaluate(background, dominant, recessive, path, column)
        summaries.append(result)
        stem = column.removeprefix("bonamame_").removesuffix("_score")
        with gzip.open(output_dir / f"{stem}_ranks.tsv.gz", "wt", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
    (output_dir / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")
    fields = ["model", "top1", "top5", "top10", "top20", "top100", "mrr", "median_rank", "mean_rank", "auc_top20"]
    with (output_dir / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in summaries:
            writer.writerow({"model": result["model"], **{name: result["overall"][name] for name in fields[1:]}})

    ordered = sorted(summaries, key=lambda item: item["overall"]["mrr"], reverse=True)
    names = [item["model"].removeprefix("bonamame_").removesuffix("_score").replace("_", "\n") for item in ordered]
    x = range(len(ordered))
    figure, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    axes[0].bar(x, [item["overall"]["top20"] for item in ordered], color="#355C7D")
    axes[0].set_ylabel("Fraction of causal spike-ins in top 20")
    axes[0].set_ylim(0, 1)
    axes[1].bar(x, [item["overall"]["mrr"] for item in ordered], color="#C06C84")
    axes[1].set_ylabel("Mean reciprocal rank")
    for axis in axes:
        axis.set_xticks(list(x), names, fontsize=8)
        axis.grid(axis="y", alpha=0.2)
    figure.suptitle("VarPB 2023 full-background declared-missense ablation")
    figure.tight_layout()
    figure.savefig(output_dir / "varpb_2023_model_comparison.png", dpi=220, bbox_inches="tight")
    plt.close(figure)
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--dominant", type=Path, required=True)
    parser.add_argument("--recessive", type=Path, required=True)
    parser.add_argument("--scores-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summaries = run_suite(args.background, args.dominant, args.recessive, args.scores_dir, args.output_dir)
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
