#!/usr/bin/env python3
"""Evaluate the frozen predictor-augmented SBNN on Track 1 candidate pairs.

The checkpoint was trained on ClinVar 2023 and selected without the challenge
case. REVEL and AlphaMissense values are read from the same local OpenCRAVAT
modules used for the external generalization experiment.
"""

from __future__ import annotations

import argparse
import json
import pickle
import re
import sqlite3
from pathlib import Path

import numpy as np

from rare_disease_real_kid.checkpoint import load_checkpoint

from rare_disease_real_kid.experiments.augmented_varpb import augmented_features
from rare_disease_real_kid.models.bayesian_layers import _sigmoid
from rare_disease_real_kid.evaluation.track1_transfer import (
    SEED,
    TIERS,
    challenge_variant,
    read_pairs,
    write_submission,
)
from rare_disease_real_kid.evaluation.track1_scores import score as official_score


def parse_module_score(value) -> float | None:
    """Return the maximum numeric score from an OpenCRAVAT module cell."""
    if value is None:
        return None
    parsed = []
    for item in re.split(r"[;,|]", str(value)):
        try:
            number = float(item)
        except ValueError:
            continue
        if np.isfinite(number):
            parsed.append(number)
    return max(parsed) if parsed else None


def module_scores(database: Path, column: str, records) -> dict:
    """Look up GRCh38 predictor values for the requested challenge alleles."""
    result = {}
    with sqlite3.connect(database) as connection:
        for record in records:
            key = (record.chrom, record.pos, record.ref, record.alt)
            if key in result:
                continue
            if not re.fullmatch(r"(?:[0-9]+|X|Y|M|MT)", record.chrom):
                raise ValueError(f"Unsupported chromosome name: {record.chrom!r}")
            table = "chr" + record.chrom
            try:
                rows = connection.execute(
                    f'SELECT "{column}" FROM "{table}" '
                    "WHERE pos=? AND ref=? AND alt=?",
                    (record.pos, record.ref, record.alt),
                ).fetchall()
            except sqlite3.OperationalError:
                rows = []
            values = [parse_module_score(row[0]) for row in rows]
            available = [value for value in values if value is not None]
            result[key] = max(available) if available else None
    return result


def posterior_pair_statistics(model, x, phenotype, *, draws: int, seed: int):
    """Propagate shared posterior draws through both alleles of every pair."""
    if x.shape[0] != 2 * len(phenotype):
        raise ValueError("Expected exactly two allele rows per candidate pair")
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(draws):
        allele_probability = _sigmoid(model.forward(x, rng)[0]).reshape(-1, 2)
        pair_probability = np.sqrt(
            np.clip(allele_probability[:, 0], 1e-30, 1.0)
            * np.clip(allele_probability[:, 1], 1e-30, 1.0)
        )
        samples.append(pair_probability * phenotype)
    samples = np.asarray(samples)
    return samples.mean(axis=0), samples.std(axis=0)


def evaluate(args) -> dict:
    artifact = load_checkpoint(args.checkpoint, allow_external=True)

    audit = {}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for tier, member in TIERS.items():
        rows = read_pairs(args.ranking_zip, member)
        records = [
            challenge_variant(row, allele, 2 * index + allele)
            for index, row in enumerate(rows)
            for allele in (1, 2)
        ]
        revel = module_scores(args.revel_database, "score", records)
        alphamissense = module_scores(
            args.alphamissense_database, "am_pathogenicity", records
        )
        predictors = {
            key: (revel.get(key), alphamissense.get(key))
            for key in set(revel) | set(alphamissense)
        }
        raw = augmented_features(records, artifact["gene_statistics"], predictors)
        x = artifact["scaler"].transform(raw).astype(np.float32)
        phenotype = np.asarray([float(row["phenotype_score"]) for row in rows])
        pair_mean, pair_std = posterior_pair_statistics(
            artifact["model"], x, phenotype, draws=args.draws, seed=args.seed
        )
        output = args.output_dir / f"{tier}_augmented_pairwise_sbnn.csv"
        write_submission(output, rows, pair_mean)
        result = official_score(output.resolve())
        order = np.argsort(-pair_mean, kind="stable")
        result.update(
            top_gene=rows[int(order[0])]["gene"],
            top_pair_score=float(pair_mean[order[0]]),
            top_pair_posterior_std=float(pair_std[order[0]]),
            candidate_pairs=len(rows),
            unique_alleles=len(predictors),
            revel_alleles=int(sum(value is not None for value in revel.values())),
            alphamissense_alleles=int(
                sum(value is not None for value in alphamissense.values())
            ),
        )
        audit[tier] = result
        print(
            tier,
            result["top_gene"],
            result["rank_points"],
            result["f_max"],
            flush=True,
        )

    report = {
        "experiment": "frozen predictor-augmented SBNN transferred to Track 1",
        "checkpoint": str(args.checkpoint),
        "challenge_used_for_training_or_model_selection": False,
        "architecture": "54-200-100-1 SBNN",
        "ranks": list(artifact["config"]["ranks"]),
        "pair_rule": (
            "posterior mean of phenotype_score * "
            "sqrt(p_allele_1 * p_allele_2), with a shared weight draw"
        ),
        "draws": args.draws,
        "seed": args.seed,
        "results": audit,
    }
    (args.output_dir / "track1_augmented_sbnn_audit.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranking-zip", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--revel-database", type=Path, required=True)
    parser.add_argument("--alphamissense-database", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=128)
    parser.add_argument("--seed", type=int, default=SEED + 3)
    args = parser.parse_args()
    if args.draws < 2:
        parser.error("--draws must be at least 2")
    print(json.dumps(evaluate(args), indent=2))


if __name__ == "__main__":
    main()
