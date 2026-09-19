#!/usr/bin/env python3
"""Rebuild the figures used in the hackathon methods report from saved results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from rare_disease_real_kid.paths import FIGURES, PROJECT_ROOT

ROOT = PROJECT_ROOT


def read_json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def save(fig: plt.Figure, output: Path, name: str) -> None:
    fig.tight_layout()
    fig.savefig(output / name, dpi=220, bbox_inches="tight")
    plt.close(fig)


def analysis_funnel(output: Path) -> None:
    scoring = read_json("artifacts/results/audits/scoring_audit.json")
    annotation = read_json("artifacts/results/audits/annotation_join_audit.json")
    stages = [
        "WGS variants",
        "Heterozygous alleles",
        "Gene annotations",
        "Same-gene pairs",
        "Strict QC pairs",
        "Final shortlist",
    ]
    counts = [
        5_079_064,
        annotation["matched_heterozygous_alleles"] + annotation["unmatched_heterozygous_alleles"],
        annotation["candidate_gene_annotations"],
        scoring["all_pair_count"],
        scoring["strict_primary_pair_count"],
        6,
    ]
    colors = ["#284B63", "#356B83", "#4D8AA0", "#78A6B5", "#D69F5D", "#B85C38"]
    fig, axis = plt.subplots(figsize=(11.5, 5.8))
    y = np.arange(len(stages))
    axis.barh(y, counts, color=colors)
    axis.set_xscale("log")
    axis.set_yticks(y, stages)
    axis.invert_yaxis()
    axis.set_xlabel("Number retained (log scale)")
    axis.set_title("From 5.1 million variants to six reportable compound-heterozygous hypotheses")
    for position, count in zip(y, counts):
        axis.text(count * 1.08, position, f"{count:,}", va="center", fontsize=10)
    axis.grid(axis="x", alpha=0.2)
    save(fig, output, "analysis_funnel.png")


def system_pipeline(output: Path) -> None:
    stages = [
        ("1  Candidate generation", "FASTQ / VCF\nalignment, calling and QC", "Bonaventure pipeline", "#284B63"),
        ("2  Clinical filtering", "frequency, consequence,\nClinVar, inheritance, phase", "Bonaventure pipeline", "#356B83"),
        ("3  Molecular ranking", "gradient boost, ensemble,\nSBNN and dense BBB", "learned comparison", "#8B5E83"),
        ("4  OOD warning", "ensemble disagreement and\nBayesian posterior uncertainty", "uncertainty layer", "#B85C38"),
        ("5  Patient reasoning", "HPO similarity, disease validity,\ncompound-heterozygous logic", "combined evidence", "#4F8A5B"),
        ("6  Track 1 decision", "ranked candidate pairs, EPCR,\nand evidence trace", "challenge output", "#7A8B3A"),
    ]
    fig, axis = plt.subplots(figsize=(15, 5.2))
    axis.set_xlim(0, len(stages))
    axis.set_ylim(0, 1)
    axis.axis("off")
    for index, (title, body, owner, color) in enumerate(stages):
        x = index + 0.08
        box = FancyBboxPatch(
            (x, 0.24), 0.84, 0.54,
            boxstyle="round,pad=0.025,rounding_size=0.035",
            linewidth=1.4,
            edgecolor=color,
            facecolor=color + "18",
        )
        axis.add_patch(box)
        axis.text(x + 0.42, 0.68, title, ha="center", va="center", fontsize=10.5, fontweight="bold", color=color)
        axis.text(x + 0.42, 0.49, body, ha="center", va="center", fontsize=9)
        axis.text(x + 0.42, 0.30, owner, ha="center", va="center", fontsize=8.5, color="#555555", style="italic")
        if index < len(stages) - 1:
            axis.annotate("", xy=(index + 1.07, 0.51), xytext=(index + 0.94, 0.51), arrowprops={"arrowstyle": "->", "lw": 1.6, "color": "#555555"})
    axis.text(3, 0.92, "The submitted system: evidence is narrowed, ranked, checked for uncertainty, and explained", ha="center", fontsize=14, fontweight="bold")
    axis.text(3, 0.10, "Deterministic validity checks remain separate from learned pathogenicity and uncertainty estimates.", ha="center", fontsize=10)
    save(fig, output, "system_pipeline.png")


def track1_scores(output: Path) -> None:
    transferred = read_json(
        "artifacts/results/track1_temporal_transfer_20260906/track1_temporal_model_audit.json"
    )["results"]["strict"]
    rows = []
    labels = {
        "gradient_boost": "Gradient boost",
        "dense_mlp": "Dense MLP",
        "deep_ensemble": "Deep ensemble",
        "small_sbnn": "Small SBNN",
        "dense_bbb": "Dense BBB",
        "wide_rank50_sbnn": "Wide rank-50 SBNN",
        "wide_dense_bbb": "Wide dense BBB",
        "temporal_consensus": "Temporal consensus",
    }
    rows.extend(
        (labels[name], item["rank_points"], item["f_max"], item["f_max_threshold"])
        for name, item in transferred.items()
    )
    augmented = read_json(
        "artifacts/results/track1_augmented_sbnn_20260908/track1_augmented_sbnn_audit.json"
    )["results"]["strict"]
    rows.append(
        (
            "Augmented pairwise SBNN",
            augmented["rank_points"],
            augmented["f_max"],
            augmented["f_max_threshold"],
        )
    )
    names = [row[0] for row in rows]
    fig, (left, right) = plt.subplots(1, 2, figsize=(13.5, 6), gridspec_kw={"width_ratios": [1.2, 1]})
    y = np.arange(len(rows))
    left.barh(y, [row[1] for row in rows], color="#3D6D8A")
    left.set_yticks(y, names)
    left.invert_yaxis()
    left.set_xlim(0, 105)
    left.set_xlabel("Rank Points")
    left.set_title("All methods recover the accepted pair at rank 1")
    for i, row in enumerate(rows):
        left.text(row[1] - 2, i, f"{row[1]:.0f}", ha="right", va="center", color="white", fontweight="bold")
    thresholds = [row[3] for row in rows]
    right.barh(y, thresholds, color="#D27A45")
    right.set_yticks(y, [""] * len(y))
    right.invert_yaxis()
    right.set_xlim(0, 1.04)
    right.set_xlabel("EPCR threshold at F-max")
    right.set_title("Perfect F-max does not imply equal calibration")
    for i, value in enumerate(thresholds):
        right.text(value + 0.015, i, f"{value:.3f}", va="center", fontsize=9)
    save(fig, output, "track1_scores.png")


def temporal_model_comparison(output: Path) -> None:
    consensus = read_json("artifacts/results/clinvar_temporal_20260905/consensus_holdout.json")["holdout"]
    full_bbb = read_json("artifacts/results/full_bbb_temporal_20260906/full_bbb_results.json")["holdout"]
    wide = read_json("artifacts/results/wide_rank_temporal_20260906/wide_rank_results.json")["holdout"]
    data = {
        "Gradient boost": (consensus["Gradient boost"]["ranking"]["top1"], consensus["Gradient boost"]["f_max"], np.nan),
        "Deep ensemble": (consensus["Deep ensemble"]["ranking"]["top1"], consensus["Deep ensemble"]["f_max"], full_bbb["deep_ensemble"]["uncertainty_unseen_gene_auroc"]),
        "Small SBNN": (consensus["SBNN"]["ranking"]["top1"], consensus["SBNN"]["f_max"], 0.6020),
        "Consensus": (consensus["Selected consensus"]["ranking"]["top1"], consensus["Selected consensus"]["f_max"], np.nan),
        "Dense BBB": (full_bbb["full_bbb"]["ranking"]["top1"], full_bbb["full_bbb"]["f_max"], full_bbb["full_bbb"]["uncertainty_unseen_gene_auroc"]),
        "Wide rank-50 SBNN": (wide["wide_low_rank_sbnn"]["ranking"]["top1"], wide["wide_low_rank_sbnn"]["f_max"], wide["wide_low_rank_sbnn"]["uncertainty_unseen_gene_auroc"]),
        "Wide dense BBB": (wide["wide_dense_bbb"]["ranking"]["top1"], wide["wide_dense_bbb"]["f_max"], wide["wide_dense_bbb"]["uncertainty_unseen_gene_auroc"]),
    }
    names = list(data)
    x = np.arange(len(names))
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    metrics = [(0, "Top-1 among 100", (0.75, 0.87)), (1, "F-max", (0.89, 0.93)), (2, "Unseen-gene OOD AUROC", (0.55, 0.78))]
    colors = ["#477998", "#477998", "#8B5E83", "#7A8B3A", "#8B5E83", "#B85C38", "#8B5E83"]
    for axis, (index, title, limits) in zip(axes, metrics):
        values = [data[name][index] for name in names]
        axis.bar(x, values, color=colors)
        axis.set_xticks(x, [name.replace(" ", "\n") for name in names], fontsize=8)
        axis.set_ylim(*limits)
        axis.set_title(title)
        axis.grid(axis="y", alpha=0.2)
        for i, value in enumerate(values):
            if np.isfinite(value):
                axis.text(i, value + (limits[1] - limits[0]) * 0.018, f"{value:.3f}", ha="center", fontsize=8)
    fig.suptitle("Untouched odd-ID temporal holdout: ranking and uncertainty answer different questions", fontsize=14)
    save(fig, output, "temporal_model_comparison.png")


def parameter_tradeoff(output: Path) -> None:
    full = read_json("artifacts/results/full_bbb_temporal_20260906/full_bbb_results.json")
    wide = read_json("artifacts/results/wide_rank_temporal_20260906/wide_rank_results.json")
    points = {
        "Deep ensemble": (full["parameter_counts"]["five_member_deep_ensemble"], 0.8025, 0.7485),
        "Small SBNN": (2212, 0.8210, 0.6020),
        "Dense BBB": (full["parameter_counts"]["dense_bbb_variational_scalars"], 0.8300, 0.6665),
        "Wide rank-50 SBNN": (wide["parameter_counts"]["wide_low_rank"]["(50, 50, 1)"], 0.8435, 0.6284),
        "Wide dense BBB": (wide["parameter_counts"]["wide_dense_bbb"], 0.8135, 0.6567),
    }
    fig, axis = plt.subplots(figsize=(9, 6))
    for name, (parameters, top1, ood) in points.items():
        axis.scatter(parameters, top1, s=90 + 500 * (ood - 0.58), c=ood, vmin=0.58, vmax=0.76, cmap="viridis", edgecolor="black")
        axis.annotate(name, (parameters, top1), xytext=(6, 5), textcoords="offset points", fontsize=9)
    axis.set_xscale("log")
    axis.set_xlabel("Fitted scalar parameters (log scale)")
    axis.set_ylabel("Top-1 accuracy among 100 candidates")
    axis.set_ylim(0.79, 0.855)
    axis.set_title("More parameters improve ranking only up to a point")
    axis.grid(alpha=0.2)
    colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0.58, 0.76), cmap="viridis"), ax=axis)
    colorbar.set_label("Unseen-gene OOD AUROC")
    save(fig, output, "parameter_tradeoff.png")


def external_generalization(output: Path) -> None:
    report = read_json(
        "artifacts/results/phenopacket_hpo2023_generalization_repeated_20260906/repeated_seed_summary.json"
    )["summary"]
    names = list(report)
    labels = [name.replace("_", "\n") for name in names]
    x = np.arange(len(names))
    width = 0.36
    molecular = [report[name]["molecular_top1"]["mean"] for name in names]
    molecular_sd = [report[name]["molecular_top1"]["std"] for name in names]
    hpo = [report[name]["hpo_conditioned_top1"]["mean"] for name in names]
    hpo_sd = [report[name]["hpo_conditioned_top1"]["std"] for name in names]
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    axes[0].bar(x - width / 2, molecular, width, yerr=molecular_sd, capsize=3, label="molecular only", color="#496A81")
    axes[0].bar(x + width / 2, hpo, width, yerr=hpo_sd, capsize=3, label="molecular + HPO", color="#D17B49")
    axes[0].axhline(0.1, color="#333333", linestyle="--", linewidth=1, label="random")
    axes[0].set_xticks(x, labels, fontsize=8)
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("Top-1 accuracy")
    axes[0].set_title("All 719 temporally new solved cases")
    axes[0].legend()
    pair_molecular = [report[name]["molecular_two_variant_top1"]["mean"] for name in names]
    pair_hpo = [report[name]["hpo_conditioned_two_variant_top1"]["mean"] for name in names]
    axes[1].bar(x - width / 2, pair_molecular, width, label="molecular only", color="#496A81")
    axes[1].bar(x + width / 2, pair_hpo, width, label="molecular + HPO", color="#D17B49")
    axes[1].axhline(0.1, color="#333333", linestyle="--", linewidth=1)
    axes[1].set_xticks(x, labels, fontsize=8)
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Top-1 accuracy")
    axes[1].set_title("64 challenge-like two-variant cases")
    axes[1].legend()
    for axis in axes:
        axis.grid(axis="y", alpha=0.2)
    fig.suptitle("External Phenopacket generalization across five resampled candidate cohorts", fontsize=14)
    save(fig, output, "external_generalization.png")


def phenotype_coverage(output: Path) -> None:
    report = read_json(
        "artifacts/results/phenopacket_hpo2023_generalization_repeated_20260906/repeated_seed_summary.json"
    )["summary"]
    names = list(report)
    x = np.arange(len(names))
    width = 0.36
    known = [report[name]["hpo_conditioned_known_gene_top1"]["mean"] for name in names]
    unseen = [report[name]["hpo_conditioned_unseen_gene_top1"]["mean"] for name in names]
    fig, axis = plt.subplots(figsize=(13.5, 5.8))
    axis.bar(x - width / 2, known, width, label="causal gene represented in older phenotype cases", color="#4F8A5B")
    axis.bar(x + width / 2, unseen, width, label="causal gene absent from phenotype reference", color="#9A6B78")
    axis.set_xticks(x, [name.replace("_", "\n") for name in names], fontsize=8)
    axis.set_ylim(0, 1)
    axis.set_ylabel("Phenotype-conditioned top-1 accuracy")
    axis.set_title("Phenotype matching is powerful for known genes, but it does not solve novel genes")
    axis.grid(axis="y", alpha=0.2)
    axis.legend()
    save(fig, output, "phenotype_coverage.png")


def phenotype_robustness(output: Path) -> None:
    report = read_json(
        "artifacts/results/patient_generalization_protocol_20260907/evaluation/patient_generalization_results.json"
    )["results"]
    conditions = ["hpo_original", "hpo_drop20", "hpo_drop40", "hpo_noise20"]
    labels = ["Original", "Drop 20%", "Drop 40%", "Add 20% noise"]
    models = {
        "Wide rank-50 SBNN": "wide_rank50_sbnn",
        "Dense BBB": "dense_bbb",
        "Deep ensemble": "deep_ensemble",
        "Temporal consensus": "temporal_consensus",
    }
    colors = ["#B85C38", "#8B5E83", "#477998", "#4F8A5B"]
    fig, axis = plt.subplots(figsize=(10.5, 5.8))
    x = np.arange(len(conditions))
    for (label, model), color in zip(models.items(), colors):
        values = [report[condition][model]["top1_end_to_end"] for condition in conditions]
        axis.plot(x, values, marker="o", linewidth=2.2, label=label, color=color)
    axis.set_xticks(x, labels)
    axis.set_ylim(0.45, 0.63)
    axis.set_ylabel("Top-1 causal-hypothesis accuracy")
    axis.set_title("Moderate phenotype corruption causes limited degradation")
    axis.grid(alpha=0.22)
    axis.legend(ncol=2)
    axis.text(
        0.01, 0.02,
        "Single fixed candidate cohort; HPO mixing weights were not retuned.",
        transform=axis.transAxes,
        fontsize=9,
        color="#555555",
    )
    save(fig, output, "phenotype_robustness.png")


def varpb_full_background(output: Path) -> None:
    learned = read_json("artifacts/results/varpb_2023_opencravat_gnomad/evaluation/summary.json")
    external = read_json("artifacts/results/varpb_2023_external_predictors/evaluation/summary.json")
    augmented = read_json("artifacts/results/varpb_augmentation_20260908/evaluation/summary.json")
    rows = learned + external + augmented
    preferred_order = [
        "bonamame_small_sbnn_score",
        "bonamame_wide_rank50_sbnn_score",
        "bonamame_augmented_pairwise_sbnn_score",
        "bonamame_temporal_consensus_score",
        "bonamame_dense_bbb_score",
        "bonamame_deep_ensemble_score",
        "bonamame_gradient_boost_score",
        "bonamame_revel_score",
        "bonamame_alphamissense_score",
    ]
    by_name = {row["model"]: row for row in rows}
    selected = [by_name[name] for name in preferred_order]
    labels = [
        "Small SBNN",
        "Wide rank-50 SBNN",
        "Augmented pairwise SBNN",
        "Temporal consensus",
        "Dense BBB",
        "Deep ensemble",
        "Gradient boost",
        "REVEL",
        "AlphaMissense",
    ]
    colors = ["#8B5E83", "#B85C38", "#A13D2D", "#4F8A5B", "#8B5E83", "#477998", "#477998", "#D49A3A", "#D49A3A"]
    y = np.arange(len(selected))
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.2))
    median_ranks = [row["overall"]["median_rank"] for row in selected]
    top20 = [row["overall"]["top20"] for row in selected]
    axes[0].barh(y, median_ranks, color=colors)
    axes[0].set_xscale("log")
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Median causal-variant rank (log scale; lower is better)")
    axes[0].set_title("Rank in 107 real genome backgrounds")
    for index, value in enumerate(median_ranks):
        axes[0].text(value * 1.08, index, f"{value:g}", va="center", fontsize=9)
    axes[1].barh(y, top20, color=colors)
    axes[1].set_yticks(y, [""] * len(y))
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, 0.72)
    axes[1].set_xlabel("Fraction of 275,739 spike-in trials ranked in Top 20")
    axes[1].set_title("Coverage-preserving Top-20 performance")
    for index, value in enumerate(top20):
        axes[1].text(value + 0.012, index, f"{100 * value:.1f}%", va="center", fontsize=9)
    for axis in axes:
        axis.grid(axis="x", alpha=0.2)
    fig.suptitle("VarPB v3 reveals the gap between compact candidate reranking and full-background prioritization", fontsize=13.5)
    save(fig, output, "varpb_full_background.png")


def augmented_scalability(output: Path) -> None:
    report = read_json("artifacts/results/varpb_augmentation_20260908/training_report.json")
    counts = report["parameter_counts"]
    methods = ["Rank-50 SBNN", "Dense BBB", "Five-network ensemble"]
    values = [
        counts["augmented_sbnn_variational_scalars"],
        counts["matched_dense_bbb_variational_scalars"],
        counts["matched_five_network_ensemble_scalars"],
    ]
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8), gridspec_kw={"width_ratios": [1.0, 1.35]})
    y = np.arange(len(methods))
    axes[0].barh(y, values, color=["#A13D2D", "#8B5E83", "#477998"])
    axes[0].set_yticks(y, methods)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Stored fitted scalars")
    axes[0].set_title("Matched 54-200-100-1 architectures")
    axes[0].grid(axis="x", alpha=0.2)
    for index, value in enumerate(values):
        axes[0].text(value + 3000, index, f"{value:,}", va="center", fontsize=9)
    axes[0].set_xlim(0, 175000)

    capability_methods = ["Rank-50 SBNN", "Dense BBB", "Deep ensemble", "REVEL", "AlphaMissense"]
    capabilities = ["Posterior\nweight draws", "Predictive\ndistribution", "Epistemic\nuncertainty proxy"]
    matrix = np.asarray([
        [1, 1, 1],
        [1, 1, 1],
        [0, 1, 1],
        [0, 0, 0],
        [0, 0, 0],
    ])
    axes[1].imshow(matrix, cmap=plt.matplotlib.colors.ListedColormap(["#ECECEC", "#4F8A5B"]), vmin=0, vmax=1, aspect="auto")
    axes[1].set_xticks(np.arange(len(capabilities)), capabilities)
    axes[1].set_yticks(np.arange(len(capability_methods)), capability_methods)
    axes[1].set_title("Uncertainty available from the standard model output")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            axes[1].text(column, row, "yes" if matrix[row, column] else "no", ha="center", va="center", color="white" if matrix[row, column] else "#555555", fontweight="bold")
    fig.suptitle("The SBNN adds a compact posterior; external pathogenicity scores remain deterministic", fontsize=13.5)
    save(fig, output, "augmented_scalability_uncertainty.png")


def augmented_generalization_pipeline(output: Path) -> None:
    fig, axis = plt.subplots(figsize=(15.5, 8.2))
    axis.set_xlim(0, 15.5)
    axis.set_ylim(0, 8.2)
    axis.axis("off")

    def box(x, y, width, height, title, body, color):
        patch = FancyBboxPatch(
            (x, y), width, height,
            boxstyle="round,pad=0.03,rounding_size=0.08",
            linewidth=1.5, edgecolor=color, facecolor=color + "18",
        )
        axis.add_patch(patch)
        axis.text(x + width / 2, y + height * 0.72, title, ha="center", va="center", fontsize=10.5, fontweight="bold", color=color)
        axis.text(x + width / 2, y + height * 0.36, body, ha="center", va="center", fontsize=9)

    def arrow(x1, y1, x2, y2):
        axis.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops={"arrowstyle": "->", "lw": 1.7, "color": "#555555"})

    axis.text(7.75, 7.85, "Predictor-augmented SBNN: training, frozen transfer, and evaluation", ha="center", fontsize=15, fontweight="bold")
    axis.text(0.25, 6.95, "TRAINING AND MODEL SELECTION", fontsize=10, fontweight="bold", color="#555555")
    box(0.25, 5.45, 2.55, 1.15, "ClinVar 2023", "146,810 labeled\nmissense SNVs", "#284B63")
    box(3.15, 5.45, 2.75, 1.15, "54 features", "50 molecular features +\nREVEL/AlphaMissense + missingness", "#356B83")
    box(6.25, 5.45, 2.55, 1.15, "Leakage control", "80/20 gene-disjoint split;\nno VarPB or challenge labels", "#4F8A5B")
    box(9.15, 5.45, 2.65, 1.15, "Factorial selection", "24 rank, KL, posterior-scale,\nand prior-scale settings", "#8B5E83")
    box(12.15, 5.35, 3.05, 1.35, "Frozen SBNN", "54-200-100-1; ranks 50,50,1\nKL 0.001; shared-draw RankNet", "#A13D2D")
    for left, right in [(2.8, 3.15), (5.9, 6.25), (8.8, 9.15), (11.8, 12.15)]:
        arrow(left, 6.02, right, 6.02)

    axis.text(0.25, 4.45, "FROZEN EVALUATION", fontsize=10, fontweight="bold", color="#555555")
    box(0.25, 2.85, 2.7, 1.15, "VarPB background test", "2,577 causal variants x\n107 real genome backgrounds", "#D49A3A")
    box(3.35, 2.85, 2.55, 1.15, "Allele ranking", "posterior mean score and\nposterior uncertainty", "#A13D2D")
    box(6.30, 2.85, 2.75, 1.15, "Generalization result", "Top-20 65.70%; Top-100 86.63%\nmedian rank 14", "#4F8A5B")
    arrow(2.95, 3.42, 3.35, 3.42)
    arrow(5.90, 3.42, 6.30, 3.42)
    arrow(13.68, 5.35, 5.2, 4.0)

    box(0.25, 0.75, 2.7, 1.15, "Track 1 transfer", "two alleles per retained\ncompound-het hypothesis", "#477998")
    box(3.35, 0.75, 2.55, 1.15, "Pair integration", "mean over shared draws of\nHPO x sqrt(p1 p2)", "#8B5E83")
    box(6.30, 0.75, 2.75, 1.15, "Challenge result", "BUB1B pair at rank 1\nRank Points 100; F-max 1.000", "#4F8A5B")
    arrow(2.95, 1.32, 3.35, 1.32)
    arrow(5.90, 1.32, 6.30, 1.32)
    arrow(13.68, 5.35, 5.2, 1.90)

    axis.text(
        11.9, 2.45,
        "One frozen posterior, two complementary tests:\n"
        "genome-background transfer measures shortlist recall;\n"
        "Track 1 verifies that the revision preserves the child-level solution.",
        ha="center", va="center", fontsize=10.5,
        bbox={"boxstyle": "round,pad=0.5", "facecolor": "#F4F1EA", "edgecolor": "#777777"},
    )
    save(fig, output, "augmented_generalization_pipeline.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=FIGURES / "report")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    generated = 0
    skipped = []
    for builder in (
        system_pipeline,
        analysis_funnel,
        track1_scores,
        temporal_model_comparison,
        parameter_tradeoff,
        external_generalization,
        phenotype_coverage,
        phenotype_robustness,
        varpb_full_background,
        augmented_generalization_pipeline,
        augmented_scalability,
    ):
        try:
            builder(args.output_dir)
            generated += 1
        except (FileNotFoundError, KeyError) as exc:
            missing = exc.filename if isinstance(exc, FileNotFoundError) else str(exc)
            skipped.append(f"{builder.__name__}: {missing}")
    print(f"Generated {generated} figures in {args.output_dir}")
    for item in skipped:
        print(f"SKIPPED (required artifact unavailable): {item}")


if __name__ == "__main__":
    main()
