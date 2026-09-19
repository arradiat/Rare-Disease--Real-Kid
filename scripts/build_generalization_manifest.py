#!/usr/bin/env python3
"""Build the integrity manifest for committed generalization evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "artifacts" / "results"
CHECKPOINTS = ROOT / "artifacts" / "checkpoints"
OUTPUT = ROOT / "artifacts" / "manifests" / "generalization_artifacts.json"

RESULT_GROUPS = (
    "full_bbb_temporal_20260906",
    "patient_generalization_protocol_20260907",
    "phenopacket_hpo2023_generalization_repeated_20260906",
    "phenopacket_hpo2023_generalization_seed_20260905",
    "phenopacket_hpo2023_generalization_seed_20260906",
    "phenopacket_hpo2023_generalization_seed_20260907",
    "phenopacket_hpo2023_generalization_seed_20260908",
    "phenopacket_hpo2023_generalization_seed_20260909",
    "phenopacket_hpo_generalization_20260906",
    "phenopacket_hpo_generalization_repeated_20260906",
    "phenopacket_hpo_generalization_seed_20260905",
    "phenopacket_hpo_generalization_seed_20260906",
    "phenopacket_hpo_generalization_seed_20260907",
    "phenopacket_hpo_generalization_seed_20260908",
    "phenopacket_hpo_generalization_seed_20260909",
    "phenopacket_temporal_generalization_20260906",
    "sbnn_broad_uncertainty_sweep_20260906",
    "sbnn_factorial_fast_sweep_20260906",
    "sbnn_factorial_pareto_complete_20260906",
    "sbnn_factorial_sweep_20260906",
    "sbnn_hyperparameter_sweep_20260905",
    "sbnn_hyperparameter_sweep_pareto_20260906",
    "sbnn_hyperparameter_sweep_uncertainty_20260906",
    "sbnn_run_20260905",
    "varpb_2023_declared_missense",
    "varpb_2023_external_predictors",
    "varpb_2023_opencravat_gnomad",
    "varpb_augmentation_20260908",
    "track1_augmented_sbnn_20260908",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(path: Path, group: str) -> dict[str, object]:
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "group": group,
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def main() -> None:
    files: list[dict[str, object]] = []
    for group in RESULT_GROUPS:
        directory = RESULTS / group
        if not directory.is_dir():
            raise FileNotFoundError(f"Missing generalization result directory: {directory}")
        files.extend(
            record(path, group)
            for path in sorted(directory.rglob("*"))
            if path.is_file() and path.name != ".DS_Store"
        )

    checkpoint = CHECKPOINTS / "augmented_pairwise_sbnn" / "augmented_sbnn.pkl"
    if not checkpoint.is_file():
        raise FileNotFoundError(f"Missing augmented SBNN checkpoint: {checkpoint}")
    files.append(record(checkpoint, "augmented_pairwise_sbnn_checkpoint"))
    files.sort(key=lambda item: str(item["path"]))

    payload = {
        "schema_version": 1,
        "description": (
            "Integrity inventory for the committed Phenopacket/HPO, VarPB, "
            "predictor-augmented SBNN, and frozen Track 1 transfer artifacts."
        ),
        "file_count": len(files),
        "total_size_bytes": sum(int(item["size_bytes"]) for item in files),
        "files": files,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(files)} records to {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
