"""Verify the repository and its committed scientific artifacts."""

from __future__ import annotations

import hashlib
import importlib
import json
import re
from pathlib import Path

from rare_disease_real_kid.checkpoint import load_checkpoint
from rare_disease_real_kid.evaluation.submission_audit import audit
from rare_disease_real_kid.evaluation.track1_scores import score
from rare_disease_real_kid.models.registry import MODEL_REGISTRY, SELECTED_SUBMISSION_MODEL
from rare_disease_real_kid.paths import CHECKPOINTS, MANIFESTS, PROJECT_ROOT, SUBMISSION


CANONICAL_SUBMISSION_SHA256 = "308eb684858fc9bfc97cd8b81f0689914408f78239a1b42563419498157510eb"
EXPECTED_MODELS = {
    "gradient_boost",
    "dense_mlp",
    "deep_ensemble",
    "small_sbnn",
    "dense_bbb",
    "wide_rank50_sbnn",
    "wide_dense_bbb",
    "temporal_consensus",
    "augmented_pairwise_sbnn",
}

EXPECTED_AUGMENTED_CHECKPOINT_KEYS = {
    "augmented_features",
    "config",
    "gene_statistics",
    "model",
    "predictor_columns",
    "scaler",
    "seed",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def markdown_link_errors() -> list[str]:
    errors = []
    pattern = re.compile(r"\[[^\]]*\]\(([^)]+)\)")
    for document in sorted([PROJECT_ROOT / "README.md", *PROJECT_ROOT.joinpath("docs").glob("*.md"), PROJECT_ROOT / "submission/README.md", PROJECT_ROOT / "data/README.md"]):
        for target in pattern.findall(document.read_text(encoding="utf-8")):
            if "://" in target or target.startswith(("#", "mailto:")):
                continue
            local = target.split("#", 1)[0]
            if local and not (document.parent / local).resolve().exists():
                errors.append(f"{document.relative_to(PROJECT_ROOT)} -> {target}")
    return errors


def verify_repository() -> dict[str, object]:
    checks: dict[str, object] = {}
    errors: list[str] = []

    manifest_path = MANIFESTS / "original_file_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    checks["original_file_count"] = manifest["original_file_count"]
    if manifest["original_file_count"] != 131:
        errors.append("Original manifest does not contain 131 files")
    missing_destinations = sorted(
        record["destination_path"]
        for record in manifest["files"]
        if not (PROJECT_ROOT / record["destination_path"]).exists()
    )
    checks["missing_migrated_destinations"] = missing_destinations
    errors.extend(f"Missing migrated destination: {path}" for path in missing_destinations)
    final_hash_errors = []
    for record in manifest["files"]:
        destination = PROJECT_ROOT / record["destination_path"]
        if destination.exists() and record.get("final_sha256") != sha256(destination):
            final_hash_errors.append(record["destination_path"])
    checks["final_destination_hash_errors"] = sorted(set(final_hash_errors))
    errors.extend(f"Final destination hash mismatch: {path}" for path in sorted(set(final_hash_errors)))

    backup = PROJECT_ROOT.with_name(manifest["source_snapshot"])
    backup_errors = []
    if backup.exists():
        for record in manifest["files"]:
            source = backup / record["original_path"]
            if not source.exists() or sha256(source) != record["sha256"]:
                backup_errors.append(record["original_path"])
    checks["backup_available"] = backup.exists()
    checks["backup_hash_errors"] = backup_errors
    errors.extend(f"Backup hash mismatch: {path}" for path in backup_errors)

    generalization_manifest_path = MANIFESTS / "generalization_artifacts.json"
    generalization_manifest = json.loads(
        generalization_manifest_path.read_text(encoding="utf-8")
    )
    missing_generalization_artifacts = []
    generalization_hash_errors = []
    generalization_size = 0
    for record in generalization_manifest["files"]:
        artifact = PROJECT_ROOT / record["path"]
        if not artifact.is_file():
            missing_generalization_artifacts.append(record["path"])
            continue
        generalization_size += artifact.stat().st_size
        if sha256(artifact) != record["sha256"]:
            generalization_hash_errors.append(record["path"])
    checks["generalization_artifact_count"] = len(generalization_manifest["files"])
    checks["generalization_artifact_bytes"] = generalization_size
    checks["missing_generalization_artifacts"] = missing_generalization_artifacts
    checks["generalization_hash_errors"] = generalization_hash_errors
    if generalization_manifest["file_count"] != len(generalization_manifest["files"]):
        errors.append("Generalization artifact manifest file count is inconsistent")
    if generalization_manifest["total_size_bytes"] != generalization_size:
        errors.append("Generalization artifact manifest byte count is inconsistent")
    errors.extend(
        f"Missing generalization artifact: {path}"
        for path in missing_generalization_artifacts
    )
    errors.extend(
        f"Generalization artifact hash mismatch: {path}"
        for path in generalization_hash_errors
    )

    root_python = sorted(path.name for path in PROJECT_ROOT.glob("*.py"))
    checks["root_python_files"] = root_python
    if root_python:
        errors.append(f"Loose Python files remain at repository root: {root_python}")

    checks["model_registry"] = sorted(MODEL_REGISTRY)
    if set(MODEL_REGISTRY) != EXPECTED_MODELS:
        errors.append("Model registry does not contain the nine documented models")
    for specification in MODEL_REGISTRY.values():
        importlib.import_module(specification.module)
    if SELECTED_SUBMISSION_MODEL != "wide_rank50_sbnn":
        errors.append("Selected submission model changed")

    checkpoint = CHECKPOINTS / "temporal" / "trained_models.pkl"
    models = load_checkpoint(checkpoint)
    checks["checkpoint_groups"] = sorted(models)

    augmented_checkpoint = (
        CHECKPOINTS / "augmented_pairwise_sbnn" / "augmented_sbnn.pkl"
    )
    augmented_artifact = load_checkpoint(augmented_checkpoint)
    augmented_keys = set(augmented_artifact) if isinstance(augmented_artifact, dict) else set()
    checks["augmented_checkpoint_keys"] = sorted(augmented_keys)
    if augmented_keys != EXPECTED_AUGMENTED_CHECKPOINT_KEYS:
        errors.append("Augmented pairwise SBNN checkpoint has an unexpected structure")

    submission = SUBMISSION / "model1_bonamame_wide-rank50-sbnn.csv"
    submission_hash = sha256(submission)
    checks["submission_sha256"] = submission_hash
    if submission_hash != CANONICAL_SUBMISSION_SHA256:
        errors.append("Canonical submission hash changed")
    schema = audit(submission)
    track1 = score(submission)
    checks["submission_schema"] = schema
    checks["submission_score"] = track1
    if not schema["ok"]:
        errors.append(f"Submission schema failed: {schema['errors']}")
    if not track1["exact_perfect_score"]:
        errors.append("Submission no longer achieves the locally reproduced exact score")

    link_errors = markdown_link_errors()
    checks["broken_markdown_links"] = link_errors
    errors.extend(f"Broken Markdown link: {item}" for item in link_errors)

    checks["ok"] = not errors
    checks["errors"] = errors
    return checks


def main() -> int:
    result = verify_repository()
    output = MANIFESTS / "verification_report.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
