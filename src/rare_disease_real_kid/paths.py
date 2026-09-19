"""One source of truth for repository and artifact paths."""

from __future__ import annotations

import os
from pathlib import Path


def find_project_root(start: Path | None = None) -> Path:
    """Find the checkout containing ``pyproject.toml``.

    ``RDKID_PROJECT_ROOT`` takes precedence so installed commands can operate on a
    separately located artifact checkout.
    """
    override = os.environ.get("RDKID_PROJECT_ROOT")
    if override:
        root = Path(override).expanduser().resolve()
        if not (root / "pyproject.toml").is_file():
            raise FileNotFoundError(f"RDKID_PROJECT_ROOT is not a project checkout: {root}")
        return root

    origin = (start or Path(__file__)).resolve()
    candidates = [origin] if origin.is_dir() else [origin.parent]
    candidates.extend(candidates[0].parents)
    for candidate in candidates:
        if (candidate / "pyproject.toml").is_file():
            return candidate
    raise FileNotFoundError(
        "Could not locate pyproject.toml. Run inside the repository or set "
        "RDKID_PROJECT_ROOT."
    )


PROJECT_ROOT = find_project_root()
ARTIFACTS = PROJECT_ROOT / "artifacts"
RESULTS = ARTIFACTS / "results"
CHECKPOINTS = ARTIFACTS / "checkpoints"
FIGURES = ARTIFACTS / "figures"
REPORTS = ARTIFACTS / "reports"
MANIFESTS = ARTIFACTS / "manifests"
SUBMISSION = PROJECT_ROOT / "submission"
DATA = Path(os.environ.get("RDKID_DATA_DIR", PROJECT_ROOT / "data")).expanduser()


def artifact_path(*parts: str) -> Path:
    return ARTIFACTS.joinpath(*parts)

