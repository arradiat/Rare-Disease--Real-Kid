"""Controlled loading for the bundled legacy model checkpoint."""

from __future__ import annotations

import importlib
import pickle
from pathlib import Path
from typing import BinaryIO, Any

from rare_disease_real_kid.paths import ARTIFACTS


LEGACY_MODULES = {
    "clinvar_temporal_experiment": "rare_disease_real_kid.experiments.clinvar_temporal",
    "sbnn_pair_ranker": "rare_disease_real_kid.models.bayesian_layers",
}


class LegacyCheckpointUnpickler(pickle.Unpickler):
    """Map the two historical project-module names to the packaged modules."""

    def find_class(self, module: str, name: str) -> Any:
        mapped = LEGACY_MODULES.get(module, module)
        if mapped != module:
            return getattr(importlib.import_module(mapped), name)
        return super().find_class(mapped, name)


def _inside(path: Path, directory: Path) -> bool:
    try:
        path.resolve().relative_to(directory.resolve())
        return True
    except ValueError:
        return False


def load_checkpoint(path: Path, *, allow_external: bool = False) -> Any:
    """Load a trusted bundled checkpoint or an explicitly authorized external one."""
    path = path.expanduser().resolve()
    if not _inside(path, ARTIFACTS) and not allow_external:
        raise ValueError(
            "Refusing to unpickle an external checkpoint. Only bundled artifacts are "
            "trusted by default; pass allow_external=True for a checkpoint you trust."
        )
    with path.open("rb") as handle:
        return LegacyCheckpointUnpickler(handle).load()


def load_checkpoint_file(handle: BinaryIO) -> Any:
    """Load from an already-authorized binary handle."""
    return LegacyCheckpointUnpickler(handle).load()

