from __future__ import annotations

import os
import pickle

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "8")

import numpy as np

from rare_disease_real_kid.checkpoint import load_checkpoint
from rare_disease_real_kid.experiments.clinvar_temporal import Variant
from rare_disease_real_kid.experiments.phenopacket_temporal import predict_all_statistics
from rare_disease_real_kid.paths import CHECKPOINTS


def _records() -> list[Variant]:
    return [
        Variant(
            variation_id=index,
            chrom="1",
            pos=1000 + index,
            ref="A",
            alt="G",
            gene=f"GENE{index % 2}",
            label=index % 2,
            consequence="missense_variant",
            variant_type="single_nucleotide_variant",
            review="",
            origin=1,
            af_esp=0.0001,
            af_exac=None,
            af_tgp=None,
            has_rs=False,
        )
        for index in range(1, 7)
    ]


def test_converted_checkpoint_matches_legacy_predictions() -> None:
    legacy = load_checkpoint(CHECKPOINTS / "temporal" / "trained_models.legacy.pkl")
    with (CHECKPOINTS / "temporal" / "trained_models.pkl").open("rb") as handle:
        canonical = pickle.load(handle)

    expected = predict_all_statistics(legacy, _records(), neutralize_gene_history=True)
    actual = predict_all_statistics(canonical, _records(), neutralize_gene_history=True)
    assert set(actual) == set(expected)
    for name in expected:
        np.testing.assert_allclose(actual[name][0], expected[name][0], rtol=0, atol=0)
        np.testing.assert_allclose(
            actual[name][1], expected[name][1], rtol=0, atol=0, equal_nan=True
        )


def test_augmented_pairwise_checkpoint_is_loadable() -> None:
    artifact = load_checkpoint(
        CHECKPOINTS / "augmented_pairwise_sbnn" / "augmented_sbnn.pkl"
    )
    assert set(artifact) == {
        "augmented_features",
        "config",
        "gene_statistics",
        "model",
        "predictor_columns",
        "scaler",
        "seed",
    }
    assert artifact["config"]["ranks"] == (50, 50, 1)
    assert artifact["seed"] == 20260908
