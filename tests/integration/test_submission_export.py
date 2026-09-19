from __future__ import annotations

import csv
import json

from rare_disease_real_kid.evaluation.submission_audit import REQUIRED, audit
from rare_disease_real_kid.evaluation.submission_export import MANIFEST, OUTPUT, export_submission
from rare_disease_real_kid.evaluation.submission_generate import generate_submission
from rare_disease_real_kid.evaluation.track1_scores import score


def test_selected_sbnn_export_is_upload_ready() -> None:
    manifest = export_submission()
    with OUTPUT.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    assert tuple(reader.fieldnames or ()) == REQUIRED
    assert len(rows) == 6
    assert rows[0]["notes"].startswith("gene=BUB1B;model=wide_rank50_sbnn;")
    assert audit(OUTPUT)["ok"]
    assert score(OUTPUT)["exact_perfect_score"]
    assert manifest["model"]["implementation"] == "src/rare_disease_real_kid/models/wide_rank50_sbnn.py"
    assert manifest["generation"]["entry_point"] == "rdkid submission generate"
    assert manifest["generation"]["model_inference"] == "rare_disease_real_kid.evaluation.track1_transfer"
    assert manifest["generation"]["submission_export"] == "rare_disease_real_kid.evaluation.submission_export"
    assert json.loads(MANIFEST.read_text(encoding="utf-8"))["sha256"] == manifest["sha256"]


def test_public_entry_point_exports_the_selected_saved_prediction() -> None:
    manifest = generate_submission()
    assert manifest["model"]["name"] == "wide_rank50_sbnn"
    assert manifest["submission"] == "submission/model1_bonamame_wide-rank50-sbnn.csv"
