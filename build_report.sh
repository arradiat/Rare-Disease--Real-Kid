#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo_dir"

python capture_environment.py
python report_visualizations.py --output-dir report/figures
python -m pytest -q
ranking_zip="${RANKING_ZIP:-../ranking.zip}"
if [[ -f "$ranking_zip" ]]; then
  python generate_track1_submission.py \
    --ranking-zip "$ranking_zip" \
    --checkpoint results/track1_temporal_transfer_20260906/trained_models.pkl
  python evaluate_augmented_sbnn_track1.py \
    --ranking-zip "$ranking_zip" \
    --checkpoint results/varpb_augmentation_20260908/augmented_sbnn.pkl \
    --revel-database data/open/opencravat_modules/annotators/revel/data/revel.sqlite \
    --alphamissense-database data/open/opencravat_modules/annotators/alphamissense/data/alphamissense.sqlite \
    --output-dir results/track1_augmented_sbnn_20260908
else
  echo "No ranking.zip found; using the committed audited SBNN prediction."
  python generate_track1_submission.py
fi
python verify_track1_exact_scores.py

cd report
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex

cd "$repo_dir"
mkdir -p output/pdf
cp report/main.pdf output/pdf/bonamame_track1_report.pdf

echo "Built output/pdf/bonamame_track1_report.pdf and regenerated all verification artifacts."
