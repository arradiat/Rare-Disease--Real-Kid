#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"
export PYTHONPATH="$repo_dir/src${PYTHONPATH:+:$PYTHONPATH}"

python -m rare_disease_real_kid.reporting.environment
python -m rare_disease_real_kid.reporting.figures
python -m pytest -q
ranking_zip="${RANKING_ZIP:-../ranking.zip}"
if [[ -f "$ranking_zip" ]]; then
  python -m rare_disease_real_kid submission generate \
    --ranking-zip "$ranking_zip" \
    --checkpoint artifacts/checkpoints/temporal/trained_models.pkl
else
  echo "No ranking.zip found; using the committed audited SBNN prediction."
  python -m rare_disease_real_kid submission generate
fi
python -m rare_disease_real_kid.evaluation.track1_scores

build_dir="$repo_dir/artifacts/reports/build"
mkdir -p "$build_dir"
cd report
pdflatex -output-directory="$build_dir" -interaction=nonstopmode -halt-on-error main.tex
cd "$build_dir"
BIBINPUTS="$repo_dir/report:" bibtex main
cd "$repo_dir/report"
pdflatex -output-directory="$build_dir" -interaction=nonstopmode -halt-on-error main.tex
pdflatex -output-directory="$build_dir" -interaction=nonstopmode -halt-on-error main.tex

cd "$repo_dir"
cp "$build_dir/main.pdf" artifacts/reports/rebuilt_from_source.pdf

echo "Built artifacts/reports/rebuilt_from_source.pdf and regenerated available verification artifacts."
