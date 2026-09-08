#!/usr/bin/env python3
"""Prepare a leakage-controlled ClinVar missense set for predictor augmentation.

Only unambiguous pathogenic/likely-pathogenic or benign/likely-benign,
biallelic SNVs whose molecular consequence includes ``missense_variant`` are
retained. VarPB records are never read by this script.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from clinvar_temporal_experiment import clean_label, parse_info


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, encoding="utf-8") if path.suffix == ".gz" else path.open(mode.removesuffix("t"), encoding="utf-8")


def filter_labeled_missense(source: Path, output: Path, manifest: Path) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    counts = {
        "input_records": 0,
        "unambiguous_labeled": 0,
        "biallelic_snv": 0,
        "labeled_missense_snv": 0,
        "positive": 0,
        "negative": 0,
    }
    with open_text(source, "rt") as incoming, open_text(output, "wt") as outgoing:
        outgoing.write("##fileformat=VCFv4.2\n")
        outgoing.write("##reference=GRCh38\n")
        outgoing.write("##source=ClinVar-20231226-labeled-missense-filter\n")
        outgoing.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for line in incoming:
            if line.startswith("#"):
                continue
            counts["input_records"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[2].isdigit():
                continue
            info = parse_info(fields[7])
            label = clean_label(info.get("CLNSIG"))
            if label is None:
                continue
            counts["unambiguous_labeled"] += 1
            ref, alt = fields[3].upper(), fields[4].upper()
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            counts["biallelic_snv"] += 1
            if "missense_variant" not in info.get("MC", "").lower():
                continue
            counts["labeled_missense_snv"] += 1
            counts["positive" if label else "negative"] += 1
            outgoing.write("\t".join(fields[:8]) + "\n")
    result = {
        "source": str(source),
        "output": str(output),
        "selection": "unambiguous ClinVar label AND biallelic SNV AND MC contains missense_variant",
        "varpb_used_for_selection": False,
        **counts,
    }
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    print(filter_labeled_missense(args.source, args.output, args.manifest))


if __name__ == "__main__":
    main()
