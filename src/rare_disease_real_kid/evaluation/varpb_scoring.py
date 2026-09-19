#!/usr/bin/env python3
"""Score an annotated VarPB variant table with every saved model family.

VarPB uses hg19 coordinates, but coordinates are identifiers rather than model
features.  Consequence, allele frequency and gene are read from an annotation
table; ClinVar review fields and gene-history features are deliberately neutral.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import pickle
from contextlib import ExitStack
from pathlib import Path

from rare_disease_real_kid.experiments.clinvar_temporal import Variant
from rare_disease_real_kid.experiments.phenopacket_temporal import predict_all_statistics
from rare_disease_real_kid.checkpoint import load_checkpoint


REQUIRED = {"Chromosome", "Position", "Ref", "Alt", "Gene", "Consequence", "PopulationAF"}


def official_chromosome(value: str) -> str:
    """Return the chromosome spelling expected by VarPB's reference script."""
    chromosome = value.strip()
    return chromosome if chromosome.startswith("chr") else f"chr{chromosome}"


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(mode.removesuffix("t"), newline="", encoding="utf-8")


def optional_frequency(value: str) -> float | None:
    if value.strip() in {"", ".", "NA", "None"}:
        return None
    result = float(value)
    if not 0 <= result <= 1:
        raise ValueError(f"PopulationAF must be in [0, 1], got {value!r}")
    return result


def as_variant(row: dict[str, str], identifier: int) -> Variant:
    ref = row["Ref"].upper()
    alt = row["Alt"].upper()
    if not ref or not alt or "," in alt:
        raise ValueError(f"Expected one normalized alternate allele: {row}")
    frequency = optional_frequency(row["PopulationAF"])
    return Variant(
        variation_id=identifier,
        chrom=row["Chromosome"].removeprefix("chr"),
        pos=int(row["Position"]),
        ref=ref,
        alt=alt,
        gene=row["Gene"] or "UNKNOWN",
        label=0,
        consequence=row["Consequence"].lower(),
        variant_type="single_nucleotide_variant" if len(ref) == len(alt) == 1 else "indel",
        review="",
        origin=1,
        af_esp=frequency,
        af_exac=None,
        af_tgp=None,
        has_rs=False,
    )


def chunks(reader, size: int):
    block = []
    for row in reader:
        block.append(row)
        if len(block) == size:
            yield block
            block = []
    if block:
        yield block


def score_table(annotated: Path, checkpoint: Path, output_dir: Path, chunk_size: int = 50000) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    models = load_checkpoint(checkpoint, allow_external=True)

    total = 0
    model_names = None
    with open_text(annotated, "rt") as input_handle, ExitStack() as stack:
        reader = csv.DictReader(input_handle, delimiter="\t")
        missing = REQUIRED.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Annotation table is missing columns: {sorted(missing)}")
        writers = {}
        uncertainty_writer = None
        for block in chunks(reader, chunk_size):
            records = [as_variant(row, total + index + 1) for index, row in enumerate(block)]
            statistics = predict_all_statistics(models, records, neutralize_gene_history=True)
            if model_names is None:
                model_names = sorted(statistics)
                for name in model_names:
                    handle = stack.enter_context((output_dir / f"bonamame_{name}_varpb.tsv").open("w", newline="", encoding="utf-8"))
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(["Chromosome", "Position", "Ref", "Alt", f"bonamame_{name}_score"])
                    writers[name] = writer
                uncertainty_handle = stack.enter_context(gzip.open(output_dir / "bonamame_varpb_uncertainty.tsv.gz", "wt", newline="", encoding="utf-8"))
                uncertainty_writer = csv.writer(uncertainty_handle, delimiter="\t")
                uncertainty_writer.writerow(["Chromosome", "Position", "Ref", "Alt", "model", "uncertainty"])
            for name in model_names:
                means, uncertainty = statistics[name]
                for row, score, spread in zip(block, means, uncertainty):
                    key = [official_chromosome(row["Chromosome"]), row["Position"], row["Ref"], row["Alt"]]
                    writers[name].writerow([*key, float(score)])
                    if uncertainty_writer is not None and spread == spread:
                        uncertainty_writer.writerow([*key, name, float(spread)])
            total += len(block)
            print(f"Scored {total:,} variants", flush=True)
    return {
        "variants": total,
        "models": model_names or [],
        "coordinate_system": "hg19 identifiers supplied by VarPB",
        "gene_history_neutralized": True,
        "clinvar_review_features_available": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotated-variants", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=50000)
    args = parser.parse_args()
    if args.chunk_size < 1:
        parser.error("--chunk-size must be positive")
    print(score_table(args.annotated_variants, args.checkpoint, args.output_dir, args.chunk_size))


if __name__ == "__main__":
    main()
