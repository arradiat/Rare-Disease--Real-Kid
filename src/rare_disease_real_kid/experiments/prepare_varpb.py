#!/usr/bin/env python3
"""Prepare and export the public VarPB benchmark for local model scoring.

The preparation stage deduplicates the 17.9 million sample-variant rows into a
sample-free hg19 VCF using a disk-backed SQLite index.  The export stage turns
an OpenCRAVAT SQLite result into the seven-column table consumed by
``score_varpb_models.py``.  ClinVar fields are intentionally not exported.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sqlite3
from pathlib import Path


BACKGROUND_COLUMNS = {"hg19_chr", "hg19_pos(1-based)", "ref", "alt", "sampleID"}
CAUSAL_COLUMNS = {"Chromosome", "Position", "ReferenceAllele", "AlternateAllele"}
OUTPUT_COLUMNS = ["Chromosome", "Position", "Ref", "Alt", "Gene", "Consequence", "PopulationAF"]
POPULATION_MODULES = ("gnomad", "thousandgenomes", "esp6500", "alfa")
FREQUENCY_TERMS = ("af", "freq", "frequency")
SO_LABELS = {
    "MIS": "missense_variant",
    "FSD": "frameshift_truncation",
    "FSI": "frameshift_elongation",
    "STG": "stop_gained nonsense",
    "SPL": "splice_region_variant",
    "INI": "inframe_insertion",
    "IND": "inframe_deletion",
    "SYN": "synonymous_variant",
    "INT": "intron_variant",
    "UT3": "3_prime_UTR_variant",
    "UT5": "5_prime_UTR_variant",
    "NMD": "non-coding_NMD_transcript_variant",
}


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="", encoding="utf-8")
    return path.open(mode.removesuffix("t"), newline="", encoding="utf-8")


def clean_variant(chrom: str, pos: str | int, ref: str, alt: str) -> tuple[str, int, str, str]:
    chromosome = str(chrom).strip().removeprefix("chr")
    position = int(pos)
    reference = str(ref).strip().upper()
    alternate = str(alt).strip().upper()
    if not chromosome or position < 1:
        raise ValueError(f"Invalid coordinate: {chrom}:{pos}")
    if not reference or not alternate or "," in alternate:
        raise ValueError(f"Expected one non-empty alternate allele: {chrom}:{pos} {ref}>{alt}")
    if not re.fullmatch(r"[ACGTN]+", reference) or not re.fullmatch(r"[ACGTN]+", alternate):
        raise ValueError(f"Non-VCF allele: {chrom}:{pos} {ref}>{alt}")
    return chromosome, position, reference, alternate


def chrom_rank(chrom: str) -> int:
    if chrom.isdigit():
        return int(chrom)
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(chrom.upper(), 1000)


def insert_rows(connection: sqlite3.Connection, rows, batch_size: int = 100_000) -> tuple[int, int]:
    seen = inserted = 0
    block = []
    for chrom, pos, ref, alt in rows:
        seen += 1
        block.append((chrom, chrom_rank(chrom), pos, ref, alt))
        if len(block) == batch_size:
            before = connection.total_changes
            connection.executemany("INSERT OR IGNORE INTO variants VALUES (?, ?, ?, ?, ?)", block)
            inserted += connection.total_changes - before
            connection.commit()
            block.clear()
            print(f"Read {seen:,} rows; retained {inserted:,} unique variants", flush=True)
    if block:
        before = connection.total_changes
        connection.executemany("INSERT OR IGNORE INTO variants VALUES (?, ?, ?, ?, ?)", block)
        inserted += connection.total_changes - before
        connection.commit()
    return seen, inserted


def background_rows(path: Path):
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = BACKGROUND_COLUMNS.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"VarPB background is missing columns: {sorted(missing)}")
        for row in reader:
            yield clean_variant(row["hg19_chr"], row["hg19_pos(1-based)"], row["ref"], row["alt"])


def causal_rows(paths: list[Path]):
    for path in paths:
        with open_text(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = CAUSAL_COLUMNS.difference(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")
            for row in reader:
                yield clean_variant(
                    row["Chromosome"], row["Position"], row["ReferenceAllele"], row["AlternateAllele"]
                )


def prepare_vcf(background: Path, causal: list[Path], output: Path, index_path: Path) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    index_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(index_path) as connection:
        connection.execute("PRAGMA journal_mode=WAL")
        connection.execute("PRAGMA synchronous=NORMAL")
        connection.execute(
            "CREATE TABLE IF NOT EXISTS variants (chrom TEXT, chrom_rank INTEGER, pos INTEGER, ref TEXT, alt TEXT, "
            "PRIMARY KEY (chrom, pos, ref, alt)) WITHOUT ROWID"
        )
        background_seen, background_added = insert_rows(connection, background_rows(background))
        causal_seen, causal_added = insert_rows(connection, causal_rows(causal))
        unique = int(connection.execute("SELECT COUNT(*) FROM variants").fetchone()[0])
        with open_text(output, "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##reference=GRCh37/hg19\n")
            handle.write("##source=VarPB-v3-union-preparation\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            query = "SELECT chrom, pos, ref, alt FROM variants ORDER BY chrom_rank, chrom, pos, ref, alt"
            for chrom, pos, ref, alt in connection.execute(query):
                handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")
    return {
        "background_rows": background_seen,
        "background_unique_added": background_added,
        "causal_rows": causal_seen,
        "causal_unique_added": causal_added,
        "union_unique_variants": unique,
        "coordinate_system": "GRCh37/hg19",
        "vcf": str(output),
        "deduplication_index": str(index_path),
    }


def table_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [str(row[1]) for row in connection.execute(f'PRAGMA table_info("{table}")')]


def first_present(columns: list[str], *names: str) -> str:
    lookup = {name.lower(): name for name in columns}
    for name in names:
        if name.lower() in lookup:
            return lookup[name.lower()]
    raise ValueError(f"OpenCRAVAT output lacks any of these columns: {names}")


def frequency_columns(columns: list[str]) -> list[str]:
    result = []
    for column in columns:
        lower = column.lower()
        if any(module in lower for module in POPULATION_MODULES) and any(term in lower for term in FREQUENCY_TERMS):
            result.append(column)
    return result


def maximum_frequency(values) -> str:
    found = []
    for value in values:
        for token in re.split(r"[;,|]", str(value or "")):
            try:
                number = float(token)
            except ValueError:
                continue
            if 0.0 <= number <= 1.0:
                found.append(number)
    return repr(max(found)) if found else "."


def consequence_label(value) -> str:
    code = str(value or "")
    return SO_LABELS.get(code, code.lower())


def export_annotations(database: Path, output: Path, fallback_vcf: Path | None = None) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    allowed = None
    if fallback_vcf is not None:
        allowed = set()
        with open_text(fallback_vcf, "rt") as source:
            for line in source:
                if line.startswith("#"):
                    continue
                chromosome, position, _, reference, alternate = line.rstrip("\n").split("\t")[:5]
                allowed.add(clean_variant(chromosome, position, reference, alternate))
    with sqlite3.connect(database) as connection:
        tables = {str(row[0]) for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "variant" not in tables:
            raise ValueError("OpenCRAVAT database has no variant table")
        columns = table_columns(connection, "variant")
        chrom = first_present(columns, "original_input__chrom", "base__chrom", "chrom")
        # OpenCRAVAT removes the VCF anchor base from indels in its internal
        # representation.  extra_vcf_info preserves the original hg19 VCF key.
        pos = first_present(columns, "extra_vcf_info__pos", "original_input__pos", "base__pos", "pos")
        ref = first_present(columns, "extra_vcf_info__ref", "original_input__ref_base", "base__ref_base", "ref_base")
        alt = first_present(columns, "extra_vcf_info__alt", "original_input__alt_base", "base__alt_base", "alt_base")
        gene = first_present(columns, "base__hugo", "hugo")
        consequence = first_present(columns, "base__so", "so")
        af_columns = frequency_columns(columns)
        selected = [chrom, pos, ref, alt, gene, consequence, *af_columns]
        quoted = ", ".join(f'"{name}"' for name in selected)
        database_rows = total = fallback = skipped_database_rows = 0
        exported = set()
        with open_text(output, "wt") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(OUTPUT_COLUMNS)
            for row in connection.execute(f'SELECT {quoted} FROM "variant"'):
                database_rows += 1
                chromosome, position, reference, alternate, hugo, so, *frequencies = row
                clean = clean_variant(chromosome, position, reference, alternate)
                if clean in exported or (allowed is not None and clean not in allowed):
                    skipped_database_rows += 1
                    continue
                writer.writerow([*clean, hugo or "UNKNOWN", consequence_label(so), maximum_frequency(frequencies)])
                exported.add(clean)
                total += 1
            if fallback_vcf is not None:
                with open_text(fallback_vcf, "rt") as source:
                    for line in source:
                        if line.startswith("#"):
                            continue
                        chromosome, position, _, reference, alternate = line.rstrip("\n").split("\t")[:5]
                        clean = clean_variant(chromosome, position, reference, alternate)
                        if clean in exported:
                            continue
                        writer.writerow([*clean, "UNKNOWN", "missense_variant", "."])
                        fallback += 1
    return {
        "annotated_variants": total,
        "fallback_variants": fallback,
        "output_variants": total + fallback,
        "database_rows": database_rows,
        "skipped_duplicate_or_noncanonical_database_rows": skipped_database_rows,
        "population_frequency_columns": af_columns,
        "clinvar_features_exported": False,
        "output": str(output),
    }


def export_declared_missense(vcf: Path, output: Path) -> dict:
    """Create a no-annotation ablation using only VarPB's missense contract.

    VarPB states that the background table contains missense variants.  This
    export therefore marks every union row as missense, while leaving gene and
    population frequency unavailable.  It is useful as a strict lower-information
    test, but it must not be reported as the fully annotated model result.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with open_text(vcf, "rt") as source, open_text(output, "wt") as destination:
        writer = csv.writer(destination, delimiter="\t")
        writer.writerow(OUTPUT_COLUMNS)
        for line in source:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt = line.rstrip("\n").split("\t")[:5]
            clean = clean_variant(chrom, pos, ref, alt)
            writer.writerow([*clean, "UNKNOWN", "missense_variant", "."])
            total += 1
    return {
        "variants": total,
        "annotation_level": "VarPB-declared missense ablation",
        "gene_available": False,
        "population_frequency_available": False,
        "output": str(output),
    }


def export_predictor_score(database: Path, vcf: Path, column: str, score_name: str, output: Path) -> dict:
    """Export one OpenCRAVAT predictor in VarPB's five-column new-tool format."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(database) as connection:
        columns = table_columns(connection, "variant")
        if column not in columns:
            raise ValueError(f"OpenCRAVAT variant table has no {column!r} column")
        chrom = first_present(columns, "original_input__chrom", "base__chrom")
        pos = first_present(columns, "extra_vcf_info__pos", "original_input__pos", "base__pos")
        ref = first_present(columns, "extra_vcf_info__ref", "original_input__ref_base", "base__ref_base")
        alt = first_present(columns, "extra_vcf_info__alt", "original_input__alt_base", "base__alt_base")
        selected = [chrom, pos, ref, alt, column]
        quoted = ", ".join(f'"{name}"' for name in selected)
        scores = {}
        for chromosome, position, reference, alternate, value in connection.execute(f'SELECT {quoted} FROM "variant"'):
            try:
                variant = clean_variant(chromosome, position, reference, alternate)
            except ValueError:
                continue
            if variant not in scores and value is not None:
                scores[variant] = float(value)
    observed = missing = 0
    with open_text(vcf, "rt") as source, output.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.writer(destination, delimiter="\t")
        writer.writerow(["Chromosome", "Position", "Ref", "Alt", score_name])
        for line in source:
            if line.startswith("#"):
                continue
            chromosome, position, _, reference, alternate = line.rstrip("\n").split("\t")[:5]
            variant = clean_variant(chromosome, position, reference, alternate)
            value = scores.get(variant)
            if value is None:
                value = 0.0
                missing += 1
            else:
                observed += 1
            writer.writerow([f"chr{variant[0]}", variant[1], variant[2], variant[3], value])
    return {
        "source_column": column,
        "score_name": score_name,
        "observed_scores": observed,
        "missing_scores_filled_with_zero": missing,
        "output": str(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare-vcf")
    prepare.add_argument("--background", type=Path, required=True)
    prepare.add_argument("--causal", type=Path, nargs="+", required=True)
    prepare.add_argument("--output", type=Path, required=True)
    prepare.add_argument("--index", type=Path, required=True)
    export = subparsers.add_parser("export-annotations")
    export.add_argument("--database", type=Path, required=True)
    export.add_argument("--output", type=Path, required=True)
    export.add_argument("--fallback-vcf", type=Path)
    ablation = subparsers.add_parser("export-declared-missense")
    ablation.add_argument("--vcf", type=Path, required=True)
    ablation.add_argument("--output", type=Path, required=True)
    predictor = subparsers.add_parser("export-predictor-score")
    predictor.add_argument("--database", type=Path, required=True)
    predictor.add_argument("--vcf", type=Path, required=True)
    predictor.add_argument("--column", required=True)
    predictor.add_argument("--score-name", required=True)
    predictor.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare-vcf":
        result = prepare_vcf(args.background, args.causal, args.output, args.index)
    elif args.command == "export-annotations":
        result = export_annotations(args.database, args.output, args.fallback_vcf)
    elif args.command == "export-declared-missense":
        result = export_declared_missense(args.vcf, args.output)
    else:
        result = export_predictor_score(args.database, args.vcf, args.column, args.score_name, args.output)
    print(result)


if __name__ == "__main__":
    main()
