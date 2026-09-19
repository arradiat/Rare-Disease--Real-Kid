"""Unified command-line interface for the project."""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections.abc import Callable

from rare_disease_real_kid.paths import PROJECT_ROOT


def _invoke(main_function: Callable[[], object], forwarded: list[str]) -> int:
    previous = sys.argv
    try:
        sys.argv = [previous[0], *forwarded]
        result = main_function()
        return int(result) if isinstance(result, int) else 0
    finally:
        sys.argv = previous


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rdkid", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    doctor = commands.add_parser("doctor", help="Check an execution environment")
    doctor.add_argument("--profile", choices=("core", "report", "wgs"), default="core")
    commands.add_parser("verify", help="Verify committed code and artifacts")

    submission = commands.add_parser("submission", help="Generate the Track 1 submission")
    submission.add_argument("action", choices=("generate",))
    submission.add_argument("args", nargs=argparse.REMAINDER)

    pipeline = commands.add_parser("pipeline", help="Run a WGS pipeline stage")
    pipeline.add_argument(
        "stage",
        nargs="?",
        choices=(
            "inspect",
            "phenotype",
            "prepare-reference",
            "align",
            "normalize",
            "phase",
            "annotate",
            "rank",
            "submit",
            "report",
            "all",
        ),
    )
    pipeline.add_argument("args", nargs=argparse.REMAINDER)

    benchmark = commands.add_parser("benchmark", help="Run a public-data benchmark")
    benchmark.add_argument("benchmark", choices=("clinvar", "phenopacket", "varpb"))
    benchmark.add_argument("args", nargs=argparse.REMAINDER)

    report = commands.add_parser("report", help="Build the methods report")
    report.add_argument("action", choices=("build",))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "doctor":
        from rare_disease_real_kid.doctor import main as doctor_main

        return doctor_main(["--profile", args.profile])
    if args.command == "verify":
        from rare_disease_real_kid.verification import main as verify_main

        return verify_main()
    if args.command == "submission":
        from rare_disease_real_kid.evaluation.submission_generate import main as submission_main

        return _invoke(submission_main, args.args)
    if args.command == "pipeline":
        from rare_disease_real_kid.pipeline.wgs import main as pipeline_main

        forwarded = [args.stage, *args.args] if args.stage else ["--help"]
        return _invoke(pipeline_main, forwarded)
    if args.command == "benchmark":
        if args.benchmark == "clinvar":
            from rare_disease_real_kid.experiments.clinvar_temporal import main as benchmark_main
        elif args.benchmark == "phenopacket":
            from rare_disease_real_kid.experiments.phenopacket_hpo import main as benchmark_main
        else:
            from rare_disease_real_kid.evaluation.varpb_suite import main as benchmark_main
        return _invoke(benchmark_main, args.args)
    if args.command == "report":
        completed = subprocess.run(
            ["bash", str(PROJECT_ROOT / "scripts" / "build_report.sh")],
            cwd=PROJECT_ROOT,
            check=False,
        )
        return completed.returncode
    return 2
