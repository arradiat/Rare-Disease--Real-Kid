"""Environment preflight for lightweight, reporting, or full-WGS workflows."""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import sys

from rare_disease_real_kid.paths import DATA, PROJECT_ROOT


PROFILES = {
    "core": {
        "modules": ("numpy", "scipy", "sklearn", "pandas"),
        "commands": (),
    },
    "report": {
        "modules": ("numpy", "scipy", "sklearn", "pandas", "matplotlib", "PIL"),
        "commands": ("pdflatex", "bibtex"),
    },
    "wgs": {
        "modules": ("pysam", "pybwa", "whatshap", "open_cravat", "docx", "openpyxl"),
        "commands": ("oc", "samtools", "whatshap"),
    },
}


def diagnose(profile: str = "core") -> dict[str, object]:
    requested = PROFILES[profile]
    modules = {name: importlib.util.find_spec(name) is not None for name in requested["modules"]}
    commands = {name: shutil.which(name) for name in requested["commands"]}
    missing = [name for name, present in modules.items() if not present]
    missing.extend(name for name, location in commands.items() if location is None)
    return {
        "profile": profile,
        "ok": not missing,
        "python": sys.version.split()[0],
        "project_root": str(PROJECT_ROOT),
        "data_dir": str(DATA),
        "data_dir_exists": DATA.exists(),
        "python_modules": modules,
        "commands": commands,
        "missing": missing,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=tuple(PROFILES), default="core")
    args = parser.parse_args(argv)
    result = diagnose(args.profile)
    print(json.dumps(result, indent=2))
    return 0 if result["ok"] else 1

