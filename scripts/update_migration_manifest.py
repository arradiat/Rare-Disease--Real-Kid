"""Add current destination hashes to the original-file migration ledger."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "artifacts/manifests/original_file_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    for record in payload["files"]:
        destination = ROOT / record["destination_path"]
        record["final_sha256"] = sha256(destination)
        record["final_size_bytes"] = destination.stat().st_size
        record["content_status"] = (
            "unchanged" if record["final_sha256"] == record["sha256"] else "modified_after_migration"
        )
    payload["finalized"] = True
    MANIFEST.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Finalized {len(payload['files'])} migration records")


if __name__ == "__main__":
    main()

