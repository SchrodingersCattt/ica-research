"""Fast, non-submitting validation for the ICA research review bundle."""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIVATE_MARKERS = re.compile(
    "|".join(
        [
            "BEGIN" + r" .* " + "PRIVATE" + r" " + "KEY",
            "pass" + r"word\s*[:=]",
            "api" + r"[_-]?" + "key" + r"\s*[:=]",
            "sec" + r"ret\s*[:=]",
        ]
    ),
    re.IGNORECASE,
)


def check_python() -> list[str]:
    failures = []
    for path in (ROOT / "workflows").rglob("*.py"):
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError as exc:
            failures.append(f"syntax: {path}: {exc}")
    return failures


def main() -> int:
    failures = check_python()
    structures = sorted((ROOT / "examples" / "structures").glob("*.cif"))
    if len(structures) < 3:
        failures.append("expected CA, LA and NEt4-CA example structures")
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or ".git" in path.parts:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        if path.suffix.lower() not in {".py", ".md", ".yml", ".yaml", ".json", ".tcl", ".bat", ".lmp", ".inp"}:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if PRIVATE_MARKERS.search(text):
            failures.append(f"private marker: {path}")
    template = ROOT / "workflows" / "finetune" / "template.example.json"
    try:
        data = json.loads(template.read_text(encoding="utf-8"))
        assert data["model"]["descriptor"]["repflow"]["n_dim"] == 256
    except Exception as exc:  # pragma: no cover - smoke-test diagnostic
        failures.append(f"training template: {exc}")
    if failures:
        print("SMOKE TEST FAILED")
        print("\n".join(failures))
        return 1
    print(f"SMOKE TEST PASSED: {len(structures)} structures; Python syntax and redaction checks OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
