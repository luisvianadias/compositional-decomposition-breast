# -*- coding: utf-8 -*-
"""Junta todos os .py da pasta em scripts.md, separados por arquivo."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "scripts.md"
SELF = Path(__file__).name

files = sorted(
    p for p in HERE.glob("*.py")
    if p.name != SELF and p.is_file()
)

with OUT.open("w", encoding="utf-8") as out:
    out.write("# Scripts\n\n")
    out.write(f"Transcrição de {len(files)} arquivos `.py` de `scripts/`.\n\n")

    for i, p in enumerate(files):
        code = p.read_text(encoding="utf-8")
        out.write("---\n\n")
        out.write(f"## {p.name}\n\n")
        out.write("```python\n")
        out.write(code.rstrip("\n"))
        out.write("\n```\n\n")

print(f"[ok] {len(files)} scripts -> {OUT}")
for p in files:
    print(f"  - {p.name}")
