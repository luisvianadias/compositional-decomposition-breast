# -*- coding: utf-8 -*-
"""Gera requirements.txt com versões pinadas do ambiente atual."""
import importlib.metadata as im

packages = [
    "scikit-learn", "numpy", "scipy", "pandas",
    "requests", "matplotlib", "rdata",
]
lines = []
for pkg in packages:
    try:
        v = im.version(pkg)
        lines.append(f"{pkg}=={v}")
    except im.PackageNotFoundError:
        print(f"  aviso: {pkg} não instalado")
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  "requirements.txt"), "w").write("\n".join(lines) + "\n")
print("requirements.txt com versões pinadas")
