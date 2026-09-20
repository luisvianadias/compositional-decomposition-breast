# -*- coding: utf-8 -*-
import os
import importlib.metadata as im

packages = ["scikit-learn", "numpy", "scipy", "pandas",
            "requests", "matplotlib", "rdata"]
lines = []
for pkg in packages:
    try:
        v = im.version(pkg)
        lines.append(f"{pkg}=={v}")
    except im.PackageNotFoundError:
        pass
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "requirements.txt")
open(out, "w").write("\n".join(lines) + "\n")
print("requirements.txt pinado")
