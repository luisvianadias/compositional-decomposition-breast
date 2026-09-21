# -*- coding: utf-8 -*-
"""Figure 2 — Capture by intrinsic subtype (METABRIC, n=1,980).
Single classifier (logistic regression). Clean, no model comparison."""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "results")
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "figures")
os.makedirs(FIG, exist_ok=True)

C_ORANGE = "#E69F00"

C_BAR = "#E69F00"
C_TEXT = "#B36B00"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 8,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
})

data_path = os.path.join(RES, "08_clinical_subtypes.json")
if not os.path.exists(data_path):
    sys.exit(f"[figure2] {data_path} ausente — rode 08_clinical_subtypes.py "
             "antes.")
cap = json.load(open(data_path))["pam50_capture"]

order = ["Basal", "Her2", "LumB", "LumA", "Normal"]
labels, n_str, vals = [], [], []

for st in order:
    v = cap.get(st, cap.get(st.replace("Basal", "Basal-like"), None))
    if v is None:
        continue
    labels.append(st)
    vals.append(v["mean"] * 100)
    n_str.append(f"n={int(v['n']):,}")

fig, ax = plt.subplots(figsize=(5.5, 3.4))
x = np.arange(len(vals))

ax.yaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)

bars = ax.bar(x, vals, 0.48, color=C_ORANGE, alpha=0.9,
              edgecolor="none", zorder=3)

for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 2.0, f"{v:.0f}%",
            ha="center", va="bottom", fontsize=7.5, fontweight="bold",
            color="#B36B00")

ax.axhline(50, color="#D95F02", ls="--", lw=0.9, alpha=0.6, zorder=2,
           label="50% (chance)")

ax.set_xticks(x)
ax.set_xticklabels(
    [f"{l}\n({n})" for l, n in zip(labels, n_str)],
    fontsize=7.8, color="#222222")
ax.set_ylabel("Capture (%)", fontsize=8.5, labelpad=6)
ax.set_ylim(0, 112)
ax.set_yticks([0, 20, 40, 60, 80, 100])

ax.legend(loc="upper right", fontsize=7.5, frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(FIG, "figure2_pam50.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "figure2_pam50.png"), bbox_inches="tight",
            dpi=300)
plt.close(fig)
print("Figure 2 salva — barras LogReg apenas, sem comparação de modelos")
