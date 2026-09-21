# -*- coding: utf-8 -*-
"""Gera Figura 1 e Figura S1 com proporções e contenção de texto rigorosas.
(versão do usuário, lendo results/ público — nested_compression.json +
05_cross_cohort_axis.json)"""

import json
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
from matplotlib.lines import Line2D

matplotlib.use("Agg")
sys.stdout.reconfigure(encoding="utf-8")

DIR = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(DIR, "..", "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "font.size": 8,
        "axes.labelsize": 8.5,
        "axes.titlesize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    }
)

C_NULL = "#999999"
C_BLUE = "#0072B2"
C_ORANGE = "#D55E00"
C_GREEN = "#009E73"
C_YELLOW = "#E69F00"
C_GRID = "#F0F0F0"

nested = json.load(open(os.path.join(DIR, "..", "results",
                                     "nested_compression.json")))
axis = json.load(open(os.path.join(DIR, "..", "results",
                                   "05_cross_cohort_axis.json")))
ablation = json.load(open(os.path.join(DIR, "..", "results",
                                       "04_adipose_ablation.json")))

K = [1, 2, 3, 5, 10, 20]
curve = [nested[f"k{k}"]["acc"] * 100 for k in K]
full_acc = nested["full"]["acc"] * 100
panel_accs = {n: nested[n]["acc"] * 100
              for n in ("panelA", "panelB", "panelC", "randDE")}
resid_acc = ablation["after"]["acc"] * 100
resid_auc = ablation["after"]["auc"]
null_draws = np.array(nested["null_acc_per_draw"]) * 100
null_med = float(np.median(null_draws))
null_q25, null_q75 = [float(q) for q in
                      np.quantile(null_draws, [0.25, 0.75])]

ladder_arms = [
    ("Null\n(B=500)", null_draws),
    ("Panel A\n(12 adip.)",
     np.array(nested["panelA"]["per_split_acc"]) * 100),
    ("Panel B\n(15 ad.-ass.)",
     np.array(nested["panelB"]["per_split_acc"]) * 100),
    ("randDE\n(matched)",
     np.array(nested["randDE"]["per_split_acc"]) * 100),
    ("Panel C\n(+stroma)",
     np.array(nested["panelC"]["per_split_acc"]) * 100),
    ("Full\n(23k genes)",
     np.array(nested["full"]["per_split_acc"]) * 100),
]


def panel_tag(ax, letter):
    ax.text(-0.12, 1.08, letter, transform=ax.transAxes, fontsize=10.5,
            fontweight="bold", va="top", ha="right")


fig = plt.figure(figsize=(7.2, 5.2))
gs = fig.add_gridspec(
    2, 3, height_ratios=[1.0, 1.05], width_ratios=[0.85, 0.95, 1.2],
    hspace=0.48, wspace=0.38, left=0.08, right=0.98, top=0.94, bottom=0.08,
)

axA = fig.add_subplot(gs[0, 0])
panel_tag(axA, "A")
axA.grid(axis="y", color=C_GRID, lw=0.6, zorder=0)
axA.axhspan(null_q25, null_q75, color=C_NULL, alpha=0.25, lw=0, zorder=1)
axA.axhline(null_med, color=C_NULL, lw=0.9, ls="--", zorder=2)
axA.axhline(full_acc, color="#333333", lw=0.8, ls=":", zorder=2)
axA.plot(K, curve, "o-", color=C_BLUE, lw=1.4, ms=4.2, mec="white",
         mew=0.5, zorder=3)
axA.set_xscale("log")
axA.set_xticks(K)
axA.set_xticklabels([str(k) for k in K])
axA.set_xlabel("Genes in nested panel ($k$)")
axA.set_ylabel("Accuracy (%)")
axA.set_ylim(80, 101)
axA.text(0.05, 0.08,
         f"Full = {full_acc:.1f}%\nNull median = {null_med:.1f}%",
         transform=axA.transAxes, fontsize=6.8, color="#444444")

axB = fig.add_subplot(gs[0, 1:])
panel_tag(axB, "B")
axB.grid(axis="y", color=C_GRID, lw=0.6, zorder=0)
colors = [C_NULL, C_BLUE, C_BLUE, C_YELLOW, C_GREEN, "#333333"]
pos = np.arange(len(ladder_arms))
data = [a * 100 if a.max() <= 1.01 else a for _, a in ladder_arms]
bp = axB.boxplot(data, positions=pos, widths=0.52, patch_artist=True,
                 showfliers=False,
                 medianprops=dict(color="black", lw=1.1),
                 whiskerprops=dict(lw=0.7, color="#444444"),
                 capprops=dict(lw=0.7, color="#444444"), zorder=2)
for patch, c in zip(bp["boxes"], colors):
    patch.set_facecolor(c)
    patch.set_alpha(0.85 if c != C_NULL else 0.45)
    patch.set_edgecolor("#222222")
    patch.set_linewidth(0.7)
axB.set_xticks(pos)
axB.set_xticklabels([n for n, _ in ladder_arms], fontsize=6.8)
axB.set_ylabel("Accuracy (%)")
axB.set_ylim(78, 100.5)

axC = fig.add_subplot(gs[1, 0])
panel_tag(axC, "C")
axC.grid(axis="x", color=C_GRID, lw=0.6, zorder=0)
y_pos = [1, 0]
for y, val, col in zip(y_pos, [full_acc, resid_acc], [C_BLUE, C_YELLOW]):
    axC.hlines(y, xmin=82, xmax=val, color=col, lw=2.0, zorder=2)
    axC.plot(val, y, "o", color=col, ms=7, mec="black", mew=0.7, zorder=3)
    axC.text(val + 0.6, y, f"{val:.1f}%", va="center", ha="left",
             fontweight="bold", fontsize=7.2)
axC.set_yticks(y_pos)
axC.set_yticklabels(["Full model", "Residual model\n(non-adipose)"],
                    fontsize=7.2)
axC.set_xlim(82, 103)
axC.set_xlabel("Accuracy (%)")
axC.text(0.06, 0.14,
         f"AUC: {ablation['before']['auc']:.3f} -> {resid_auc:.3f}",
         transform=axC.transAxes, fontsize=6.8, color="#444444")

axD = fig.add_subplot(gs[1, 1])
panel_tag(axD, "D")
axD.grid(axis="y", color=C_GRID, lw=0.6, zorder=0)
pops = ["TCGA\ntumor", "TCGA\nnormal", "GTEx\nbreast", "GTEx\nadipose"]
valsD = [axis["tcga_tumor"], axis["tcga_adjacent_normal"],
         axis["gtex_breast"], axis["gtex_adipose"]]
colsD = [C_ORANGE, C_BLUE, C_GREEN, C_YELLOW]
bars = axD.bar(pops, valsD, color=colsD, width=0.55, edgecolor="#222222",
               lw=0.6, zorder=2, alpha=0.9)
axD.axhline(0, color="#333333", lw=0.7, zorder=3)
axD.set_ylabel("Adipose score ($z$)")
axD.set_ylim(-1.6, 0.95)
for b, val in zip(bars, valsD):
    offset = 0.05 if val >= 0 else -0.12
    va = "bottom" if val >= 0 else "top"
    axD.text(b.get_x() + b.get_width() / 2, val + offset, f"{val:+.2f}",
             ha="center", va=va, fontsize=6.5, fontweight="bold")

axE = fig.add_subplot(gs[1, 2])
panel_tag(axE, "E")
axE.set_xlim(0, 10)
axE.set_ylim(0, 10)
axE.axis("off")


def node_box(x, y, w, h, text, fc):
    axE.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                                 fc=fc, ec="#333333", lw=0.7, zorder=2))
    axE.text(x + w / 2, y + h / 2, text, ha="center", va="center",
             fontsize=6.8, zorder=3, linespacing=1.2)


def conn_arrow(start, end):
    axE.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                                  mutation_scale=8, lw=0.8,
                                  color="#333333", zorder=1))


node_box(0.8, 8.2, 8.4, 1.4,
         f"Full Transcriptome: {full_acc:.1f}% "
         f"(AUC {ablation['before']['auc']:.3f})", "#EAEAEA")
conn_arrow((3.0, 8.2), (2.5, 7.0))
conn_arrow((7.0, 8.2), (7.5, 7.0))
node_box(0.2, 4.8, 4.5, 2.2,
         f"Adipose-Stromal Axis\nPanel C: {panel_accs['panelC']:.1f}%\n"
         "- Adipocyte loss\n- Stroma exp.", "#D9EAF7")
node_box(5.3, 4.8, 4.5, 2.2,
         f"Residual Signal\n{resid_acc:.1f}% (AUC {resid_auc:.3f})\n"
         "- Tumor intrinsic\n- Non-adipose", "#FBE8D6")
conn_arrow((2.45, 4.8), (2.45, 3.2))
node_box(0.2, 1.6, 4.5, 1.6,
         "Tissue Architecture\nLobular < Ductal\n(p <= 1.7e-4)", "#EBF4FA")

fig.savefig(os.path.join(FIG, "figure1.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "figure1.png"), bbox_inches="tight", dpi=300)
plt.close(fig)

figS, axS = plt.subplots(figsize=(6.2, 3.5))
panel_tag(axS, "S1")
axS.grid(axis="y", color=C_GRID, lw=0.6, zorder=0)

n, bins, _ = axS.hist(
    null_draws,
    bins=30,
    color="#D8DCE3",
    edgecolor="#6B7280",
    lw=0.5,
    density=False,
    zorder=2,
)
max_h = max(n)
axS.set_ylim(0, max_h * 1.22)
axS.set_xlim(80, 101)

lines_data = [
    (null_med, C_NULL, "--", 1.1),
    (panel_accs["panelA"], C_BLUE, ":", 1.2),
    (panel_accs["panelB"], C_BLUE, "-.", 1.2),
    (panel_accs["randDE"], C_YELLOW, "--", 1.2),
    (panel_accs["panelC"], C_GREEN, "-", 1.4),
    (full_acc, "#111111", "-", 1.4),
]
for x, c, ls, lw in lines_data:
    axS.axvline(x, color=c, lw=lw, ls=ls, zorder=3)

legend_handles = [
    Line2D([0], [0], color=C_NULL, lw=1.2, ls="--",
           label=f"Null median ({null_med:.1f}%)"),
    Line2D([0], [0], color=C_BLUE, lw=1.2, ls=":",
           label=f"Panel A ({panel_accs['panelA']:.1f}%)"),
    Line2D([0], [0], color=C_BLUE, lw=1.2, ls="-.",
           label=f"Panel B ({panel_accs['panelB']:.1f}%)"),
    Line2D([0], [0], color=C_YELLOW, lw=1.2, ls="--",
           label=f"randDE ({panel_accs['randDE']:.1f}%)"),
    Line2D([0], [0], color=C_GREEN, lw=1.4, ls="-",
           label=f"Panel C ({panel_accs['panelC']:.1f}%)"),
    Line2D([0], [0], color="#111111", lw=1.4, ls="-",
           label=f"Full model ({full_acc:.1f}%)"),
]
axS.legend(
    handles=legend_handles, loc="upper left", frameon=True,
    facecolor="white", edgecolor="#D0D0D0", framealpha=0.95,
    fontsize=7.3, borderpad=0.6, labelspacing=0.45, handlelength=2.0,
)

null_max = float(null_draws.max())
n_ge = int((null_draws >= panel_accs["panelC"]).sum())
p_emp = (n_ge + 1) / (len(null_draws) + 1)
axS.annotate(
    f"Panel C & Full\n{n_ge}/{len(null_draws)} null >= "
    f"{panel_accs['panelC']:.1f}%\n"
    f"(empirical P = {p_emp:.3f}; null max = {null_max:.1f}%)",
    xy=(98.8, max_h * 0.75),
    xytext=(93.5, max_h * 0.95),
    fontsize=6.8, color="#222222", ha="center", va="center",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="#FAFAFA",
              edgecolor="#CCCCCC", lw=0.6),
    arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-0.15",
                    color="#444444", lw=0.8),
    zorder=4,
)

axS.set_xlabel("Accuracy (%) - 10-gene random panels (B = 500 draws)")
axS.set_ylabel("Count")
figS.tight_layout()
figS.savefig(os.path.join(FIG, "figureS1.pdf"), bbox_inches="tight")
figS.savefig(os.path.join(FIG, "figureS1.png"), bbox_inches="tight",
             dpi=300)
plt.close(figS)

print("Figuras 1 e S1 regeneradas com sucesso.")
