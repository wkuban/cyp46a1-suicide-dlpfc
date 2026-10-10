"""Figures 2-5 from the rebuilt RT-qPCR analysis.

Figure 2  CYP46A1 and CYP2D6 by group.
Figure 3  CYP19A1 and the redox partners POR and FDX1.
Figure 4  Summary of the five analysed estimates with 95% confidence intervals.
Figure 5  Pathway schematic.

CYP11A1 and CYP17A1 are absent from Figures 3 and 4 and marked as not
quantifiable in Figure 5: their signal sits at or beyond the limit of
quantification (median Ct 35.9 and 38.4, with 27% and 64% of wells returning no
Ct), so they do not enter the group comparison or the multiplicity correction.

Inputs : results/qpcr_rebuilt_relative.json, results/qpcr_statistics_rebuilt.json
Outputs: figures/Figure_2 .. Figure_5 (.tiff/.png/.pdf)
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy import stats

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_FIG = os.path.join(_ROOT, "figures")
os.makedirs(_FIG, exist_ok=True)
_res = lambda f: os.path.join(_ROOT, "results", f)
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 9,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8,
                     "ytick.major.width": 0.8, "pdf.fonttype": 42})

D = json.load(open(_res("qpcr_rebuilt_relative.json")))
R = json.load(open(_res("qpcr_statistics_rebuilt.json")))
CC, SC, INK, MUTED = "#2F6DB5", "#C7782A", "#222222", "#6b6b6b"
NAMES = {"FDX1": "FDX1 (adrenodoxin)"}
lab2 = lambda t: ("%g" % (2.0 ** t)) if t >= -3 else "1/%d" % (2 ** -t)


def save(fig, name):
    fig.savefig(os.path.join(_FIG, name + ".tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(os.path.join(_FIG, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(_FIG, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def panel(ax, gene, letter):
    c = np.log2(D[gene]["raw_ctrl"]); s = np.log2(D[gene]["raw_sui"])
    rng = np.random.default_rng(7)
    for i, (x, col) in enumerate([(c, CC), (s, SC)]):
        ax.scatter(np.full(len(x), i) + rng.uniform(-0.14, 0.14, len(x)), x, s=24,
                   facecolor=col, edgecolor="white", linewidth=0.6, zorder=3)
        m = x.mean()
        h = stats.t.ppf(.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
        ax.hlines(m, i - 0.27, i + 0.27, color=INK, lw=1.8, zorder=4)
        ax.errorbar(i, m, yerr=h, color=INK, capsize=4, lw=1, capthick=1, zorder=4)
    allv = np.concatenate([c, s])
    lo, hi = np.floor(allv.min()) - 0.5, np.ceil(allv.max()) + 1.7
    ax.set_ylim(lo, hi)
    ticks = [t for t in range(int(np.floor(lo)), int(np.ceil(hi)) + 1)]
    step = 2 if len(ticks) > 9 else 1
    ticks = [t for t in ticks if t % step == 0 and lo <= t <= hi]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.set_yticklabels([lab2(t) for t in ticks])
    ax.axhline(0, color="#9a9a9a", lw=0.6, ls=(0, (3, 3)), zorder=1)
    ax.set_xlim(-0.6, 1.6); ax.set_xticks([0, 1])
    ax.set_xticklabels(["Control\n(n = %d)" % len(c), "Suicide\n(n = %d)" % len(s)])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(NAMES.get(gene, gene), fontsize=10, fontweight="bold")
    ax.set_ylabel("Relative mRNA expression\n(2$^{\\mathregular{-ΔΔCt}}$, "
                  "log$_{\\mathregular{2}}$ scale)")
    r = R[gene]; ytop = allv.max() + 0.75
    ax.plot([0, 0, 1, 1], [ytop - 0.2, ytop, ytop, ytop - 0.2], color=INK, lw=0.8)
    lab = ("p = %.3f" % r["p"]) if r["p"] >= 0.001 else "p < 0.001"
    ax.text(0.5, ytop + 0.12, lab + ("*" if r["p"] < 0.05 else ""), ha="center",
            va="bottom", fontsize=8.5, color=INK)
    ax.text(-0.34, 1.07, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom")


fig, axs = plt.subplots(1, 2, figsize=(5.2, 3.1))
for ax, g, l in zip(axs, ["CYP46A1", "CYP2D6"], "AB"):
    panel(ax, g, l)
fig.tight_layout(w_pad=1.8); save(fig, "Figure_2")

fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.2))
for ax, g, l in zip(axs, ["CYP19A1", "POR", "FDX1"], "ABC"):
    panel(ax, g, l)
    if g != "CYP19A1":
        ax.set_ylabel("")
fig.tight_layout(w_pad=1.6); save(fig, "Figure_3")

ORDER = ["CYP46A1", "CYP2D6", "CYP19A1", "POR", "FDX1"]
fig, ax = plt.subplots(figsize=(6.6, 2.6))
y = np.arange(len(ORDER))[::-1]
for yi, g in zip(y, ORDER):
    r = R[g]; col = SC if r["p"] < 0.05 else "#7a7a7a"
    ax.plot([np.log2(r["fold_change_lo"]), np.log2(r["fold_change_hi"])], [yi, yi],
            color=col, lw=1.6, solid_capstyle="round")
    ax.plot(np.log2(r["fold_change"]), yi, "o", ms=7, color=col, mec="white",
            mew=0.8, zorder=3)
ax.axvline(0, color="#9a9a9a", lw=0.8, ls=(0, (3, 3)))
ax.set_yticks(y); ax.set_yticklabels(ORDER, fontstyle="italic")
xt = [-4, -3, -2, -1, 0, 1, 2, 3]
ax.set_xticks(xt); ax.set_xticklabels([lab2(t) for t in xt]); ax.set_xlim(-4.3, 3.2)
ax.set_xlabel("Fold change, suicide vs control "
              "(95% CI, log$_{\\mathregular{2}}$ scale)")
ax.spines[["top", "right", "left"]].set_visible(False)
ax.tick_params(axis="y", length=0)
for xpos, head in ((1.02, "n (C/S)"), (1.17, "p"), (1.29, "q (FDR)")):
    ax.text(xpos, 1.06, head, transform=ax.transAxes, fontsize=8,
            fontweight="bold", ha="left")
for yi, g in zip(y, ORDER):
    r = R[g]
    ax.annotate("%d/%d" % (r["nc"], r["ns"]), xy=(1.02, yi),
                xycoords=("axes fraction", "data"), fontsize=8, va="center")
    ax.annotate("%.3f" % r["p"], xy=(1.17, yi), xycoords=("axes fraction", "data"),
                fontsize=8, va="center",
                fontweight="bold" if r["p"] < 0.05 else "normal")
    ax.annotate("%.3f" % r["q_fdr"], xy=(1.29, yi),
                xycoords=("axes fraction", "data"), fontsize=8, va="center")
ax.annotate("lower in suicide", xy=(-0.2, -0.95), fontsize=7.5, color=MUTED, ha="right")
ax.annotate("higher in suicide", xy=(0.2, -0.95), fontsize=7.5, color=MUTED, ha="left")
ax.set_ylim(-1.2, len(ORDER) - 0.5)
fig.tight_layout(); save(fig, "Figure_4")

# Figure 5: pathway schematic, numbers taken from the rebuilt statistics
fig, ax = plt.subplots(figsize=(8.6, 4.7))
ax.set_xlim(0, 120); ax.set_ylim(0, 62); ax.axis("off")
DOWN, NS, ND = SC, "#8a8a8a", "#bdbdbd"
fc = lambda g: "FC %.2f (%.2f–%.2f), p = %.3f" % (
    R[g]["fold_change"], R[g]["fold_change_lo"], R[g]["fold_change_hi"], R[g]["p"])


def box(x, y, w, h, text, fcol="white", ec="#444444", ls="-", fs=7.6, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.3,rounding_size=1.0",
                                fc=fcol, ec=ec, lw=1, ls=ls))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color=INK)


def arrow(x1, y1, x2, y2, col="#444444"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=9, lw=1.1, color=col,
                                 shrinkA=0, shrinkB=0))


ax.text(24, 59, "Monoamines & endocannabinoids", ha="center", fontsize=9,
        fontweight="bold", color=INK)
ax.text(64, 59, "Cholesterol elimination", ha="center", fontsize=9,
        fontweight="bold", color=INK)
ax.text(99, 59, "Neurosteroidogenesis", ha="center", fontsize=9,
        fontweight="bold", color=INK)
for i, (a_, b_) in enumerate([("Tyramine", "Dopamine"),
                              ("5-Methoxytryptamine", "Serotonin"),
                              ("Anandamide", "Oxidized metabolites")]):
    yy = 49 - i * 7
    box(1, yy, 20, 4, a_); box(27, yy, 20, 4, b_); arrow(21.6, yy + 2, 26.4, yy + 2, DOWN)
box(17, 24, 14, 4.2, "CYP2D6", fcol="#f6e3cf", ec=DOWN, fs=8.5, bold=True)
ax.text(24, 21.6, "↓ " + fc("CYP2D6"), ha="center", va="center", fontsize=7, color=DOWN)
ax.plot([24, 24], [28.6, 34.6], color=DOWN, lw=1, ls=(0, (2, 2)))
box(53, 49, 22, 4, "Cholesterol")
box(53, 29, 22, 4, "24S-Hydroxycholesterol")
box(53, 18.5, 22, 4, "Efflux across the BBB", ec="#9a9a9a")
arrow(64, 48.4, 64, 33.6, DOWN); arrow(64, 28.4, 64, 23.1)
box(66, 42, 13, 4.2, "CYP46A1", fcol="#f6e3cf", ec=DOWN, fs=8.5, bold=True)
ax.text(74.5, 40.2, "↓ " + fc("CYP46A1").replace(", p", "\np"), ha="center",
        va="top", fontsize=7, color=DOWN)
ax.text(64, 15.2, "CYP46A1 is neuron-specific", ha="center", fontsize=7,
        color=MUTED, style="italic")
for yy, t in [(49, "Cholesterol"), (39, "Pregnenolone"),
              (29, "DHEA / androgens"), (19, "Estrogens")]:
    box(90, yy, 18, 4, t)
for y1, n, note in [(48.4, "CYP11A1", "not quantifiable"),
                    (38.4, "CYP17A1", "not quantifiable"),
                    (28.4, "CYP19A1", "n.s.")]:
    arrow(99, y1, 99, y1 - 5.0)
    ax.text(100.6, y1 - 2.5, n, fontsize=7.4, color=MUTED, va="center",
            fontweight="bold")
    ax.text(100.6, y1 - 4.4, note, fontsize=6.6, color=MUTED, va="center")
ax.text(97.2, 44.6, "FDX1 n.s.", fontsize=7, color=MUTED, ha="right", va="center")
box(85.5, 7.5, 27, 4.6, "CYP11B1 · CYP11B2 · CYP21A2", ec=ND, ls="--", fs=6.8)
ax.text(99, 4.6, "not detected in DLPFC", ha="center", fontsize=7, color=MUTED)
box(1, 8.6, 82, 4.4, "POR (P450 oxidoreductase), electron donor for CYP2D6, "
    "CYP46A1, CYP17A1 and CYP19A1: n.s.", fcol="#f4f4f4", ec="#9a9a9a", fs=7.0)
for x, fcol, ec, ls, text in ((1, "#f6e3cf", DOWN, "-", "lower in suicide (nominal; q > 0.05)"),
                              (38, "#f4f4f4", NS, "-", "no difference (n.s.)"),
                              (62, "white", ND, "--", "not detected / not quantifiable")):
    ax.add_patch(FancyBboxPatch((x, 0.8), 3, 2.2, boxstyle="round,pad=0.2",
                                fc=fcol, ec=ec, ls=ls))
    ax.text(x + 4.5, 1.9, text, fontsize=7.2, va="center")
save(fig, "Figure_5")
print("zapisano Figure_2 .. Figure_5 w", _FIG)
