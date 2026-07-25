# make_um_consistency_fig.py -- clean vector figure for the filter-consistency
# robustness sweep, from um_consistency_results.mat.
import numpy as np
import scipy.io as sio
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "font.size": 8,
    "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.spines.top": False, "axes.linewidth": 0.8, "legend.frameon": False,
})
VERM = "#D55E00"; BLUE = "#0072B2"; GREY = "#999999"; BLACK = "#000000"

d = sio.loadmat(ROOT / "um_consistency_results.mat", squeeze_me=True, struct_as_record=False)
sw, ref = d["sweep"], d["ref"]
x = np.asarray(sw.scale, float)
coll = np.asarray(sw.coll, float); lo = np.asarray(sw.lo, float); hi = np.asarray(sw.hi, float)
sig = np.asarray(sw.sigma, float); le = np.asarray(sw.locerr, float)
ref_coll = float(ref.coll)

fig, axL = plt.subplots(figsize=(4.8, 3.1))

# shade overconfident region (filt/true < 1)
axL.axvspan(x.min()*0.8, 1.0, color="#f2e6e0", alpha=0.5, zorder=0)
axL.axvline(1.0, color=GREY, ls=":", lw=1.0, zorder=1)

# left axis: collision rate with asymmetric Wilson CI
axL.errorbar(x, coll, yerr=[coll - lo, hi - coll], fmt="o-", color=BLUE,
             lw=1.6, ms=5, mfc=BLUE, capsize=3, elinewidth=1.0, zorder=3,
             label="cov-aware collisions (95 % CI)")
axL.axhline(ref_coll, color=BLUE, ls="--", lw=1.1, zorder=2,
            label=f"fixed-margin ref ({ref_coll:.0f} %)")
axL.set_ylabel("Collision rate (%)", color=BLUE)
axL.tick_params(axis="y", labelcolor=BLUE)
axL.set_ylim(-1, max(hi.max(), ref_coll) * 1.12)
axL.set_xscale("log")
axL.set_xlabel("filter-assumed / true noise variance  (overconfident ←  1  → conservative)")
axL.set_xticks(x)
axL.set_xticklabels([f"{v:g}" for v in x])
axL.grid(True, which="major", axis="y", lw=0.3, alpha=0.4)

# right axis: reported sigma vs true localization error
axR = axL.twinx()
axR.spines["top"].set_visible(False)
axR.plot(x, sig, "s-", color=VERM, lw=1.4, ms=4, mfc=VERM, label="reported $\\sigma$")
axR.plot(x, le, "^:", color=BLACK, lw=1.4, ms=4, mfc="none", label="true loc. error")
axR.set_ylabel("reported $\\sigma$ / true loc. error (m)")

# combined legend
h1, l1 = axL.get_legend_handles_labels()
h2, l2 = axR.get_legend_handles_labels()
axL.legend(h1 + h2, l1 + l2, loc="upper center", ncol=2, fontsize=6.3,
           handlelength=1.6, columnspacing=1.2)
axL.set_title("Robustness to filter (in)consistency", fontsize=8.5)
fig.tight_layout()
fig.savefig(OUT / "um_consistency.pdf", bbox_inches="tight")
fig.savefig(OUT / "um_consistency_clean.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("saved um_consistency.pdf / um_consistency_clean.png")
