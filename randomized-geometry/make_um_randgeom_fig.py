# make_um_randgeom_fig.py -- randomized-geometry transfer result.
# (a) sample random geometries, (b) pooled collision rates (Wilson CI),
# (c) per-geometry: fixed margin collides everywhere, cov-aware nowhere.
import numpy as np
import scipy.io as sio
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
mpl.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "pdf.fonttype":42,"font.size":8,"axes.labelsize":8,"axes.titlesize":8.5,
    "legend.fontsize":7,"xtick.labelsize":7.5,"ytick.labelsize":7.5,
    "axes.spines.top":False,"axes.spines.right":False,
    "axes.linewidth":0.8,"legend.frameon":False})
BLUE="#0072B2"; VERM="#D55E00"; GREEN="#009E73"; GREY="#999999"; PURPLE="#CC79A7"; BLACK="#000000"

d = sio.loadmat(ROOT/"um_randgeom_results.mat", squeeze_me=True, struct_as_record=False)
geom = np.atleast_1d(d["geom"]); COLL = np.atleast_2d(d["COLL"]); modes = [str(m) for m in d["modes"]]
NG = int(d["NG"]); M = int(d["M"]); base = d["base"]
wps = np.atleast_2d(base.wps); start = np.ravel(base.x0_nom)[:2]
path = np.column_stack([start[:,None], wps, start[:,None]])

def wil(k,n):
    z=1.96;p=k/n;den=1+z**2/n;ctr=(p+z**2/(2*n))/den
    hw=z*np.sqrt(p*(1-p)/n+z**2/(4*n**2))/den
    return 100*max(0,ctr-hw),100*min(1,ctr+hw)

fig, ax = plt.subplots(1, 3, figsize=(11.0, 3.3))

# (a) sample geometries
cols = [BLUE, VERM, GREEN, PURPLE, "#E69F00", "#000000"]
ax[0].plot(path[0], path[1], "-", color=GREY, lw=1.2, alpha=0.7, zorder=1, label="patrol path")
ax[0].plot(wps[0], wps[1], "*", ms=12, mfc=GREEN, mec="k", mew=0.5, ls="none", zorder=4)
ax[0].plot(start[0], start[1], "o", ms=6, mfc="k", mec="k", ls="none", zorder=4)
idx = np.linspace(0, NG-1, 6).astype(int)
for j,gi in enumerate(idx):
    g = geom[gi]; lm = np.atleast_2d(g.lm); ob = np.ravel(g.obs)
    c = cols[j]
    ax[0].add_patch(Circle((ob[0],ob[1]), ob[2], facecolor=c, edgecolor="k", lw=0.4, alpha=0.45, zorder=3))
    ax[0].plot(lm[0], lm[1], "^", ms=6, mfc=c, mec="k", mew=0.4, ls="none", alpha=0.8, zorder=3)
ax[0].set_aspect("equal"); ax[0].grid(True, lw=0.3, alpha=0.4)
ax[0].set_xlim(-2.2,2.2); ax[0].set_ylim(-1.7,2.1)
ax[0].set_xlabel("x (m)"); ax[0].set_ylabel("y (m)")
ax[0].set_title(f"(a) 6 of {NG} random geometries")
ax[0].legend(loc="lower left")

# (b) pooled collision rate
tot = NG*M
rates = [100*np.sum(COLL[:,i])/tot for i in range(len(modes))]
cis = [wil(int(np.sum(COLL[:,i])), tot) for i in range(len(modes))]
xlab = ["oracle","odom\n+fixed","slam\n+fixed","slam\n+cov"]
bcol = [GREY, VERM, PURPLE, BLUE]
xb = np.arange(len(modes))
ax[1].bar(xb, rates, 0.62, color=bcol)
for i,(r,(lo,hi)) in enumerate(zip(rates,cis)):
    ax[1].plot([i,i],[lo,hi], color="k", lw=1.0)
    ax[1].text(i, hi+0.8, f"{r:.1f}%", ha="center", fontsize=7.5)
ax[1].set_xticks(xb); ax[1].set_xticklabels(xlab)
ax[1].set_ylabel("pooled collision rate (%)")
ax[1].set_ylim(0, max(rates)*1.2)
ax[1].set_title(f"(b) pooled over {tot} trials ($\\gamma$ frozen)")
ax[1].grid(True, axis="y", lw=0.3, alpha=0.4)

# (c) per-geometry: fixed collides everywhere, cov-aware nowhere
sf = COLL[:, modes.index("slam_fixed")]/M*100
sc = COLL[:, modes.index("slam_cov")]/M*100
order = np.argsort(-sf)
xg = np.arange(NG)
ax[2].bar(xg, sf[order], 0.8, color=PURPLE, label="slam + fixed")
ax[2].plot(xg, sc[order], "o", color=BLUE, ms=4, label="slam + cov-aware (all 0)")
ax[2].set_xlabel("geometry (sorted by fixed-margin collisions)")
ax[2].set_ylabel("collision rate (%)")
ax[2].set_title("(c) per-geometry")
ax[2].legend(loc="upper right"); ax[2].grid(True, axis="y", lw=0.3, alpha=0.4)

fig.suptitle("Randomized-geometry transfer: $\\gamma=2$ (tuned once) stays collision-free across unseen layouts", fontsize=9, y=1.03)
fig.tight_layout()
fig.savefig(OUT/"um_randgeom.pdf", bbox_inches="tight")
fig.savefig(OUT/"um_randgeom.png", dpi=300, bbox_inches="tight"); plt.close(fig)
print("saved um_randgeom")
