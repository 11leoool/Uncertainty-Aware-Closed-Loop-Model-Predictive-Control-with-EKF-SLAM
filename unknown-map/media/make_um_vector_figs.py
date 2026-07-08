# make_um_vector_figs.py -- submission-grade vector figures for the unknown-map
# paper (10_unknown_map), regenerated from the MATLAB .mat results.
# Okabe-Ito palette, editable TrueType text (pdf.fonttype 42), tight bboxes.
import numpy as np
import scipy.io as sio
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent   # code\Unknown map
OUT = Path(__file__).resolve().parent

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

# Okabe-Ito
BLACK  = "#000000"
ORANGE = "#E69F00"
SKY    = "#56B4E9"
GREEN  = "#009E73"
BLUE   = "#0072B2"
VERM   = "#D55E00"
PURPLE = "#CC79A7"
GREY   = "#999999"

def load(name):
    return sio.loadmat(ROOT / name, squeeze_me=True, struct_as_record=False)

def savefig(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {stem}.pdf/.png")

def draw_scene(ax, cfg, obs=None, r_sense=None):
    lm = np.atleast_2d(cfg.lm)
    wps = np.atleast_2d(cfg.wps)
    if obs is not None:
        ax.add_patch(Circle((obs[0], obs[1]), obs[2], facecolor="#e8b0b0",
                            edgecolor="none", zorder=1))
    ax.plot(lm[0], lm[1], "^", ms=7, mfc="#F0E442", mec="k", mew=0.6,
            ls="none", zorder=3, label="landmarks")
    ax.plot(wps[0], wps[1], "*", ms=9, mfc=GREEN, mec="k", mew=0.4,
            ls="none", zorder=3, label="waypoints")
    ax.set_aspect("equal")
    ax.grid(True, lw=0.3, alpha=0.5)
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")

# ---------------------------------------------------------------- patrol pair
d = load("mc_um4_results.mat")
res, cfg = d["results"], d["cfg"]

fig, ax = plt.subplots(figsize=(3.35, 2.9))
draw_scene(ax, cfg)
for mode, col, ls, lab in [("oracle", BLACK, "-", "oracle"),
                           ("odom", VERM, "-.", "odometry"),
                           ("slam", BLUE, "-", "SLAM")]:
    tr = getattr(res, mode).sample.traj_true
    ax.plot(tr[0], tr[1], color=col, ls=ls, lw=1.2, label=lab, zorder=2)
ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.2)
savefig(fig, "um_patrol_traj")

fig, ax = plt.subplots(figsize=(3.35, 2.9))
T = 0.05
so = res.odom.sample.sig_t
ss = res.slam.sample.sig_t
ax.plot(np.arange(1, so.size + 1) * T, so, color=VERM, ls="-.", lw=1.2,
        label="odometry (dead reckoning)")
ax.plot(np.arange(1, ss.size + 1) * T, ss, color=BLUE, lw=1.2, label="SLAM")
ax.set_xlabel("time (s)")
ax.set_ylabel(r"$\sigma_{\max}$ (m)")
ax.grid(True, lw=0.3, alpha=0.5)
ax.legend(loc="upper left")
savefig(fig, "um_patrol_sigma")

# --------------------------------------------------------- static obstacle pair
d = load("um_obstacle_results.mat")
res, cfg = d["results"], d["cfg"]
obs = np.ravel(cfg.obs)
si = int(np.flatnonzero(np.ravel(res.slam_fixed.collided))[0]) \
     if np.any(res.slam_fixed.collided) else 0

fig, ax = plt.subplots(figsize=(3.35, 2.9))
draw_scene(ax, cfg, obs=obs)
for mode, col, ls, lab in [("oracle", BLACK, "-", "oracle"),
                           ("odom_fixed", VERM, "-.", "odom + fixed"),
                           ("slam_fixed", PURPLE, "--", "SLAM + fixed"),
                           ("slam_cov", BLUE, "-", "SLAM + cov-aware")]:
    tr = getattr(res, mode).sample[si].traj_true
    ax.plot(tr[0], tr[1], color=col, ls=ls, lw=1.1, label=lab, zorder=2)
ax.legend(loc="upper left", handlelength=1.7, borderaxespad=0.2)
savefig(fig, "um_obs_traj")

rc = res.slam_cov.sample[si]
fig, ax = plt.subplots(figsize=(3.35, 2.9))
T = 0.2
t = np.arange(1, rc.sig_t.size + 1) * T
ax.plot(t, rc.sig_t, color=BLUE, lw=1.2, label=r"$\sigma_{\max}$ (robot)")
ax.plot(t, rc.delta_t, color=BLUE, ls="--", lw=1.2, label="applied inflation")
oe = np.array(rc.obserr_t, dtype=float)
ax.plot(t, oe, color=BLACK, ls=":", lw=1.2, label="obstacle est. error")
ax.set_xlabel("time (s)"); ax.set_ylabel("(m)")
ax.grid(True, lw=0.3, alpha=0.5)
ax.legend(loc="upper right")
savefig(fig, "um_obs_margin")

# ------------------------------------------------------------- static gamma sweep
d = load("um_gamma_ablation_results.mat")
sw = d["sweep"]
fig, ax = plt.subplots(figsize=(4.6, 2.9))
ax.plot(sw.gamma, sw.coll, "-o", color=VERM, lw=1.4, ms=4, mfc=VERM)
ax.set_xlabel(r"$\gamma$")
ax.set_ylabel("Collision rate (%)", color=VERM)
ax.tick_params(axis="y", labelcolor=VERM)
ax.set_ylim(-0.7, None)
ax2 = ax.twinx()
ax2.spines["top"].set_visible(False)
ax2.plot(sw.gamma, sw.path, "-s", color=BLUE, lw=1.4, ms=4, mfc=BLUE)
ax2.set_ylabel("Mean path length (m)", color=BLUE)
ax2.tick_params(axis="y", labelcolor=BLUE)
ax.grid(True, lw=0.3, alpha=0.5)
savefig(fig, "um_gamma_sweep")

# ------------------------------------------------------------ dynamic gamma sweep
d = load("um_gamma_dyn_results.mat")
sw = d["sweep"]
fig, ax = plt.subplots(figsize=(4.6, 2.9))
ax.plot(sw.gamma, sw.coll, "-o", color=VERM, lw=1.4, ms=4, mfc=VERM)
ax.set_xlabel(r"$\gamma$")
ax.set_ylabel("Collision rate (%)", color=VERM)
ax.tick_params(axis="y", labelcolor=VERM)
ax.set_ylim(-2, None)
ax2 = ax.twinx()
ax2.spines["top"].set_visible(False)
ax2.plot(sw.gamma, sw.path, "-s", color=BLUE, lw=1.4, ms=4, mfc=BLUE)
ax2.set_ylabel("Mean path length (m)", color=BLUE)
ax2.tick_params(axis="y", labelcolor=BLUE)
ax.grid(True, lw=0.3, alpha=0.5)
savefig(fig, "um_gamma_dyn")

# ------------------------------------------- dynamic disk-snapshot side-by-side
d = load("um_dynamic_results.mat")
res, cfg = d["results"], d["cfg"]
r_rob = 0.15
r_obs = float(cfg.obs_r)
si = int(np.flatnonzero(np.ravel(res.static.collided))[0]) \
     if np.any(res.static.collided) else 0

def snapshot_panel(ax, sample, cfg, title, nsnap=7):
    # Simultaneous robot/obstacle snapshots share a progress colour; the
    # obstacle disk is the translucent twin and turns red on overlap.
    tr, ot = sample.traj, sample.otraj
    K = tr.shape[1]
    idx = np.unique(np.round(np.linspace(0, K - 1, nsnap)).astype(int))
    # always include the closest-approach instant so a collision is visible
    sep = np.hypot(tr[0, :] - ot[0, :], tr[1, :] - ot[1, :])
    idx = np.unique(np.append(idx, int(np.argmin(sep))))
    cmap = plt.get_cmap("viridis")
    lm = np.atleast_2d(cfg.lm); wps = np.atleast_2d(cfg.wps)
    ax.plot(lm[0], lm[1], "^", ms=6, mfc="#F0E442", mec="k", mew=0.5, ls="none")
    ax.plot(wps[0], wps[1], "*", ms=8, mfc=GREEN, mec="k", mew=0.4, ls="none")
    ax.plot(tr[0], tr[1], color=GREY, lw=0.7, alpha=0.8, zorder=1)
    ax.plot(ot[0], ot[1], color="#c9a0a0", lw=0.7, ls="--", alpha=0.9, zorder=1)
    for k in idx:
        frac = k / max(K - 1, 1)          # colour by true time fraction
        col = cmap(frac)
        overlap = np.hypot(tr[0, k] - ot[0, k], tr[1, k] - ot[1, k]) < (r_rob + r_obs)
        oc = "#d62728" if overlap else col
        ax.add_patch(Circle((ot[0, k], ot[1, k]), r_obs, facecolor=oc,
                            edgecolor=col, lw=0.8, alpha=0.35 if not overlap else 0.85,
                            zorder=2))
        ax.add_patch(Circle((tr[0, k], tr[1, k]), r_rob, facecolor=col,
                            edgecolor="k", lw=0.3, alpha=0.9, zorder=3))
    ax.set_xlim(-2.15, 2.15); ax.set_ylim(-1.0, 1.9)
    ax.set_aspect("equal")
    ax.grid(True, lw=0.3, alpha=0.5)
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
    ax.set_title(title, fontsize=8)

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1), sharey=True)
snapshot_panel(axes[0], res.static.sample[si], cfg,
               "(a) frozen estimate, fixed margin")
snapshot_panel(axes[1], res.cv_cov.sample[si], cfg,
               "(b) covariance-aware margin")
axes[1].set_ylabel("")
fig.tight_layout()
savefig(fig, "um_dyn_snapshots")

print("all figures done")
