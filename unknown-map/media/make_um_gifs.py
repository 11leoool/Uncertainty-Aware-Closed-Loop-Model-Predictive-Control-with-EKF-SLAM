# make_um_gifs.py -- supplementary GIF animations for the unknown-map dynamic
# obstacle study (same trial for both strategies): robot and obstacle disks,
# sensing-range circle, trails; obstacle flashes red on physical overlap.
import numpy as np
import scipy.io as sio
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
})

GREEN = "#009E73"; GREY = "#999999"; BLUE = "#0072B2"; RED = "#d62728"
PINK = "#e8b0b0"

d = sio.loadmat(ROOT / "um_dynamic_results.mat", squeeze_me=True,
                struct_as_record=False)
res, cfg = d["results"], d["cfg"]
r_rob, r_obs = 0.15, float(cfg.obs_r)
r_sense = float(cfg.r_sense)
T = 0.2
si = int(np.flatnonzero(np.ravel(res.static.collided))[0]) \
     if np.any(res.static.collided) else 0

def make_gif(sample, title, fname, step=2, fps=8):
    tr, ot = sample.traj, sample.otraj
    K = tr.shape[1]
    frames = list(range(0, K, step))

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    lm = np.atleast_2d(cfg.lm); wps = np.atleast_2d(cfg.wps)
    ax.plot(lm[0], lm[1], "^", ms=8, mfc="#F0E442", mec="k", mew=0.6, ls="none")
    ax.plot(wps[0], wps[1], "*", ms=11, mfc=GREEN, mec="k", mew=0.5, ls="none")
    te = getattr(sample, "traj_est", None)      # believed pose (ghost), if stored
    trail_r, = ax.plot([], [], color=BLUE, lw=1.0, alpha=0.7)
    trail_o, = ax.plot([], [], color="#c9a0a0", lw=1.0, ls="--", alpha=0.8)
    disk_r = Circle((tr[0, 0], tr[1, 0]), r_rob, facecolor=BLUE,
                    edgecolor="k", lw=0.5, zorder=4)
    disk_o = Circle((ot[0, 0], ot[1, 0]), r_obs, facecolor=PINK,
                    edgecolor="k", lw=0.5, zorder=3)
    ring = Circle((tr[0, 0], tr[1, 0]), r_sense, facecolor="none",
                  edgecolor=GREY, ls=":", lw=1.0, zorder=2)
    ghost = Circle((tr[0, 0], tr[1, 0]), r_rob, facecolor="none",
                   edgecolor=BLUE, ls="--", lw=1.0, alpha=0.9, zorder=5)
    for p in (disk_o, disk_r, ring, ghost):
        ax.add_patch(p)
    ghost.set_visible(te is not None)
    txt = ax.text(0.02, 0.98, "", transform=ax.transAxes, va="top", fontsize=8)
    ax.set_xlim(-2.15, 2.15); ax.set_ylim(-1.1, 2.0)
    ax.set_aspect("equal"); ax.grid(True, lw=0.3, alpha=0.5)
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
    ax.set_title(title, fontsize=10)

    collided_flag = {"hit": False}

    def update(k):
        trail_r.set_data(tr[0, :k + 1], tr[1, :k + 1])
        trail_o.set_data(ot[0, :k + 1], ot[1, :k + 1])
        disk_r.center = (tr[0, k], tr[1, k])
        disk_o.center = (ot[0, k], ot[1, k])
        ring.center = (tr[0, k], tr[1, k])
        if te is not None:
            ghost.center = (te[0, k], te[1, k])   # believed pose (dashed)
        d_ro = np.hypot(tr[0, k] - ot[0, k], tr[1, k] - ot[1, k])
        in_range = d_ro <= r_sense
        overlap = d_ro < (r_rob + r_obs)
        if overlap:
            collided_flag["hit"] = True
        disk_o.set_facecolor(RED if collided_flag["hit"] else
                             (PINK if in_range else "#f2dcdc"))
        disk_o.set_alpha(1.0 if in_range or collided_flag["hit"] else 0.55)
        status = "COLLISION" if collided_flag["hit"] else \
                 ("obstacle in range" if in_range else "obstacle out of range")
        txt.set_text(f"t = {k * T:5.1f} s   {status}")
        return trail_r, trail_o, disk_r, disk_o, ring, txt

    anim = FuncAnimation(fig, update, frames=frames, blit=False)
    anim.save(OUT / fname, writer=PillowWriter(fps=fps), dpi=110)
    plt.close(fig)
    print(f"saved {fname} ({len(frames)} frames)")

make_gif(res.static.sample[si],
         "Frozen estimate, fixed margin (collision)",
         "um_dyn_static_collision.gif")
make_gif(res.cv_cov.sample[si],
         "Covariance-aware margin (safe; dashed = believed pose)",
         "um_dyn_cvcov_safe.gif")
print("done")
