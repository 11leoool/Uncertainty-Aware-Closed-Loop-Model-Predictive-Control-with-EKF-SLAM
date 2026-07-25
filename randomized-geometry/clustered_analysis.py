# clustered_analysis.py -- canonical analysis for the randomized-geometry study.
#
# Trials are CLUSTERED within geometries, so an interval computed as though they
# were independent understates uncertainty. This produces the intervals reported
# in the manuscript:
#
#   * naive trial-level Wilson interval  (reference only)
#   * cluster bootstrap over GEOMETRIES  (resampled with replacement)
#   * geometry-level Wilson interval     (used when a count is zero, where the
#                                         bootstrap is degenerate)
#   * design effect, quantifying departure from independence
#   * geometry-level sign test, which assumes no independence across trials
#
# Prefers the v2 dataset (disjoint screening/evaluation seeds + continuous
# collision scoring). For v2 the CONTINUOUS scoring is primary, because checking
# clearance only at sampled states can miss an inter-sample dip.
#
# Reproduce:  python clustered_analysis.py

import numpy as np
import scipy.io as sio
from scipy.stats import binomtest
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE / "um_randgeom_v2_results.mat"
V1 = HERE / "um_randgeom_results.mat"
MAT = V2 if V2.exists() else V1
IS_V2 = MAT == V2

B, SEED = 20000, 7
d = sio.loadmat(MAT, squeeze_me=True, struct_as_record=False)
modes = [str(m) for m in d["modes"]]
NG = int(d["NG"])
M = int(d["M_EVAL"]) if IS_V2 else int(d["M"])
N = NG * M
rng = np.random.default_rng(SEED)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    den = 1 + z * z / n
    ctr = (p + z * z / (2 * n)) / den
    hw = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, ctr - hw), min(1.0, ctr + hw))


def report(COLL, label):
    print(f"\n### {label}")
    print(f"{'mode':<14}{'rate':>10}{'naive (trial)':>22}{'cluster bootstrap':>24}"
          f"{'geometry-level':>22}{'DEFF':>7}")
    for i, m in enumerate(modes):
        k = COLL[:, i].sum()
        rate = k / N
        lo_n, hi_n = wilson(k, N)
        bs = np.empty(B)
        for b in range(B):
            pick = rng.integers(0, NG, NG)
            bs[b] = COLL[pick, i].sum() / (NG * M)
        lo_c, hi_c = np.quantile(bs, [0.025, 0.975])
        gfail = int((COLL[:, i] > 0).sum())
        lo_g, hi_g = wilson(gfail, NG)
        pg = COLL[:, i] / M
        var_naive = rate * (1 - rate) / N
        var_clust = pg.var(ddof=1) / NG if NG > 1 else 0.0
        deff = (var_clust / var_naive) if var_naive > 0 else float("nan")
        print(f"{m:<14}{int(k):>4}/{N:<5}{100*rate:>5.2f}%  "
              f"[{100*lo_n:>5.2f},{100*hi_n:>6.2f}]     "
              f"[{100*lo_c:>5.2f},{100*hi_c:>6.2f}]    "
              f"{gfail:>2}/{NG} [{100*lo_g:>5.2f},{100*hi_g:>6.2f}] {deff:>6.1f}")
    i_cov, i_fix = modes.index("slam_cov"), modes.index("slam_fixed")
    b10 = int((COLL[:, i_fix] > COLL[:, i_cov]).sum())
    b01 = int((COLL[:, i_cov] > COLL[:, i_fix]).sum())
    p_sign = binomtest(b10, b10 + b01, 0.5).pvalue if (b10 + b01) > 0 else 1.0
    print(f"  geometry-level sign test, slam_cov vs slam_fixed: cov better in "
          f"{b10}/{NG}, worse in {b01}, tied {NG-b10-b01}  ->  p = {p_sign:.3g}")


print(f"Dataset: {MAT.name}  ({'v2: disjoint seeds + continuous scoring' if IS_V2 else 'v1'})")
print(f"{NG} geometries x {M} evaluation seeds = {N} trials/mode")
print(f"Cluster bootstrap: {B} resamples of GEOMETRIES (seed {SEED})")

if IS_V2:
    report(np.atleast_2d(d["COLLC"]).astype(float), "PRIMARY: continuous (segment) collision scoring")
    report(np.atleast_2d(d["COLL"]).astype(float), "reference: sampled-state scoring")
else:
    report(np.atleast_2d(d["COLL"]).astype(float), "sampled-state scoring")

print("\nNotes:")
print("* A zero count makes the cluster bootstrap degenerate ([0,0]); quote the")
print("  geometry-level Wilson bound there instead.")
print("* DEFF >> 1 confirms trials are not independent within a geometry.")
