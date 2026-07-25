# clustered_analysis.py -- canonical analysis for the randomized-geometry study.
#
# The 600 trials are CLUSTERED within 30 geometries, so an interval computed as
# though they were independent understates uncertainty. This script produces the
# intervals reported in the manuscript's randomized-geometry table and figure:
#
#   * naive trial-level Wilson interval  (shown only for reference)
#   * cluster bootstrap over GEOMETRIES  (resampled with replacement)
#   * geometry-level Wilson interval     (used when a count is zero, where the
#                                         bootstrap is degenerate)
#   * design effect, quantifying how far from independent the trials are
#   * geometry-level sign test, which assumes no independence across trials
#
# Reproduce:  python clustered_analysis.py
# Input:      um_randgeom_results.mat  (committed alongside this script)

import numpy as np
import scipy.io as sio
from scipy.stats import binomtest
from pathlib import Path

MAT = Path(__file__).resolve().parent / "um_randgeom_results.mat"
B = 20000          # bootstrap resamples
SEED = 7           # fixed for reproducibility

d = sio.loadmat(MAT, squeeze_me=True, struct_as_record=False)
COLL = np.atleast_2d(d["COLL"]).astype(float)     # [NG x nmodes] collisions per geometry
modes = [str(m) for m in d["modes"]]
NG, M = int(d["NG"]), int(d["M"])
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


print(f"Randomized-geometry study: {NG} geometries x {M} seeds = {N} trials/mode")
print(f"Cluster bootstrap: {B} resamples of GEOMETRIES (seed {SEED})\n")
print(f"{'mode':<14}{'rate':>8}{'naive (trial)':>22}{'cluster bootstrap':>24}"
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

    print(f"{m:<14}{100*rate:>7.2f}%  [{100*lo_n:>5.2f},{100*hi_n:>6.2f}]     "
          f"[{100*lo_c:>5.2f},{100*hi_c:>6.2f}]    "
          f"{gfail:>2}/{NG} [{100*lo_g:>5.2f},{100*hi_g:>6.2f}] {deff:>6.1f}")

print("\nNotes:")
print("* A zero count makes the cluster bootstrap degenerate ([0,0]); the")
print("  geometry-level Wilson bound is the honest interval to quote there.")
print("* DEFF >> 1 confirms trials are not independent within a geometry, so the")
print("  naive trial-level interval is anti-conservative for the fixed-margin arms.")

i_cov, i_fix = modes.index("slam_cov"), modes.index("slam_fixed")
b10 = int((COLL[:, i_fix] > COLL[:, i_cov]).sum())
b01 = int((COLL[:, i_cov] > COLL[:, i_fix]).sum())
p_sign = binomtest(b10, b10 + b01, 0.5).pvalue if (b10 + b01) > 0 else 1.0
print(f"\nGeometry-level sign test, slam_cov vs slam_fixed:")
print(f"  cov better in {b10}/{NG} geometries, worse in {b01}, "
      f"tied in {NG-b10-b01}  ->  p = {p_sign:.3g}")
