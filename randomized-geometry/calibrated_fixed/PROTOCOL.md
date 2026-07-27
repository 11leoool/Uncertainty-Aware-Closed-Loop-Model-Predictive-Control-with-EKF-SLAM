# Frozen protocol: calibrated fixed-margin baseline across randomized geometries

**Status: FROZEN before any run. No calibration or evaluation outcome exists.**

## Question

The randomized-geometry headline compares the covariance-aware margin against
the small fixed buffer delta0 = 0.06 m, which the static ablation already showed
is under-margined. This experiment asks: does a **single fixed constant,
calibrated on disjoint training layouts and frozen**, transfer across the same
30 evaluation geometries as well as the covariance-aware margin does?

## Interpretation, declared in advance

- If the calibrated constant performs comparably: the transfer contribution is
  framed as **automatic online margin generation and transfer** (the framing the
  corrected manuscript already uses), not superior adaptivity.
- If the covariance-aware margin is better (paired, two-way inference): that is
  reportable as evidence for adaptivity under transfer.
- Collision counts this small may be inconclusive; an inconclusive outcome is
  reported as inconclusive.

## Design

1. **Calibration (training) phase — all streams disjoint from every prior study:**
   - geometry sampling stream `rng(3030)`, N_TRAIN = 10 layouts, same
     `sample_geom` sampler and oracle-feasibility screen as the published run;
   - screening noises `rng(9500)` (3 seeds, never evaluated);
   - calibration noises `rng(7100)` stream, M_CAL = 20 per layout;
   - run `slam_cov` (gamma = 2 frozen, as published) on all 200 trials;
   - **c = mean over trials of (per-trial mean of applied inflation delta_t,
     conditional on delta_t > 0)** — the active-step conditional mean, equal
     trial weighting, per the project estimand rule.
2. **Freeze c** — recorded below before any evaluation trial runs.
3. **Evaluation:** replay the published v2 stream setup verbatim (screening 9001,
   evaluation 4242, geometry 2024) to regenerate the identical 30 geometries and
   20 evaluation noises; run `slam_fixed` with `fixed_extra = c` over the full
   30 x 20 grid.
   - **Replay validation (mandatory):** re-run plain `slam_fixed` (delta0) on the
     9 cells (g,k) in {1..3}x{1..3} and assert `min_clear_cont` matches the
     stored `um_randgeom_v2_results.mat` TRIAL records bit-exactly. If this
     fails, the replay is not the published experiment and the evaluation stops.
4. **Analysis:** paired per (geometry, seed) against the stored `slam_cov` and
   `slam_fixed` records; continuous collision scoring primary; two-way
   geometry-and-seed bootstrap (B = 20000, `rng(11)`) for rate differences;
   path length over all trials; worst clearance. No p-value is promoted to a
   claim without the two-way interval.

- Frozen `c`: **0.2731 m** (per-trial conditional means: median 0.2762,
  IQR 0.1729, min 0.1400, max 0.6588, n = 200; mean active fraction 0.875).
  Recorded 2026-07-28, before any evaluation trial ran. Note for context: this
  is ~4.5x the delta0 = 0.06 m baseline and ~1.9x the fixed-scene matched
  constant (0.146 m) — the randomized layouts demand larger margins than the
  single published scene, which is itself informative.

## Outcome (recorded 2026-07-28, after the run)

Replay validation PASSED (9/9 cells bit-exact). Calibrated constant 1/600
collisions vs the adaptive margin's 3/600; paired two-way rate difference
+0.33 pp [-0.67, +1.83], spanning zero. Per the pre-declared rules: the two are
COMPARABLE, the transfer advantage over the delta0 baselines is carried by
margin size, and the covariance-aware controller's demonstrated role is
generating that margin automatically (the constant was obtained by running the
adaptive controller on training layouts and summarising it).
