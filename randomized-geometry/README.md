# Randomized-geometry transfer — static obstacle

Self-contained archive of the randomized-geometry Monte Carlo. Answers the reviewers'
top concern directly: γ = 2 was tuned on **one** fixed geometry; does it transfer to
layouts it never saw? This is paper-grade (not just thesis) — it closes the
"in-sample / single-geometry" critique with data.

## Design
Per geometry, the obstacle position and all three landmark positions are randomized
under feasibility constraints (obstacle within contact distance of a random path leg
so avoidance is genuinely required; landmarks ≥ 0.45 m from the obstacle, spread apart,
≥ 1 within sensing range of the start). An **oracle feasibility filter** then keeps only
geometries the ideal-state controller solves collision-free — so every tested layout is
provably solvable and any collision elsewhere is an estimation/margin failure, not an
impossible scenario. **γ is frozen at 2** (the value tuned on the original geometry);
no per-geometry retuning. 30 feasible geometries × 20 shared seeds = 600 trials/mode.

## Result — γ transfers perfectly

| Mode | Pooled collision % (Wilson) | Worst clearance | Geometries with ≥1 collision |
|---|---|---|---|
| oracle (clairvoyant, fixed δ₀) | 1.33% [0.68, 2.61] | −0.013 | 5 / 30 |
| odom + fixed δ₀ | 36.0% [32.3, 39.9] | −0.283 | **30 / 30** |
| slam + fixed δ₀ | 31.0% [27.4, 34.8] | −0.180 | **30 / 30** |
| **slam + covariance-aware** | **0.00% [0.00, 0.64]** | **+0.025** | **0 / 30** |

- The covariance-aware margin, with γ frozen, is **collision-free across all 30 unseen
  random geometries** (0 of 600 trials; worst true clearance still positive at +2.5 cm).
- Fixed margins collide in **every single geometry** (30/30), pooling to 31–36%.
- Note the covariance-aware margin even edges out the clairvoyant oracle (0% vs 1.3%):
  the oracle uses the small fixed δ₀ = 0.06 m, whereas the covariance-aware margin is
  larger, so it clears the tightest random layouts where the oracle's small buffer
  occasionally grazes. (Not magic — a larger, uncertainty-scaled margin vs a small
  constant one.)

**Takeaway:** the operating point tuned on one geometry generalizes; the covariance-
aware margin's collision-free behaviour is not an artifact of the geometry it was tuned
on. This directly closes the single-geometry / in-sample-tuning objection.

## Files
`run_um_randgeom.m` (sampler + oracle feasibility filter + 4-mode run), deps
`mc_run_trial_um_obs.m`, `mc_build_mpc_dyn.m`, `mc_ekf_step.m`, `static_obs_ekf.m`;
results `um_randgeom_results.*`; `make_um_randgeom_fig.py` → `um_randgeom.pdf`
(sample geometries, pooled bars, per-geometry breakdown).

## Reproduce
```matlab
run_um_randgeom           % ~1.5-2 h (sampling + 2400 trials)
```
```
python make_um_randgeom_fig.py
```

Unlike the mismatch studies, **this one is paper-relevant** — it closes a named reviewer
concern and could be added to the manuscript as a transfer/generalization result.
