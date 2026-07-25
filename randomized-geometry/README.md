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
solvable on the screening seeds; the oracle still collides on other seeds of accepted layouts, so acceptance is evidence of solvability rather than a guarantee. A collision elsewhere is attributable to estimation or margin rather than to an
impossible scenario. **γ is frozen at 2** (the value tuned on the original geometry);
no per-geometry retuning. 30 feasible geometries × 20 shared seeds = 600 trials/mode.

## Result — γ transfers, as a large reduction rather than elimination

**These are the v2 numbers.** An earlier version of this study (v1) reported 0/600 for
the covariance-aware margin. That result was invalidated by two protocol defects, both
fixed in `run_um_randgeom_v2.m`: screening realizations were reused as evaluation
realizations, and collisions were scored only at sampled states. See the correction
commit for detail. **Do not cite the v1 numbers.**

Primary scoring is continuous (minimum distance from each travelled segment to the
obstacle); intervals are cluster bootstrap over geometries, from `clustered_analysis.py`.

| Mode | Collisions | Rate (cluster-bootstrap 95% CI) | Worst clearance | Geometries with ≥1 |
|---|---|---|---|---|
| oracle (clairvoyant, fixed δ₀) | 2 / 600 | 0.33% [0.00, 1.00] | −0.003 | 1 / 30 |
| odom + fixed δ₀ | 179 / 600 | 29.8% [22.3, 37.3] | −0.280 | 26 / 30 |
| slam + fixed δ₀ | 161 / 600 | 26.8% [20.3, 33.7] | −0.240 | 29 / 30 |
| **slam + covariance-aware** | **3 / 600** | **0.50% [0.00, 1.33]** | **−0.069** | **2 / 30** |

- The covariance-aware margin **does not eliminate collisions**: 3 of 600 trials, in 2 of
  the 30 geometries, worst penetration 6.9 cm.
- It reduces the collision rate by roughly fifty-fold relative to a fixed margin
  (0.50% vs 26.8–29.8%), and is better in **29 of 30 geometries**, never worse, tied in
  one (geometry-level sign test, p = 3.7 × 10⁻⁹).
- The design effect is 3.7–4.4 for the fixed-margin arms, so trials are clearly not
  independent within a geometry; pooled Wilson intervals over 600 trials would be
  anti-conservative.

**Takeaway:** the operating point tuned on one geometry transfers to unseen layouts in
the sense of a large, consistent reduction — not as a guarantee of zero collisions. This
addresses the single-geometry / in-sample-tuning objection while being explicit about
what the evidence does and does not support.

**Known limitation of the scoring.** Continuous scoring interpolates a straight chord
between consecutive sampled poses; the robot actually follows a constant-twist arc. The
chord is a good approximation at these step sizes and does not affect the covariance-aware
penetrations (5.9–6.9 cm, far larger than any chord–arc discrepancy), but it could matter
for millimetre-scale grazes such as the oracle's.

## Files
**Current (v2, use these):** `run_um_randgeom_v2.m` (disjoint screening/evaluation
streams + continuous collision scoring), results `um_randgeom_v2_results.{mat,txt}`,
analysis `clustered_analysis.py`.

**Superseded (v1, retained for provenance only — numbers invalidated):**
`run_um_randgeom.m` (sampler + oracle feasibility filter + 4-mode run), deps
`mc_run_trial_um_obs.m`, `mc_build_mpc_dyn.m`, `mc_ekf_step.m`, `static_obs_ekf.m`;
results `um_randgeom_results.*`; `make_um_randgeom_fig.py` (v1 figure script;
the committed `um_randgeom.pdf` is regenerated from v2 data)
(sample geometries, pooled bars, per-geometry breakdown).

## Reproduce
```matlab
run_um_randgeom_v2        % ~30 min (sampling + screening + 2400 evaluation trials)
```
```
python clustered_analysis.py
```

Unlike the mismatch studies, **this one is paper-relevant** — it closes a named reviewer
concern and could be added to the manuscript as a transfer/generalization result.
