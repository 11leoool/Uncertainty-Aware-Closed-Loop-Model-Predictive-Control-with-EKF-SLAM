# Robustness to filter mis-calibration (filter-consistency sweep)

Self-contained archive of the filter-consistency Monte Carlo used in the
manuscript's *Robustness to Filter Mis-calibration* subsection. It asks how the
covariance-aware safety margin behaves when the filter's noise model is wrong.

## Design
On the sensed static-obstacle scenario, the filter's **assumed** process and
measurement noise is scaled by a factor while the **true** noise is drawn from
the original, fixed distribution on the same seeds. The ratio of assumed to true
variance sweeps the filter from overconfident (ratio < 1) through consistent
(ratio = 1) to conservative (ratio > 1). γ is frozen at 2, so only the filter's
self-assessment changes — this is genuine (in)consistency, not a post-hoc
rescaling of the reported covariance. `N_MC = 50` per rung.

## Result — graceful degradation, wide safe band
| Assumed/True | Regime | Collision % (Wilson) | Worst clear. | Reported σ | True loc. err |
|---|---|---|---|---|---|
| 0.10 | 10× overconfident | 2.0% [0.4, 10.5] | −0.016 | 0.022 | 0.077 |
| 0.25 | 4× overconfident | 0.0% [0.0, 7.1] | +0.008 | 0.034 | 0.076 |
| 0.50 | 2× overconfident | 0.0% [0.0, 7.1] | +0.037 | 0.048 | 0.076 |
| 1.00 | consistent | 0.0% [0.0, 7.1] | +0.073 | 0.067 | 0.076 |
| 2.00 | 2× conservative | 0.0% [0.0, 7.1] | +0.129 | 0.095 | 0.075 |
| 4.00 | 4× conservative | 0.0% [0.0, 7.1] | +0.206 | 0.133 | 0.075 |

Covariance-blind fixed-margin reference: **14.0% [7.0, 26.2]** on the same obstacle.

- The covariance-aware margin stays collision-free (0/50) down to a **fourfold
  overconfident** filter, and admits only a single shallow graze (2%, −1.6 cm) at
  **tenfold** overconfidence — versus 14% for the covariance-blind fixed margin.
- Worst clearance degrades **smoothly**; when the adaptive term collapses, the
  constant base buffer δ₀ = 0.06 m sets a floor (graceful, not abrupt, failure).
- The **true** localization error is essentially flat (~0.076 m): the EKF point
  estimate is robust to noise mis-tuning; only the covariance (which the margin
  consumes) is sensitive.

**Takeaway:** the margin inherits the honesty of its filter, but the fixed floor
and the point estimate's insensitivity leave a broad band of safe operation before
overconfidence bites.

## Files
`run_um_consistency.m` (sweep driver; `cfg.filt_scale` scales the filter's assumed
P/Q/R), deps `mc_run_trial_um_obs.m`, `mc_build_mpc_dyn.m`, `mc_ekf_step.m`,
`static_obs_ekf.m`; results `um_consistency_results.{txt,mat}`;
`make_um_consistency_fig.py` → `um_consistency.pdf`.

## Reproduce
```matlab
run_um_consistency        % ~30-40 min (6 rungs x 4 modes x 50 seeds)
```
```
python make_um_consistency_fig.py
```
