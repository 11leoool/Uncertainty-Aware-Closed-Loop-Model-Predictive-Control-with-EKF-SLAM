# Uncertainty-Aware Closed-Loop MPC with EKF-SLAM

MATLAB/CasADi code for the paper *"Uncertainty-Aware Closed-Loop Model Predictive
Control with EKF-SLAM for Safe Navigation of Nonholonomic Mobile Robots."*

> ### ⚠ Interim correction notice (2026-07-26)
>
> The **dynamic matched-margin ablation** (`unknown-map/um_dyn_ablation_results.txt`;
> the dynamic row of the manuscript's ablation table) is **superseded — do not
> cite**. Its reported "mean applied inflation 4.4419 m" is a logging artifact: a
> stale value was carried over the ~72% of the mission during which no obstacle
> constraint existed. Independently, the fixed comparator was derived from an
> all-mission mean despite being applied only during active tracking, so the two
> arms were not exposure-matched. The conclusion that *margin size rather than
> adaptivity* explains the safety results is **suspended** — not reversed — until
> a corrected rerun under a disjoint calibration/evaluation protocol completes.
> Details: [`unknown-map/SUPERSEDED_dyn_ablation.md`](unknown-map/SUPERSEDED_dyn_ablation.md).
> Frozen rerun protocol: [`unknown-map/DYN_ABLATION_ESTIMAND.md`](unknown-map/DYN_ABLATION_ESTIMAND.md).
>
> The previously linked `paper.pdf` contained the superseded table and has been
> **removed**; a corrected build will be restored here. A narrative correction is
> also pending: the stored trials show the static obstacle passage is *not* the
> landmark-poor, peak-uncertainty stretch the text describes (2 of 3 landmarks
> are in range there in all 50 trials, and pose uncertainty is ≈55% of its peak).
>
> **Unaffected and re-verified against the stored data:** the headline collision
> results (fixed margins 36–88% vs. 0 observed for the covariance-aware margin),
> both γ-sweeps, the randomized-geometry transfer, the mis-calibration sweep,
> and the N_MC=100 robustness runs.

A differential-drive (unicycle) robot is controlled by a constrained nonlinear MPC
(NMPC) that is fed pose estimates from an EKF-SLAM filter. Beyond the usual
certainty-equivalence coupling, the **EKF-SLAM covariance is fed back into the
controller** to size obstacle keep-out constraints online: the keep-out radius grows
when the robot is uncertain and tightens when it is confident.

## Current version: unknown-map SLAM (`unknown-map/`)

The paper's current setting: **the environment is unknown before it is detected.**
Landmarks are initialized from their first in-range observation (finite sensing
range, start-frame convention); obstacles are *discovered*, estimated online by a
2-state EKF (static) or a coasting constant-velocity EKF with track management
(dynamic); the margin is
`delta = delta0 + gamma * sqrt(lambda_max(P_xy + Sigma_obs))` with both covariances
live estimates. Key results (N_MC=50, confirmed at N_MC=100):

- **Tracking + mapping (return-to-start patrol, ~7 m):** SLAM cuts terminal error
  **48%** vs. dead reckoning (0.064 vs 0.122 m) while building the map from scratch
  to **4 cm**.
- **Sensed static obstacle:** fixed margins collide in 36% (odom) / 14% (SLAM);
  the covariance-aware margin recorded **no collision (0 of 50 and 0 of 100 trials) for <1% extra
  path**.
- **Dynamic obstacle, intermittent visibility:** freeze 88%, CV+fixed 44%,
  **cov-aware 0 of 50** (matching the oracle's collision count) at an honest detour cost
  (path 15.4 vs 7.0 m). Without track management the coasting covariance
  blockades the workspace — safety without mission completion.
- **Perception–control coupling (new finding):** the safety detours degrade the
  localization that waypoint precision depends on — waypoints are attained in the
  *belief* frame (mean true miss 0.11 m) while safety remains true-frame because
  the margin scales with exactly the covariance the detour grows.
- **Ablation (static; descriptive):** a fixed constant set to the cov-aware
  controller's *all-mission* mean inflation (0.146 m — about 7% below its
  active-track mean, so the arms were not exactly matched) also recorded zero
  collisions. The **dynamic** matched row is superseded — see the correction
  notice above. What the ablations demonstrate is *self-tuning*: the controller
  finds a working margin online with no calibration sweep.

Animations (dashed circle = believed pose):
[safe run (cv_cov)](unknown-map/media/um_dyn_cvcov_safe.gif) ·
[collision run (frozen estimate)](unknown-map/media/um_dyn_static_collision.gif)

## Earlier surveyed-landmark studies (Monte-Carlo, 50 trials)

- **Tracking:** closing the loop with EKF-SLAM cuts terminal tracking error ~39% vs.
  dead-reckoning (odometry), approaching the ideal-state (oracle) bound.
- **Safety:** with a *fixed* obstacle margin, both SLAM and odometry collide in ~26–28%
  of trials. The **covariance-aware margin recorded no SLAM collision (0 of 50)** at no
  accuracy cost — an outcome neither a fixed margin nor odometry can match.
- **Dynamic obstacle:** a moving obstacle is tracked by a constant-velocity EKF; the
  margin is inflated by the *combined* robot-localization and obstacle-prediction
  covariance. This yields **0 of 50 collisions** where ignoring the motion (48%) or a
  fixed margin (36%) fail.
- **Ablation (size vs. adaptivity):** a *fixed* margin matched to cv_cov's *mean* margin
  records the same zero count (0/50). Per the correction notice above, the
  matched-margin methodology is under corrected evaluation and the
  size-not-adaptivity conclusion is **suspended**; the demonstrated value of the
  covariance-aware shaping is that it **self-tunes** the margin online (no
  hand-tuned constant). The equal-safety path-efficiency comparison is withdrawn
  pending exactly matched arms.
- **Robustness (N_MC=100):** all studies re-run at 100 trials reproduce the findings with
  tighter Wilson 95% CIs — the zero-collision arms are bounded at ≤3.7% upper 95%, the ablation
  conclusion holds, and the γ-sweeps become monotone. See `mc-m100-robustness/mc100_results.txt`.
- **Real-time:** ~4–7 ms per control step (p95 ≤ 15 ms), an order-of-magnitude margin;
  the EKF and the covariance term are negligible.

| Free-space tracking | Stage A (fixed) vs Stage B (cov-aware) | γ sweep (safety–efficiency) |
|---|---|---|
| ![tracking](figures/no_obs_error.png) | ![collision](figures/obs_stageB_collision.png) | ![sweep](figures/gamma_sweep.png) |

### Dynamic (moving) obstacle

Same trial, same obstacle motion — only the avoidance strategy differs. Under a fixed
margin the robot and obstacle disks overlap (collision); the covariance-aware margin
detours and stays clear:

![dynamic side-by-side](figures/dyn_sidebyside.png)

| Dynamic γ sweep | Animations |
|---|---|
| ![dyn sweep](figures/dyn_gamma_sweep.png) | [safe run (cv_cov)](figures/dyn_anim_cvcov_safe.gif) · [collision run (static)](figures/dyn_anim_static_collision.gif) |

## Repository structure

```
unknown-map/       CURRENT: unknown-map framing (patrol, sensing range, sensed
                   obstacles, dynamic + track management, sweeps, ablation, N_MC=100,
                   timing; media/ has figure+GIF generation scripts and outputs)
non-obstacle/      Free-space point-stabilization study (oracle / odom / slam)
obstacle-stage-a/  Static obstacle avoidance with a fixed safety margin
obstacle-stage-b/  Static obstacle avoidance with the covariance-aware chance constraint
gamma-sweep/       Static safety vs. efficiency trade-off over the chance factor gamma
dynamic-obstacle/  Moving obstacle: CV-EKF tracker + time-varying chance constraint
ablation-adaptivity/  Margin size vs. adaptivity ablation (matched-mean control)
mc-m100-robustness/   All studies re-run at N_MC=100 (Wilson CIs) as a robustness check
randomized-geometry/  gamma frozen at 2, transferred across 30 randomized feasible
                   geometries (600 trials, disjoint screening/eval seeds, continuous
                   scoring): 0.50% collisions vs 27-30% for fixed margins
filter-consistency/   Robustness to filter mis-calibration: no collisions observed to 4x
                   overconfidence, graceful degradation, vs 14% cov-blind reference
legacy/            Original single-run prototype (kept for reference)
figures/           Figures used in the README / paper
```

Each experiment folder is self-contained (it carries its own copy of `mc_ekf_step.m`).

| Run this | Produces |
|---|---|
| `non-obstacle/run_montecarlo.m` | `mc_results.mat`, tracking-error + trajectory figures |
| `obstacle-stage-a/run_montecarlo_obs.m` | collision-rate + trajectory figures (fixed margin) |
| `obstacle-stage-b/run_montecarlo_obs.m` | same, with `cfg.cov_aware = true` |
| `gamma-sweep/run_gamma_sweep.m` | `gamma_sweep.mat`, trade-off figure |
| `dynamic-obstacle/run_montecarlo_dyn.m` | 4-strategy collision study (oracle/static/cv_fixed/cv_cov) |
| `dynamic-obstacle/run_gamma_sweep_dyn.m` | dynamic safety–efficiency trade-off |
| `dynamic-obstacle/fig_side_by_side.m`, `make_dyn_media.m` | paper figure + animations |
| `ablation-adaptivity/run_ablation.m` | margin size-vs-adaptivity ablation (static + dynamic) |
| `mc-m100-robustness/run_mc100.m` | all 7 studies at N_MC=100 with Wilson 95% CIs (`mc100_results.txt`) |
| `randomized-geometry/run_um_randgeom.m` | gamma-transfer over 30 randomized geometries (`um_randgeom_results.txt`) |
| `randomized-geometry/run_um_randgeom_v2.m` | corrected rerun: disjoint screening/eval seeds + continuous collision scoring (`um_randgeom_v2_results.txt`) |
| `randomized-geometry/clustered_analysis.py` | clustered intervals + geometry-level sign test (the analysis reported in the paper) |
| `filter-consistency/run_um_consistency.m` | filter mis-calibration sweep (`um_consistency_results.txt`) |
| `unknown-map/paired_reanalysis_v2.m` | paired re-analysis of the shared-seed studies (`paired_reanalysis_v2.txt`) |
| `*/time_perf*.m` | per-step timing benchmark |

### Paired analysis of the shared-seed studies

Every mode in the unknown-map studies sees the identical noise realization at a
given trial index, so those studies are paired data and are analysed as such in
`unknown-map/paired_reanalysis_v2.m` (no new simulation; it re-reads the stored
per-trial outcomes). Collisions use exact McNemar tests on discordant pairs;
continuous outcomes use paired bootstrap confidence intervals with **sign-flip
permutation** p-values. Two multiplicity families are declared before the
results are read: the confirmatory continuous outcomes within each study, and
the four collision tests pooled across both studies. All four collision tests
survive Holm correction (adjusted p = 4.5e-13, 1.4e-6, 1.5e-5, 0.016). The
resampling is seeded (`rng(31,'twister')`), so the intervals and p-values
reproduce exactly.

Three points are worth reading before using these numbers:

- Waypoint completion is recomputed from the stored **true** trajectory. The
  logged `wp_reached` field is the executive's stage index, which also advances
  on timeout, so it returns 3 for every mode and must not be used as a
  completion metric. On the corrected metric the dynamic covariance-aware
  controller reaches **fewer** waypoints than `cv_fixed` (-0.48 of 3,
  permutation p = 0.0006, survives Holm) - it trades task progress for
  clearance, and that is reported rather than smoothed over.
- Path length is given over all trials ("simulated path", which includes
  post-collision travel by comparators that collide in up to 44% of trials). It
  is also given over the collision-free subset of pairs, but that row is
  **exploratory and carries no p-value**: conditioning on both controllers
  avoiding a collision selects on an outcome the controller itself affects, so
  the surviving pairs are not exchangeable and a test on them would have no
  protected error rate.
- The randomized-geometry design reuses the same 20 noise realizations across
  all 30 geometries, so it is **crossed** (geometry x seed), not merely
  clustered. Intervals resample geometries and seeds independently. The
  oracle comparison is +0.17 pp, 95% CI [-1.33, +1.67]; this is a failure to
  detect a difference and is **not** evidence of equivalence, since no
  equivalence margin was declared in advance.

### Covariance adequacy study (pre-registered, not yet run)

`coverage-study/PREREGISTRATION_v2.md` is the frozen protocol for the next
experiment, which measures whether the covariance driving the safety margin
actually covers the robot-obstacle relative error. It is committed and tagged
(`prereg-coverage-v2`) **before** the run, and it states in advance what result
would narrow the paper's claims. No coverage results exist in this repository at
that tag.

`coverage-study/PREREGISTRATION.md` is the superseded v1, kept unchanged at tag
`prereg-coverage-v1` for provenance. Review of v1 found that its primary
decision rule called coverage "adequate" whenever the interval merely *contained*
0.95 - the non-significance-as-equivalence error that the TOST in the same
document existed to prevent. v2 replaces it with one-sided non-inferiority
against a declared 0.93 boundary, adds an explicit inconclusive category, ties
the TOST to a 90% interval inside [0.93, 0.97], justifies the sample size from a
power calculation with a pilot-driven re-estimation rule, and specifies
time-point eligibility, numerical handling and the strata definitions. The pause
for review is what surfaced this, before any coverage number existed.

## Requirements

- **MATLAB** (developed on a recent release; uses `wrapToPi`/`wrapTo2Pi`).
- **CasADi** for MATLAB, v3.5.5 (<https://web.casadi.org/>).

Each runnable script begins with an `addpath(...)` to CasADi — **edit that path** to
point at your CasADi install, e.g.:

```matlab
addpath('C:\path\to\casadi-windows-matlabR2016a-v3.5.5');
import casadi.*
```

Then, in MATLAB, `cd` into an experiment folder and run the corresponding
`run_*` script.

## Method summary

- **Unicycle model** with a range–bearing measurement model.
- **EKF-SLAM** with exact velocity-motion prediction and a sequential per-landmark
  measurement update; landmarks are surveyed (known prior, refined online).
- **NMPC** via multiple-shooting transcription, solved with CasADi + IPOPT, with
  actuator-limit and disk obstacle-avoidance constraints.
- **Covariance-aware chance constraint:** the obstacle keep-out margin is
  `delta = delta0 + gamma * sqrt(lambda_max(Sigma_xy))`, where `Sigma_xy` is the
  EKF-SLAM position-covariance block.
- **Dynamic obstacle:** a constant-velocity EKF (`cv_tracker_step.m`) tracks the
  moving obstacle; the horizon constraint uses the predicted track, and the margin
  becomes `delta_k = delta0 + gamma * sqrt(lambda_max(P_xy + Sigma_obs,k))`, combining
  robot-localization and obstacle-prediction uncertainty (the latter growing with
  look-ahead). A soft slack keeps the QP feasible in tight encounters.

## Notation (code &harr; paper)

The code keeps its own self-consistent variable names; they map to the paper symbols
as follows (paper notation follows the standard Kalman convention: **Q** = process,
**R** = measurement):

| Code | Paper | Meaning |
|---|---|---|
| `X` | `x_hat` | EKF-SLAM state estimate (augmented) |
| `Sigma` | `P` | state covariance |
| `P.M` | `Q_u` | input/control-noise covariance `diag(sigma_v^2, sigma_omega^2)` |
| `P.Q` | `R` | measurement-noise covariance `diag(sigma_d^2, sigma_phi^2)` |
| `P.L` | `M` | number of landmarks |
| `P.dt`, `T` | `dt` | sampling period |
| `cfg.lm` | `(l_x,i, l_y,i)` | surveyed landmark coordinates |
| `cfg.xs` | `x*` | reference posture |
| `cfg.gamma` | `gamma` | chance-constraint factor |
| `cfg.safe_buffer` | `delta_0` | fixed keep-out buffer |
| `Q`, `R` (in `run_*`/`mc_build_*`) | `W_x`, `W_u` | MPC stage-cost weights |
| `N` | `N` | MPC prediction horizon |

Note: the parameter struct `P` in `mc_ekf_step.m` is **not** the covariance — the
covariance is `Sigma` (paper `P`).

## Attribution

The MPC/CasADi single-shooting and obstacle-avoidance formulations build on the
open-source NMPC-for-mobile-robots tutorial by M. W. Mehrez
(*Stabilizing NMPC of wheeled mobile robots using open-source real-time software*,
ICAR 2013; and the associated MATLAB workshop materials).

## Citation

If you use this code, please cite the accompanying paper (details to be added upon
publication).

## License

MIT — see [LICENSE](LICENSE).

## Contact

Dun Liu — leoliudun0818@gmail.com
