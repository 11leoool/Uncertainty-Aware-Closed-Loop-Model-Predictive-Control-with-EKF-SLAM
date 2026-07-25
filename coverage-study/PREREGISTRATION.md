# Pre-registration: covariance adequacy and coverage study

**Status: FROZEN. No coverage outcome has been computed at the time of writing.**

This document is the analysis plan for the covariance-adequacy experiment. It is
committed and tagged before the instrumented run is executed, so that the
outcome-dependent choices (which levels, which strata, which tests, which
multiplicity family) are fixed in advance rather than selected after seeing the
numbers. The tag that freezes it is recorded in **§10**.

---

## 1. The question

The controller inflates its keep-out radius by a term proportional to the
estimated uncertainty:

```
delta(k) = d0 + gamma * sqrt( lambda_max( P_xy(k) + P_obs(k) ) )
```

with `d0 = 0.06 m` and `gamma = 2.0` (`unknown-map/mc_run_trial_um_obs.m:77`).

Two approximations are embedded here, and this experiment exists to measure
whether they hold empirically:

1. **The relative covariance is approximated as a sum.** The quantity that
   governs collision is the robot-obstacle *relative* position error. Its
   covariance is `P_xy + P_obs - C - C'`, where `C` is the robot-obstacle
   cross-covariance. The implementation uses `P_xy + P_obs` and therefore
   **omits the cross-covariance**. Because the obstacle estimate is derived from
   measurements taken in the robot frame, `C` is not zero, and its omission is
   not conservative by any argument we can currently make. Its sign in practice
   is unknown to us.
2. **The covariance is reduced to a single scalar.** `lambda_max` collapses a
   2x2 covariance to an isotropic radius. The direction that matters for
   collision is the robot-to-obstacle line of sight, whose variance is
   `u' Sigma u <= lambda_max`. This particular step is conservative *given* a
   correct `Sigma`; it is not conservative if `Sigma` itself is optimistic.

**What this experiment tests:** the empirical adequacy of these two
approximations, i.e. whether the covariance the controller actually uses covers
the error it actually incurs, at stated nominal rates.

**What this experiment does not test, and will not be described as testing:**
a formally calibrated chance constraint. No probabilistic safety, feasibility,
or recursive-consistency guarantee follows from this study. A satisfactory
result licenses the statement that the scheduled margin was empirically
well-calibrated in the evaluated conditions, and nothing stronger.

## 2. What is already known (disclosure)

A genuine pre-registration must state prior information rather than imply the
question is being approached blind. Two things are already on record:

- The mis-calibration sweep (`filter-consistency/um_consistency_results.txt`)
  reports, at the consistent setting (filter/true ratio 1.00), a mean claimed
  `sigma` of 0.0672 against a mean **true** localization error of 0.0759. This
  compares a standard deviation to a mean error magnitude and is therefore
  **not** a coverage statement, but it is the reason we consider under-coverage
  a live possibility rather than a formality.
- No relative-position error, cross-covariance, Mahalanobis radius, or coverage
  proportion has been computed in any study to date. The instrumentation
  required to compute them (§5) does not yet exist in the codebase.

Given the first item, we pre-commit to reporting under-coverage as a finding
about the method if it occurs, not as a tuning problem to be resolved by
adjusting `gamma` until coverage looks correct. Any post hoc `gamma` change will
be reported as exploratory and will not be presented as the pre-registered
result.

## 3. Estimands

Let, at logged time `k` of trial `j`:

- `e_rel = (p_rob_true - p_obs_true) - (p_rob_hat - p_obs_hat)` (2-vector)
- `Sigma_hat_rel = P_xy + P_obs` (the approximation actually used)
- `m^2 = e_rel' * inv(Sigma_hat_rel) * e_rel`

Under the working Gaussian assumption `m^2 ~ chi^2_2`, so nominal `1-a` coverage
corresponds to `m^2 <= chi2inv(1-a, 2)`.

**Per-trial coverage** is `c_j(a) = ` fraction of logged samples in trial `j`
with `m^2 <= chi2inv(1-a, 2)`. The estimand is the **mean of `c_j` over trials**.
Time samples are nested within trials and are strongly autocorrelated; they will
**not** be treated as independent observations at any point in this analysis.

### Primary outcome

**95% relative-position coverage**, i.e. mean per-trial `c_j(0.05)`, evaluated
separately in the static and the dynamic scenario.

### Supporting outcomes (reported with CIs; not in the primary family)

- Relative-position coverage at nominal 50% and 90%.
- Robot-pose coverage: same construction with `e_rob` against `P_xy`.
- Obstacle coverage: same construction with `e_obs` against `P_obs`.
- Robot-pose NEES (3-DoF) and measurement NIS, as time-averaged per-trial means.
- **Operating-point exceedance**: fraction of samples with
  `e_rel' * u > gamma * sqrt(lambda_max(Sigma_hat_rel))`, where `u` is the unit
  robot-to-obstacle direction. Under a correct covariance this is bounded above
  by `P(N(0,1) > 2) = 2.28%`, since `u' Sigma u <= lambda_max`. This is the
  single number that most directly describes the deployed margin, and the bound
  is one-sided and conservative by construction.
- The omitted cross-covariance `C` itself, estimated from the paired errors
  across trials, reported as a correlation. This quantifies approximation (1)
  directly and is the diagnostic we expect to be most informative if coverage
  fails.

### Pre-specified strata (descriptive; each reported with trial-level CIs)

- Landmark-rich versus landmark-poor phase, split at the median number of
  landmarks in view across the pooled evaluation trials.
- Obstacle visible versus coasting (no current detection).

Strata are **not** added to the primary family and no stratum will be promoted to
a headline result. They exist to characterise where any failure is located.

## 4. Hypotheses and decision rules

For each scenario, let `C95` be the mean per-trial 95% relative coverage with a
trial-level bootstrap 95% CI.

- **Adequate**: the CI contains 0.95, or lies entirely above it.
- **Under-covering**: the CI lies entirely below 0.95. The margin is optimistic;
  the manuscript must then describe `gamma = 2.0` as an empirically chosen
  inflation factor and must not describe the covariance as calibrated.
- **Over-covering**: CI entirely above 0.95. Reported as conservatism, with the
  path-length cost already measured in the paired analysis as its price.

**Equivalence.** Because "no detected difference" is not evidence of adequacy, we
pre-declare an equivalence margin of **+/- 2 percentage points** on coverage
(i.e. 0.93 to 0.97) and will run a TOST at that margin. The margin is chosen as
the largest deviation we would still be willing to call well-calibrated for a
simulation study of this size; it is declared here, before any coverage number
exists, precisely so it cannot be chosen to fit the result.

**Multiplicity.** The primary family is **m = 2**: 95% relative-position coverage
in the static scenario and in the dynamic scenario. Holm correction is applied
within that family. Supporting outcomes and strata are reported with confidence
intervals and are explicitly labelled as not multiplicity-controlled.

## 5. Instrumentation to be collected

Per time step, per trial:

- true and estimated robot pose; `P` (full), and `P_xy`
- true and estimated obstacle position; `P_obs`
- relative position error and `Sigma_hat_rel = P_xy + P_obs`
- innovations and innovation covariance `S` (for NIS)
- number of landmarks in view
- detection / coasting / reacquisition / track-drop state
- covariance conditioning: eigenvalues of `P_xy`, `P_obs`, `Sigma_hat_rel`
- MPC slack (nonzero frequency and magnitude), solver status, solve time,
  deadline outcome, and the exact fallback command applied on failure

The controller-failure accounting is collected in the same run but analysed and
reported separately; it is descriptive and carries no hypothesis test.

## 6. Design and seeds

- Scenarios: static sensed obstacle, and dynamic (constant-velocity) obstacle.
- **100 evaluation trials per scenario**, paired across modes by shared seed, as
  in the existing studies.
- **Evaluation seed stream: `rng(70125)`.** This is disjoint from every seed used
  in any previous study in this repository (the `rng(2024)` family) and from the
  pilot stream below.
- **Pilot / debugging seed stream: `rng(99001)`, excluded from all frozen
  analysis.** Instrumentation will be debugged on the pilot stream only. If any
  evaluation trial is inspected during debugging, the entire evaluation stream is
  discarded and re-drawn with a new documented seed, and that event is recorded
  in §9.

## 7. Analysis procedure

1. Compute per-trial coverage proportions; never pool raw time samples.
2. Cluster bootstrap over **trials**, B = 20000, seed 11, resampling whole
   trials with their full time series intact.
3. Report the design effect (ratio of the trial-cluster variance to the naive
   per-sample variance) so the cost of the nesting is visible.
4. Holm-correct the primary family (m = 2).
5. Run the pre-declared TOST at +/- 2 pp.
6. Report all strata and supporting outcomes with CIs, unadjusted, labelled.

## 8. What would falsify the current framing

The manuscript currently justifies the scheduled margin on the grounds that the
covariance represents the relevant uncertainty. That framing does not survive if:

- 95% relative coverage is materially below 0.95 in either scenario, **or**
- operating-point exceedance materially exceeds 2.28%, **or**
- the omitted cross-covariance is large enough that `P_xy + P_obs` is optimistic
  rather than conservative in the line-of-sight direction.

In any of those cases the reported contribution narrows to: an uncertainty-scaled
margin that reduced collisions in simulation, with an explicitly approximate and
empirically under-calibrated uncertainty model. We commit to that narrowing here,
in advance, so that the result is not renegotiated after the fact.

## 9. Deviations

Any departure from this protocol is appended below with its date and reason,
before the affected analysis is run. An empty section means the protocol was
followed as written.

*(none)*

## 10. Freeze record

This protocol is frozen by the git tag `prereg-coverage-v1`. The tagged commit
contains this file and no coverage results. Any coverage artifact in this
repository must post-date that tag; if it does not, the pre-registration claim
for that artifact is void.
