# Pre-registration v2: covariance adequacy and coverage study

**Status: FROZEN, superseding v1. No coverage outcome has been computed.**

This amends `PREREGISTRATION.md` (tagged `prereg-coverage-v1`), which is
preserved unchanged for provenance. Nothing here was written with knowledge of
any coverage result: the instrumentation of §5 does not yet exist, and no
relative-position error, Mahalanobis radius or coverage proportion has been
computed in any study to date.

## 0. Why v2 exists

Review of v1 found a decision rule that contradicted the rest of the document.
v1 §4 declared coverage **adequate** whenever the confidence interval *contained*
0.95. That is the non-significance-as-equivalence error - the very inference the
TOST in the next paragraph existed to prevent. A wide, uninformative interval
would have been scored as success.

v2 replaces that rule. It also fixes seven smaller gaps: the boundary for
inadequacy was tied to 0.95 rather than to the declared margin; there was no
inconclusive category; the primary question was left ambiguous between
equivalence and non-inferiority; the sample size was asserted rather than
justified; time-point eligibility was undefined; several implementation details
(controller mode, singular covariances, heading wrapping, exact test procedure)
were unspecified; and the landmark stratum used a data-dependent median split.

Sections 1-3 of v1 (the question, the disclosure of prior information, and the
estimands) stand unchanged and are not restated in full here.

## 1. The primary question, decided

v1 left open whether the primary claim is two-sided calibration equivalence or
one-sided non-inferiority. **The primary claim is one-sided non-inferiority
against 0.93.**

The rationale is what the covariance is used for. The margin
`delta = d0 + gamma*sqrt(lambda_max(P_xy + P_obs))` exists to keep the robot
away from the obstacle. Under-coverage means the true relative error exceeds
what the margin budgeted for, which is a safety failure. Over-coverage means the
margin is larger than it needed to be, which costs path length - a real cost,
already measured in the paired analysis, but not a safety failure. The two
directions are therefore not symmetric in consequence, and a two-sided
equivalence test would treat them as if they were.

Two-sided calibration equivalence is retained as a **secondary** outcome,
because "the covariance is well calibrated" is a stronger and more interesting
statement than "it is not dangerously optimistic", and the paper should be able
to say which one it has earned.

## 2. Decision rules

Let `c_j` be the per-trial 95% relative-position coverage proportion (§3 of v1)
and let `C = mean(c_j)` over trials, with trial-level cluster bootstrap
intervals (B = 20000, seed 11).

**Primary - non-inferiority to 0.93** (one-sided, alpha = 0.05):

| Outcome | Rule |
|---|---|
| **Adequate** | lower limit of the one-sided 95% bootstrap interval > 0.93 |
| **Materially inadequate** | upper limit of the one-sided 95% interval < 0.93 |
| **Inconclusive** | the interval spans 0.93 |

"Materially inadequate" is defined relative to the declared 0.93 boundary, not
to whether an interval happens to sit below 0.95. A point estimate of 0.94 with
an interval inside [0.93, 0.97] is *not* a failure; it is a well-calibrated
margin with a small optimistic bias.

**Secondary - calibration equivalence** (TOST at +/- 2 pp): equivalence is
declared only if the **90% two-sided bootstrap interval falls entirely within
[0.93, 0.97]**. This is the interval-inclusion form of TOST at alpha = 0.05 per
side and is exactly equivalent to requiring both one-sided tests to reject.
Reporting the 90% interval, not the 95% one, is what makes the two procedures
agree; v1 did not say which interval to use.

**Tertiary - conservatism**: if the 90% interval lies entirely above 0.97, the
margin is reported as conservative, with the path-length cost quoted from the
paired analysis as its price.

An **inconclusive** primary result will be reported as inconclusive. It will not
be described as adequacy, and it will not be converted into a claim by
re-running with a different `gamma`; any such re-run is exploratory and labelled
so.

## 3. Sample size, justified

For the TOST secondary outcome with margin `delta = 0.02` and true difference
zero, the required trial count is approximately

```
n  >=  (z_{1-alpha} + z_{1-beta})^2 * sigma^2 / delta^2
```

where `sigma` is the **between-trial** SD of `c_j`. At alpha = 0.05, 80% power:

| Formula | Constant | n = 100 supports sigma up to |
|---|---|---|
| `(z_0.95 + z_0.84)^2` = 6.18 | standard TOST approximation | **0.080** |
| `(z_0.95 + z_0.90)^2` = 8.57 | conservative variant | **0.068** |

So 100 trials per scenario is adequate **only if** the between-trial coverage SD
is about 0.08 or less. That is an assumption, not a fact, and it is exactly the
quantity the pilot exists to measure.

**Sample-size rule, declared now:** the pilot stream (§6) will report the
observed `sigma_hat` of `c_j`. The evaluation `n` is then set to

```
n_eval  =  max( 100,  ceil( 6.18 * sigma_hat^2 / 0.02^2 ) )
```

rounded up to the next multiple of 25, capped at 400 per scenario. This number
is **frozen and recorded in §9 before any evaluation trial is run**, and is
computed from pilot data only. The non-inferiority primary is less demanding
than the TOST, so a sample size adequate for TOST is adequate for it.

If `sigma_hat` implies `n_eval` above the cap, the study proceeds at 400 and the
secondary TOST is reported as underpowered by design, with the achieved power
stated.

## 4. Which time points are eligible

v1 did not define this, which would have let the analysis silently mix periods
when no obstacle estimate existed. A sample `(j, k)` enters the primary analysis
if and only if, at that step:

1. the obstacle track has been **initialised** (a detection has occurred and
   `P_obs` has been set from a measurement, not from its prior), **and**
2. the track has **not been dropped** (the staleness timer has not fired), **and**
3. `Sigma_hat_rel` is numerically usable by the rule in §5.

Everything else is excluded from the primary estimand and **accounted for
separately**, not silently discarded. The following are reported per scenario as
descriptive quantities:

- fraction of steps before first detection (no obstacle estimate exists);
- fraction of steps in coasting (track alive, no current detection);
- fraction of steps after a track drop;
- fraction excluded by the numerical rule in §5.

Coasting steps are **inside** the primary analysis - they are precisely where an
obstacle covariance is most likely to be optimistic - and are additionally
reported as their own stratum (§7).

## 5. Implementation details, specified

**Primary controller mode.** The covariance-aware controller in each scenario:
`slam_cov` for the static sensed obstacle, `cv_cov` for the dynamic obstacle.
These are the modes whose keep-out margin is actually driven by
`P_xy + P_obs`. Comparator modes are not part of the coverage estimand; a
fixed-margin controller has no covariance to assess.

**Singular and ill-conditioned covariance.** Before inverting `Sigma_hat_rel`,
compute its eigenvalues. A sample is excluded from the primary analysis if the
smaller eigenvalue is `< 1e-12` m^2 or the condition number exceeds `1e8`. The
inverse is computed by symmetric eigendecomposition rather than `inv`, and
`Sigma_hat_rel` is symmetrised as `(S + S')/2` first. Excluded counts are
reported per scenario (§4); if exclusions exceed 1% of eligible samples in
either scenario, the numerical rule and the covariance conditioning are reported
as a finding in their own right, since that would itself indicate a filter
problem.

**Heading.** The relative-position estimand is planar and does not involve
heading. Heading enters only the supporting robot-pose NEES, which is 3-DoF; the
heading component of the error is wrapped to `(-pi, pi]` before forming the
quadratic form. Failing to wrap would produce spurious `2*pi` errors at the
wrap boundary and inflate NEES.

**TOST.** Implemented as interval inclusion: compute the 90% two-sided
cluster-bootstrap percentile interval for `C` and declare equivalence iff it
lies within [0.93, 0.97]. No normal approximation is used.

**Holm.** Applied within the primary family only (m = 2: the static and dynamic
non-inferiority tests). Implementation: sort the two p-values ascending, adjust
as `p_(i)_adj = max over k <= i of ( (m - k + 1) * p_(k) )`, cap at 1. Supporting
outcomes and strata are reported unadjusted and labelled as such.

**Bootstrap.** Whole trials are resampled with replacement, each carrying its
full time series; per-trial `c_j` is recomputed from the resampled set. Time
samples are never resampled independently. B = 20000, seed 11. The design effect
relative to a naive per-sample interval is reported.

## 6. Seeds

Unchanged from v1: evaluation stream `rng(70125)`, disjoint from every previous
study (the `rng(2024)` family); pilot stream `rng(99001)`, excluded from all
frozen analysis. The pilot is used for instrumentation debugging **and** for the
`sigma_hat` that sets `n_eval` (§3). If any evaluation trial is inspected before
the sample size is frozen, the evaluation stream is discarded and re-drawn with
a new documented seed, recorded in §9.

## 7. Strata

Reported descriptively with trial-level intervals, unadjusted, not promoted to
headline results:

- **Landmark visibility: 0-1 versus 2-3 landmarks in view.** This replaces v1's
  median split, which was data-dependent and therefore not a pre-registrable
  definition. The fixed cut is chosen because two range-bearing landmarks are
  the minimum for a well-constrained planar pose, so the boundary corresponds to
  a change in observability rather than to a quantile of whatever the run
  produced.
- **Obstacle visible versus coasting**, as defined in §4.

## 8. What would falsify the current framing

Unchanged from v1 in substance, restated against the corrected boundaries. The
manuscript's justification of the scheduled margin does not survive if:

- the primary non-inferiority test returns **materially inadequate** in either
  scenario, or
- operating-point exceedance materially exceeds 2.28%, or
- the omitted robot-obstacle cross-covariance is large enough that
  `P_xy + P_obs` is optimistic in the line-of-sight direction.

In any of those cases the contribution narrows to: an uncertainty-scaled margin
that reduced collisions in simulation, with an explicitly approximate and
empirically under-calibrated uncertainty model. An **inconclusive** result
narrows it differently: the margin's calibration would be reported as untested
at the achieved precision, which is a weaker statement than either adequacy or
inadequacy and must not be written as the former.

## 9. Deviations and frozen quantities

To be completed **before** the corresponding step runs.

- Pilot `sigma_hat` (static): *pending*
- Pilot `sigma_hat` (dynamic): *pending*
- Frozen `n_eval` (static): *pending*
- Frozen `n_eval` (dynamic): *pending*
- Deviations: *none*

## 10. Freeze record

v1 is frozen at tag `prereg-coverage-v1` and is **not** rewritten. This
amendment is frozen at tag `prereg-coverage-v2`. Both tags contain this
directory with no coverage results. Any coverage artifact in this repository
must post-date `prereg-coverage-v2`; if it does not, the pre-registration claim
for that artifact is void.
