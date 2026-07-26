# INVALID RESULT — DO NOT CITE

**`um_dyn_ablation_results.txt` (dynamic matched-mean ablation) is invalid.**
The corrected experiment has not yet been run. Until it is, no conclusion may be
drawn from that file or from the dynamic row of Table 8 in the manuscript.

## What is wrong

The file reports:

```
cv_cov (gamma=2)       : ... mean applied inflation 4.4419 m
fixed-matched (4.442 m): coll 0.0% [0.0, 7.1] | min clear +0.7227 | path 16.231
```

The figure **4.4419 m is not the mean applied inflation.** In
`mc_run_trial_um_dyn.m` the variable `infl_now` is assigned only inside the
live-track branch and is never reset when the obstacle track is absent or has
been dropped. The logging line then copies the last surviving value into every
subsequent timestep. Because the tracker's covariance grows monotonically while
coasting, the frozen value is the *peak* of the last coasting episode, and it is
averaged over the ~71.6% of the mission during which no obstacle constraint
existed at all. Roughly 84% of the reported figure is stale carry-over from
steps with no active constraint.

A sanity check that should have caught this: 4.44 m of inflation plus the 0.30 m
keep-out gives a 4.74 m exclusion radius in a 4x4 m workspace — larger than the
entire arena.

**A second, independent defect:** even with the logging repaired, the comparator
was mis-designed. The fixed margin was derived from an *all-mission* mean (steps
with no active constraint counted as contributing) despite `fixed_extra` being
applied only while the fixed arm's own track is active. With a ~28% active
fraction, the arms were therefore not exposure-matched under any reading of the
number. The corrected rerun uses the active-track conditional mean, estimated on
a **calibration seed stream disjoint from the evaluation seeds**
(`DYN_ABLATION_ESTIMAND.md`), and reports the realised exposure of both arms.

## What this invalidates

- the `fixed-matched (4.442 m)` row: its constant, its `0.0%` collision cell and
  its `16.231` path length;
- the manuscript's claim that a size-matched constant is equally safe in the
  **dynamic** scenario;
- the conclusion that "the static conclusion therefore extends", i.e. that
  **margin size rather than adaptivity** explains the safety result. That claim
  is **suspended**, not reversed: no genuinely matched dynamic experiment exists
  yet, in either direction.

## What is NOT affected

Verified against the stored artifacts and unchanged:

- the headline dynamic collision result (88% / 44% / 0% / 0%);
- the static sensed-obstacle results (0% / 36% / 14% / 0%);
- both gamma sweeps;
- the randomized-geometry transfer study;
- the filter-miscalibration sweep;
- the M=100 robustness runs.

The `cv_cov` arm itself was never affected: `infl_t` is instrumentation only and
never reaches the solver. The controller behaved correctly; only the measurement
of what it applied was wrong.

## The static ablation

`fixed-matched (0.146 m)` in the static scenario is **not** affected by the
stale-variable defect, but it shares a design limitation: the constant is
matched to the covariance-aware controller's *all-mission* mean inflation
(including steps with no active constraint), while the constant itself is
applied only while the constraint is active. Measured from the stored data, the
static constraint is active 93.6% of the mission, so the all-mission mean
(0.1461 m) is ~7% below the active-track mean (0.1562 m). The static row is
therefore retained as a **descriptive** result — "a fixed constant set to the
all-mission mean also recorded 0/50 collisions" — and not as a matched-mean
experiment supporting a causal claim.

## Status

A corrected dynamic experiment is being prepared under a frozen estimand and
logging contract. This notice will be replaced by the corrected artifact. The
invalid file is retained rather than deleted so that anyone who already has it
can identify what they have.
