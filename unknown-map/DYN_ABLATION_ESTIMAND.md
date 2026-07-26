# Frozen estimand and logging contract — dynamic matched-margin ablation

**Status: FROZEN before the corrected rerun. No corrected result exists yet.**

Resetting the stale variable is not sufficient. The experiment also needs a
definition of what "matched" means, fixed before the run, because the obvious
choices disagree by a factor of about four.

## 1. The estimand

**Target: the covariance-aware reference arm's mean node-0 applied inflation,
conditional on the obstacle constraint being active, computed per trial and then
averaged with equal trial weight.**

- *Conditional on active*, because the constant `fixed_extra` is itself applied
  only while the fixed arm's own track is active. Matching an all-mission mean
  (which includes inactive zeros) to a conditionally-applied constant
  systematically under-matches the fixed arm: with a ~28% active fraction, the
  fixed arm would receive roughly a quarter of the mission-integrated exposure.
- *Node 0*, because that is the inflation actually applied at the constraint the
  controller acts on this step.
- *Per trial, then averaged*, with equal weight per trial rather than pooling
  over all active timesteps. Collision outcomes are trial-weighted, so the
  matching quantity should be too; pooling would let long trials dominate.

## 2. Naming — and what cannot be claimed

The comparator is to be described as:

> a fixed margin **matched to the covariance-aware reference arm's active-track
> mean**

It must **not** be described as exposure-matched, or as applying the same
average keep-out. It cannot be: the margin changes the trajectory, which changes
detection and tracking duration, which changes exposure. Matching is a fixed
point that cannot be reached by choosing a constant in advance.

The residual mismatch is therefore **measured and reported**, not assumed away.

## 2b. Calibration/evaluation separation

The matching constant `c` must not be estimated from the same trials used to
evaluate the two arms — that is a mild form of the seed leakage that produced
the retracted randomized-geometry zero.

1. **Calibration stream:** run the covariance-aware reference arm on a fresh,
   disjoint seed stream (`rng(81000 + t)`, t = 1..50, disjoint from every seed
   family used in this repository) and compute
   `c = mean over trials of (per-trial mean of infl0_t | active)`.
2. **Freeze `c`** — recorded in this file before any evaluation trial runs.
3. **Evaluate** both arms on the original paired evaluation seeds, untouched by
   calibration.

Project-wide rule, adopted from this defect onward:

> Any parameter estimated from simulated data must be frozen using a
> calibration stream disjoint from the evaluation stream.

- Frozen `c`: *pending calibration run*

## 3. Logging contract

The rerun must record, per trial and per timestep:

- `active_t` — whether an obstacle constraint existed at MPC construction time
- track status: not-yet-detected / detected / coasting / dropped
- `infl0_t` — the exact node-0 inflation applied, or **NaN when inactive**
  (never a carried-over value)
- `P_xy` and `P_obs` exactly as used to build the constraint, at construction
  time — not the post-update covariances
- MPC slack: frequency and magnitude per step

Assertion, to run in every trial:

```matlab
assert(all(isnan(infl0_t(~active_t))), 'inflation recorded while inactive');
assert(mean_inflation_all <= active_frac * peak_inflation + 1e-12);
```

The second is a pure arithmetic invariant: an all-mission mean cannot exceed the
active fraction times the peak. It would have caught the original defect with no
instrumentation at all.

## 4. Reported outputs (all of them, both arms)

Fixed in advance so none can be selected after the fact:

- target reference-arm conditional mean (the matching quantity)
- **realised** conditional mean, both arms
- active-track fraction, both arms
- mission-integrated inflation, both arms
- median, IQR and maximum inflation, both arms — not the mean alone
- slack frequency and magnitude. A nominal ~2.6 m inflation in a 4x4 m
  workspace is only interpretable through the soft constraint's realised slack;
  without it the number is not physically meaningful.
- paired collision discordances, exact two-sided McNemar, risk difference and
  interval

## 5. Interpretation rules, declared now

- With only a handful of discordant pairs this study is **not powered** to
  separate the arms. A 2-versus-0 discordance gives exact two-sided McNemar
  **p = 0.50**. Such an outcome is **inconclusive** and will be reported as
  inconclusive — not as evidence that either controller is safer.
- The global claim that **margin size rather than adaptivity** explains the
  safety result is **suspended** until a genuinely matched experiment exists. It
  is not replaced by its converse.
- If the realised conditional means of the two arms differ by more than 10%, the
  experiment is reported as **approximately matched, with the direction and size
  of the mismatch stated**, and no causal conclusion is drawn from it.

## 6. Note on a better experiment

The matched-margin design may be structurally unable to answer the question it
was built for, for the reason in section 2: exposure cannot be held fixed while
adaptivity is varied. A design that avoids the problem entirely is a
**safety-versus-efficiency frontier**: sweep the fixed margin over a range,
plot (path length, collision rate) for each constant, and overlay the
covariance-aware operating point.

- If the covariance-aware point lies **on** the frontier traced by constants,
  adaptivity buys convenience only — it finds a good constant without tuning.
- If it lies **inside** the frontier — fewer collisions at equal path, or
  shorter path at equal safety — adaptivity buys something a constant cannot.

This requires no matching, so it cannot be undermined by the trajectory-feedback
problem. The machinery already exists (`cfg.fixed_extra`, and the gamma-sweep
harness). It is recorded here as the recommended successor experiment, not as
part of this rerun.
