# Frozen experiment protocol — N3

Protocol version: 1.0  
Status: frozen before production runs. SHA-256 is recorded in `checkpoints/N3_protocol.md` and each run output.

## Research questions

1. At equal objective-call budget, what corrected vectorial-NL/DU/LE outcomes do NSGA-II search and uniform parameter random search produce?
2. Does FPA change NSGA-II outcomes under paired seeds and identical initial populations/evaluation budgets?
3. How do the map-based generator's S-box metrics compare with a controlled x-only 1D logistic-map generator and uniform random permutations?

No expected FPA or coupling advantage is assumed. Outcomes remain exploratory absent adequate precision; results are simulation of generator/metric computation, not hardware or complete-cipher measurements.

## Metrics and objective directions

Primary optimizer objectives, ordered in minimization space: `(-LE, -NL(S), DU)`, where NL(S) is the minimum Boolean nonlinearity over 255 nonzero output masks; LE and NL are maximized, DU minimized. Post-hoc: SAC mean, BIC-SAC absolute deviation, LAP probability bias. Report selected gamma, coordinate-NL diagnostic, evaluation count and front size as descriptive only. Do not report entropy as informative, NPCR/UACI, hardware runtime/energy or side-channel data.

The FPA paired comparison endpoints are selected `LE`, vectorial `NL_min`, and `DU` (primary family of three, Holm-adjusted exact paired sign-flip tests). SAC, BIC-SAC deviation and LAP are secondary descriptive outcomes with estimates/intervals; do not make confirmatory claims on them. The optimizer and FPA are judged using selected compromise objectives and hypervolume of the nondominated set.

For compromise selection, objective values are transformed to maximization form `(LE, NL, -DU)`, normalized per nondominated front by min-max scaling with zero-span columns set to zero, and the candidate minimizing Euclidean distance to ideal `(1,1,1)` is selected. This is a deterministic front-relative compromise selector, not a curvature-based knee estimator.

Run output is reproducible from the recorded seed/config and includes source/config hashes, actual evaluation count, method variant, metric-definition identifier, selected candidate and per-evaluation records. The initial 80 candidates and initial-front compromise are retained separately within each NSGA-II run as a within-run reference without extra fitness calls.

Crowding distance is computed from the shared ordered minimization vector `(-LE, -NL, DU)`: sort each objective column ascending, assign infinite distance to boundary items, and add normalized adjacent-neighbour gaps for interior candidates; constant columns contribute zero.

## Methods and equal evaluation budget

A candidate is one valid generator parameter tuple evaluated once by `evaluate_chromosome`. Evaluation count is the actual number of evaluator calls; the current code has no fitness cache. The common comparison budget is 940 actual calls per seed/method, calculated from the current legacy FPA schedule: 80 initial + 20 generations × 40 offspring + floor(20/3) × 10 FPA reevaluations. FPA reevaluations are charged against the budget. Every method must stop at exactly 940 calls; if a budget ends mid-generation, truncate the generated candidate list deterministically and do not exceed the budget.

- NSGA-II-only: population 80, rank/crowding tournaments, uniform per-gene crossover 0.5, Gaussian mutation SD 0.01, elitist survivor selection; stop at 940 calls.
- NSGA-II+FPA: same algorithm, seed and initial population as its paired NSGA-II-only run. Every third generation, perturb up to ten first-front members by FPA; each reevaluation consumes one call from the shared 940-call budget. Do not refill after budget exhaustion.
- Uniform parameter random search: sample phase genes uniformly in `[0,1)`, `mu` uniformly in `[3.57,4)`, `gamma` uniformly in `[0,.5)` using the run seed; evaluate 940 independent candidates and retain the nondominated subset.
- Initial population: report first 80 NSGA-II initial candidates as a within-run reference; do not count this as a separate method or add extra calls.
- 1D comparator: x-only recurrence `x <- (mu*x*(1-x)) mod 1`; sample phase gene `p_x` uniformly `[0,1)` and transform by `x0=sin²(pi*p_x)` to match the 2D phase-to-state rule; sample `mu` uniformly `[3.57,4)` under the matched seed; burn 500 steps; use the resulting x state for each of 255 Fisher–Yates swaps from i=255 down to 1 with `j=floor(x*(i+1))`, advance one logistic step per swap. For optimization baseline, evaluate 940 sampled parameter pairs and retain nondominated S-box metric records. This is a generator-family comparison; no gamma coordinate exists in 1D.
- Random-permutation reference: draw 940 independent uniform permutations of `0..255` per seed as a descriptive distribution, not as an equal-cost optimizer or cipher baseline; use the same metric implementation and report distribution/quantiles, no optimizer superiority inference.

Random seeds: paired FPA/NSGA-II seeds `104729, 130363, 155921, 196613, 228017, 262147, 300007, 350377, 400009, 450001, 500009, 550007, 600011, 650011, 700001, 750019, 800011, 850009, 900001, 950009, 1000003, 1050001, 1100009, 1150003, 1200007, 1250009, 1300003, 1350001, 1400009, 1450003`. Other methods use those same identifiers as independent per-method seeds. This is 30 pairs; no replacement based on outcomes. A pilot timing of 10 2D candidate evaluations averaged 0.045 seconds per candidate; four 940-call search variants × 30 seeds is estimated at 1.41 hours, excluding analysis and random-permutation metrics. Record actual wall time per run; do not reduce the seed set silently.

## Statistical analysis

Primary reporting: all per-seed outcomes, mean, sample SD, paired mean difference, and 95% paired bootstrap percentile interval with 100,000 resamples and fixed RNG seed 20261019. Compute the exact two-sided paired sign-flip p-value by meeting-in-the-middle enumeration: form signed sums for the first and last 15 paired differences; sort one half and count all pairs whose absolute combined sum is at least the observed absolute sum (tolerance 1e-12). This enumerates all `2^30` assignments without storing the full combination matrix. Apply Holm correction to the 3 primary endpoints. Treat CIs/tests as exploratory; report estimates/uncertainty without superiority/equivalence assertions if intervals include zero or precision is weak. Secondary metric outcomes receive descriptive CIs only, not inferential claims.

## Hypervolume

Transform objectives to fixed monotone minimization coordinates: `z_LE = 0.5 - atan(LE)/pi`, `z_NL=(128-NL)/128`, `z_DU=DU/256`. These are bounded in `[0,1)` for finite LE and valid 8-bit permutation NL/DU. Fixed hypervolume reference point is `(1.1,1.1,1.1)`, strictly worse than all transformed fronts. Non-finite values or coordinates not below reference block the run; do not clamp. Compute exact 3D minimization hypervolume per seed and method. Report per-seed values and mean/SD; no hypervolume inferential test is primary.

## Reproducibility and stopping

Record source hash, protocol hash, canonical config hash, seed, variant, evaluation_count, metric definition, Python/NumPy/SciPy versions, start/end timestamps, and output hash for every run. New outputs only under `output/sbox-revision/2026-10-19-inventory-01/results/`; old files immutable. Do not start production unless metric tests pass, pilot schemas validate, every method deterministically reproduces same-seed outputs, all exact-budget assertions pass, and reference point is valid. If 30-seed cost exceeds the available local budget, stop before production and ask for approval of a reduced, explicitly exploratory seed plan; do not silently reduce.
