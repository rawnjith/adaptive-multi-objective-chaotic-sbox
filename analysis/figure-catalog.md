# Figure Catalog — N5

Figures are diagnostic summaries from the frozen 30-seed simulation grid. Seed is the unit of analysis. Bootstrap bars denote percentile 95% intervals over seeds (or paired seed differences). Random permutations are summarized within seed to avoid treating nested draws as independent.

## Figure QA log

- Renderer source: `make_figures.py`; analysis source hash and results-summary hash recorded in `figures/manifest.json`.
- Exports: vector PDF and editable SVG.
- Automated checks: all source values read from `results_summary.json` / immutable v3 run outputs; no plotted aggregate is hand-entered.
- Visual inspection: rasterized previews saved in `figures/previews/`; inspect those before treating figures as final. Do not modify source/data during this QA pass.

## Figure 1 — Objective distributions

- File: `figures/figure-01-objective-comparison.pdf` (plus SVG).
- Purpose: compare selected LE, vectorial NL_min and DU across four search methods.
- Data: all 30 selected outputs per method in `run-v3/results/*.json`.
- Display: box summaries, all per-seed points, and mean with bootstrap 95% CI.
- Caption requirements: state n=30 seeds/method, mean/CI markers, LE/NL higher-is-better and DU lower-is-better.
- Observation to check: distribution overlap and run variability, not just mean ordering.
- Interpretation boundary: descriptive; no tests across random search or 1D baseline.

## Figure 2 — Paired FPA effect and hypervolume

- File: `figures/figure-02-paired-hv.pdf` (plus SVG).
- Purpose: show paired per-seed differences for FPA vs NSGA-II and per-seed Pareto hypervolume.
- Data: paired differences and HV from the frozen summary.
- Display: all 30 deltas and zero line; per-seed HV distributions with means/95% bootstrap intervals.
- Caption requirements: clarify FPA−NSGA-II direction, bootstrap interval, exact sign-flip and Holm family (three primary outcomes), fixed transforms/reference `(1.1,1.1,1.1)`.
- Observation to check: LE delta distribution vs zero, NL/DU ties or shifts, HV distributions.
- Interpretation boundary: Holm-adjusted LE p=0.0574; no superiority/equivalence conclusion. HV has no inferential test.

## Figure 3 — Component metrics and permutation reference

- File: `figures/figure-03-component-reference.pdf` (plus SVG).
- Purpose: place optimizer selections beside generic uniform random permutation distributions.
- Data: 30 optimizer seed values per method; 940 permutations per seed aggregated to per-seed means for random-permutation reference.
- Display: distributions, seed-level optimizer points and random-permutation seed means.
- Caption requirements: optimizer n=30; permutation n=940 within each of 30 seeds, plotted as seed means; metric definitions/directions.
- Observation to check: whether selected generated boxes separate from generic random-permutation seed means in these component metrics.
- Interpretation boundary: random permutation is only a reference distribution, not a cryptographic comparator; do not use nested sample count as independent n.
