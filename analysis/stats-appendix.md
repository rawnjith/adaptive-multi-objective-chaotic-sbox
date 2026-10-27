# Statistics Appendix — N5

## Design and analysis unit

- 30 fixed seeds; all 150 method/seed cells completed at 940 actual candidate/metric calls.
- Unit: one optimizer run/seed. FPA and NSGA-II share identical initial populations; other methods use seed-matched identifiers with method-specific RNG streams.
- Methods: NSGA-II-only, NSGA-II+FPA, 2D uniform parameter random search, 1D x-only logistic parameter random search, and 940 uniform random permutations per seed (descriptive reference only).
- Primary endpoint family for FPA vs NSGA-II: selected largest LE (maximize), vectorial NL_min over 255 nonzero output masks (maximize), and DU (minimize). Holm adjustment is within this 3-endpoint family.
- Secondary component metrics: SAC, BIC-SAC deviation, LAP bias. Descriptive only; no inferential tests claimed.

## Descriptive selected outcomes (mean, sample SD; n=30)

| Method | LE | Vectorial NL_min | DU | SAC | BIC-SAC dev. | LAP bias | Hypervolume |
|---|---:|---:|---:|---:|---:|---:|---:|
| NSGA-II-only | 0.678827 (0.006804) | 95.333 (1.322) | 10.000 (0) | 0.502686 (0.004565) | 0.035477 (0.002134) | 0.127604 (0.005163) | 0.722917 (0.006259) |
| NSGA-II + FPA | 0.674967 (0.004968) | 95.733 (0.868) | 10.067 (0.365) | 0.502726 (0.005217) | 0.035131 (0.001751) | 0.126042 (0.003392) | 0.724905 (0.004541) |
| 2D parameter random search | 0.662293 (0.013124) | 95.400 (1.404) | 10.133 (0.507) | 0.499788 (0.004811) | 0.035598 (0.001623) | 0.127344 (0.005486) | 0.718932 (0.006896) |
| 1D logistic random search | 0.662846 (0.017624) | 95.000 (1.145) | 10.067 (0.365) | 0.501278 (0.004653) | 0.035152 (0.001415) | 0.128906 (0.004471) | 0.716873 (0.005934) |

Bootstrap intervals for each method and metric are in `results_summary.json`.

## Primary paired contrast: FPA minus NSGA-II-only

| Metric | Mean difference | Paired bootstrap 95% CI | Exact sign-flip p | Holm-adjusted p | Cohen dz |
|---|---:|---:|---:|---:|---:|
| LE | -0.003860 | [-0.006803, -0.000855] | 0.019117 | 0.057352 | -0.4565 |
| Vectorial NL_min | +0.400 | [-0.133, +1.000] | 0.280273 | 0.560547 | +0.2484 |
| DU | +0.0667 | [0.000, +0.200] | 1.000000 | 1.000000 | +0.1826 |

Direction matters: lower LE is worse because the optimizer maximizes LE; higher DU is worse because DU is minimized. The raw LE p-value is below .05, but not after the prespecified Holm adjustment. This is not evidence of equivalence or no effect.

## Random-search and 1D descriptive contrasts

Paired seed identifiers are shared, but these methods do not share a matched initial population or a one-to-one candidate pairing with NSGA-II. The `results_summary.json` therefore labels these contrasts descriptive-only and gives per-seed differences and effect sizes but no sign-flip p-values or confidence intervals. Do not present these as confirmatory hypothesis tests.

## Random-permutation reference

There are 940 permutations within each of 30 seeds (28,200 total samples). They are nested within seed and summarize a reference distribution; the 28,200 draws are not treated as 28,200 independent experimental units. Aggregate per-seed permutation statistics if making uncertainty statements.

- Vectorial NL_min: mean 92.6813 across samples, range 80–98.
- DU: mean 11.3438, range 10–20.
- LAP probability bias: mean 0.137964, range 0.117188–0.1875.

## Statistical procedures

- Descriptives: mean and sample SD across seeds.
- Intervals: 95% percentile bootstrap for seed-level means and paired mean differences; 100,000 resamples, NumPy Generator seed 20261019 (paired contrasts use deterministic metric-specific offset seeds).
- Primary test: exact two-sided sign-flip test, meeting-in-the-middle enumeration of all 2^30 assignments with tolerance 1e-12.
- Multiplicity: Holm correction across the three prespecified FPA primary endpoints.
- Effect size: paired Cohen dz (mean paired difference / sample SD of paired differences); if SD is zero, implementation returns zero for zero difference and signed infinity otherwise.
- Assumption check: Shapiro-Wilk on paired LE, NL_min and DU differences; interpreted diagnostically, not as a gate to switch analyses post hoc. LE differences were not rejected as non-normal at conventional levels (W=0.9578, p=0.2724); NL and DU differences were strongly non-normal (W=0.7763, p=2.49e-5; W=0.1796, p=7.77e-12). The exact sign-flip test does not assume independent groups; seeds are the independent units, and pairing is by shared initial population.

## Limitations

- A single fixed 30-seed design estimates this implementation/protocol, not a universal optimizer property.
- Hypervolume depends on the prespecified monotone transforms and fixed reference point; it is descriptive, not a primary hypothesis test.
- 1D and 2D searches differ in dimensionality and cannot alone establish a mechanistic coupling benefit.
- Random permutations are not a cipher or matched optimizer baseline.
- No complete cipher, device runtime/energy, side-channel analysis, or deployment evaluation exists.
- Existing `depositfinal_*` outputs use coordinate-only NL and are excluded from N5 inference.
