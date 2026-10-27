# Statistical summary — corrected 30-seed revision

The authoritative analysis is `output/sbox-revision/2026-10-19-inventory-01/analysis/analysis_revision.py`, whose inputs are the immutable valid files under `run-v3/results/`. The frozen run grid completed 30 seeds × five methods (150 cells), each at exactly 940 objective evaluations. Each seed is one unit of analysis; the 940 random permutations inside a seed are nested descriptive samples and are not treated as independent runs.

## Selected outcomes (mean, sample SD; n=30)

| Method | LE | Vectorial NL_min | DU | SAC | BIC-SAC deviation | LAP bias | Hypervolume |
|---|---:|---:|---:|---:|---:|---:|---:|
| NSGA-II-only | 0.678827 (0.006804) | 95.333 (1.322) | 10.000 (0) | 0.502686 (0.004565) | 0.035477 (0.002134) | 0.127604 (0.005163) | 0.722917 (0.006259) |
| NSGA-II + FPA | 0.674967 (0.004968) | 95.733 (0.868) | 10.067 (0.365) | 0.502726 (0.005217) | 0.035131 (0.001751) | 0.126042 (0.003392) | 0.724905 (0.004541) |
| 2D parameter random search | 0.662293 (0.013124) | 95.400 (1.404) | 10.133 (0.507) | 0.499788 (0.004811) | 0.035598 (0.001623) | 0.127344 (0.005486) | 0.718932 (0.006896) |
| 1D Logistic-map random search | 0.662846 (0.017624) | 95.000 (1.145) | 10.067 (0.365) | 0.501278 (0.004653) | 0.035152 (0.001415) | 0.128906 (0.004471) | 0.716873 (0.005934) |

## Primary paired contrast: FPA minus NSGA-II-only

| Metric | Mean difference | Paired bootstrap 95% CI | Exact sign-flip p | Holm-adjusted p | Cohen dz |
|---|---:|---:|---:|---:|---:|
| LE | -0.003860 | [-0.006803, -0.000855] | 0.019117 | 0.057352 | -0.4565 |
| Vectorial NL_min | +0.400 | [-0.133, +1.000] | 0.280273 | 0.560547 | +0.2484 |
| DU | +0.0667 | [0.000, +0.200] | 1.000000 | 1.000000 | +0.1826 |

LE is maximized, so the negative FPA contrast is directionally worse; DU is minimized, so the positive contrast is directionally worse. No primary endpoint has Holm-adjusted p < 0.05. Nonsignificance is not evidence of equivalence.

## Analysis procedures and limits

Intervals are 95% percentile bootstrap intervals over seed-level means or paired differences, using 100,000 resamples with NumPy Generator seed 20261019 (metric-specific deterministic offsets for paired contrasts). The primary paired test is a two-sided exact sign-flip test enumerating all 2^30 assignments using meet-in-the-middle; Holm correction covers the prespecified LE/NL_min/DU family. Hypervolume is reported descriptively under the fixed transforms/reference in `protocol.md`, without a confirmatory test.

Random-search and 1D comparisons are descriptive, not confirmatory matched tests. The random-permutation reference averages 940 permutations within each seed (28,200 total), but seeds—not individual nested permutations—are the independent units. The study is one fixed protocol on a constrained generator family; it establishes neither optimizer superiority nor a coupling mechanism benefit. It contains no complete cipher, hardware, runtime/energy, side-channel, or deployment evaluation. The like-for-like comparison requested with recent chaotic S-box designs (review comment 7(b)) remains unperformed and is acknowledged in the response.
