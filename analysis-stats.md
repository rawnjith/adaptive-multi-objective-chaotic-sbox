# Results / Statistics Appendix — authoritative paired_statistics.json

## Experimental design
Paired NSGA-II versus NSGA-II + periodic FPA runs. Both variants use seeds 42, 100, 1234, 2025, 7777; population 80; 20 generations; Gaussian mutation SD 0.01; phase genes `p_x,p_y∈[0,1)` mapped by `x0=sin²(πp_x)`, `y0=sin²(πp_y)`; `mu∈[3.57,4.0)`, `gamma∈[0,0.5)`. FPA perturbation affects up to ten members from the first nondominated front every three generations. The authoritative outputs are `output/reproducibility-deposit/results/depositfinal_{fpa,nsga_only}_seedN.json`.

## Descriptive and paired results
Values are mean (sample SD), n=5. Delta = FPA minus NSGA-II-only. Percentile bootstrap CI is from 100,000 resamples of five paired differences using NumPy Generator seed 20261007. The test is an exact two-sided paired sign-flip test over 32 sign assignments; Holm correction is across six endpoints.

| Metric | FPA mean (SD) | NSGA-II-only mean (SD) | Paired delta | Bootstrap 95% CI | Exact p | Holm p |
|---|---:|---:|---:|---:|---:|---:|
| LE | 0.672463 (0.002947) | 0.675409 (0.007806) | -0.002947 | [-0.008543, 0.003321] | 0.4375 | 1.000 |
| NLmin | 101.6 (0.894) | 102.4 (0.894) | -0.8 | [-1.6, 0.0] | 0.5000 | 1.000 |
| SAC | 0.500635 (0.003447) | 0.498242 (0.003991) | 0.002393 | [-0.001074, 0.006055] | 0.3750 | 1.000 |
| BIC-SAC deviation | 0.034180 (0.000789) | 0.034849 (0.001615) | -0.000670 | [-0.001995, 0.000656] | 0.4375 | 1.000 |
| DU | 10.0 (0.0) | 10.0 (0.0) | 0.0 | [0.0, 0.0] | 1.0000 | 1.000 |
| LAP probability bias | 0.142188 (0.006536) | 0.135938 (0.008908) | 0.006250 | [-0.004688, 0.017188] | 0.5000 | 1.000 |

No Holm-adjusted test supports a difference. Five pairs are underpowered; the tests do not show equivalence. The manuscript makes no FPA superiority claim.

Final nondominated-front sizes are FPA [6,7,7,4,6] (mean 6.0) and NSGA-II-only [4,5,3,6,4] (mean 4.4). This count is descriptive only and is not a measure of front quality. The ten selected gamma values span approximately 0.00045 to 0.49159; the selected solutions include near-zero coupling, so 1D-versus-2D comparison remains important.

## Implementation corrections
- Benettin tangent-vector LE computation now accumulates `log(norm)` before renormalization; uncoupled mu=4 sanity estimate approaches `ln(2)`.
- The modulo-one coupled map uses simultaneous old-state updates consistently in trajectory generation, LE estimation, and Fisher–Yates key generation.
- Optimizer phase genes map to valid initial states through `sin²(πp)`; phase, mu, and gamma genes wrap periodically.
- Seed and no-FPA CLI behavior fixed and provenance records variant and seed.
- Parent selection uses binary tournament by NSGA-II rank then crowding; FPA perturbation acts only on the first nondominated front.
- Fitness has three independent objectives (LE, NLmin, DU); runtime proxy removed. LAP now returns probability bias, not absolute correlation.
- AES S-box generated and checked against the known first four entries, permutation property, NLmin=112, DU=4, and LAP=0.0625.

## Limits and external blockers
No complete cipher/key schedule is specified, so NPCR/UACI are omitted. No physical MCU time/energy/side-channel data exist. No matched baselines other than AES 8-bit are included; the AES row is not a cipher-security comparison. A full archive-quality metric, larger-seed study, and 1D baseline remain future work. Local code/data release is ready but has no public DOI; the Data Availability statement must be replaced after deposit. Full visual PDF inspection and final response-letter reconciliation remain necessary.
