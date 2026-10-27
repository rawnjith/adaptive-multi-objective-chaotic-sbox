# N5 figure contract

## Figure 1 — Selected objective comparison

Core conclusion: Across the 30 prespecified seeds, the selected LE, standard vectorial NL and DU distributions overlap substantially across NSGA-II, FPA, parameter random search and the 1D comparator.

- Archetype: quantitative grid, per-seed points + mean and bootstrap 95% CI
- Panels: (a) LE, (b) vectorial NL_min, (c) DU
- Evidence hierarchy: paired FPA-vs-NSGA seed points first; random-search and 1D references secondary.
- Statistics: every seed is shown; intervals are seed-level bootstrap CI for means; no unpaired significance claims.
- Reviewer risks: DU is nearly/fully discrete and overlaps; avoid bar-only summaries and avoid hiding seed-level data.

## Figure 2 — FPA paired differences and front-quality

Core conclusion: The paired FPA changes to selected LE/NL/DU are small and uncertain after correcting the nonlinearity objective, while hypervolume describes front quality separately from front size.

- Archetype: asymmetric mixed quantitative figure
- Panels: (a) paired seed-level delta in LE/NL/DU with zero reference; (b) per-seed hypervolume across methods.
- Error bars: bootstrap 95% CI across seed-level observations; paired intervals for contrasts.
- Reviewer risks: identify exact objective transform and HV reference in caption; show all 30 seeds without overlap.

## Figure 3 — Generator controls and random-permutation reference

Core conclusion: Search-based generators and the 1D/random-permutation references occupy overlapping component-metric ranges; no claim of cryptographic strength follows.

- Archetype: distribution comparison, one panel per metric family (NL, DU, LAP)
- Source: seed-level method outputs and 940 random permutations per seed.
- Caption: state n=30 optimizer seeds and n=28,200 random permutations; distributions are descriptive.
- Reviewer risks: permutation samples are nested within 30 seeds, not 28,200 independent replications. Summarize per seed or plot within-seed summaries, do not infer using permutations as independent n.
