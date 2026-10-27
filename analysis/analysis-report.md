# Analysis Report — Corrected Chaotic S-box Study

Status: evidence summary from N4 v3. All values are simulated optimizer/generator evaluations; no hardware or whole-cipher results are implied.

## Question

With the standard vectorial S-box nonlinearity objective (minimum over 255 nonzero output masks), how do NSGA-II, NSGA-II+FPA, equal-budget uniform parameter search, and a controlled 1D logistic-map search compare? Is any observed FPA change supported after pairing by seed and adjusting the three prespecified primary endpoints?

## Data and unit of analysis

The frozen grid completed 30 seeds × five methods at exactly 940 evaluation calls per method/seed. The analysis uses the seed as the unit of analysis. FPA and NSGA-II share the initial population for each seed. Uniform random-permutation references contain 940 draws nested within each seed; those draws are summarized descriptively within seed and are not treated as 28,200 independent inferential observations.

The run grid passed `verify_run_grid.py`: 150/150 completed, no failed/pending entries, protocol hash `662824ab0364f19f9f7e61ca7ddd69ad568780dbd38e5305786b4625a4e8d551`, simulator hash `04775ce3631504457859ed6a8ec0a79e837406785899b4ee62239f31469683cf`. The verifier recomputed selected vectorial NL from S-box tables, checked bijections, hashes, call counts, front hypervolume, and paired initial populations.

## Main observations

Across the 30 selected solutions, NSGA-II-only had mean vectorial NL 95.33 (sample SD 1.32), DU 10.00 (SD 0), and LE 0.67883 (SD 0.00680). NSGA-II+FPA had mean vectorial NL 95.73 (SD 0.87), DU 10.07 (SD 0.37), and LE 0.67497 (SD 0.00497). The random-search means were NL 95.40, DU 10.13, LE 0.66229. The 1D search means were NL 95.00, DU 10.07, LE 0.66285. Full means, SDs and bootstrap intervals are in `results_summary.json`.

For the paired FPA-minus-NSGA-II contrast, mean LE difference was −0.00386, percentile paired-bootstrap 95% CI [−0.00680, −0.00085], exact sign-flip p=0.0191 and Holm-adjusted p=0.0574 across LE/NL/DU. Cohen's dz was −0.456. Vectorial NL difference was +0.40 (95% CI [−0.13, 1.00], Holm p=0.561); DU difference was +0.067 (95% CI [0.00, 0.20], Holm p=1.0). Thus the prespecified FPA family contains no Holm-adjusted p<0.05 result. The LE direction is lower under FPA, not an improvement if LE is to be maximized. No superiority or equivalence claim is warranted.

Mean hypervolume under the fixed protocol transform was 0.72292 (SD 0.00626) for NSGA-II-only, 0.72491 (SD 0.00454) for FPA, 0.71893 (SD 0.00690) for 2D random parameter search, and 0.71687 (SD 0.00593) for 1D search. Hypervolume is descriptive; no confirmatory inferential test was prespecified for it. Front size is also descriptive.

The 28,200 uniform random permutations averaged vectorial NL 92.68 (range 80–98), DU 11.34 (range 10–20), and LAP bias 0.13796 (range 0.11719–0.1875). Since these are nested samples by seed and not a matched optimizer, comparisons are distributional context only.

The 1D comparator has gamma=0 by definition. Its selected LE mean was 0.66285 and vectorial NL mean 95.0; the 2D random search means were LE 0.66229 and vectorial NL 95.4. This does not show a clear coupling benefit; the comparison is exploratory and its construction differs in available parameter dimensions.

## Supported interpretation and limits

The corrected objective produces substantially lower vectorial NL than the former coordinate-only NL headline (the old ten boxes recomputed at 90–96, while the new 30-seed NSGA/FPA means are about 95). The old results are invalid as evidence for optimizing standard vectorial NL and remain quarantined.

Within this protocol, FPA did not demonstrate a statistically supported improvement over NSGA-II on the prespecified primary endpoints after Holm correction. Search-vs-random and 2D-vs-1D contrasts are reported descriptively with seed-level differences/effect sizes and bootstrap intervals; they are not confirmatory tests. None of these comparisons establishes cryptographic security, whole-cipher security, deployment readiness, runtime or side-channel properties.

## Claim candidates

- Claim: The implementation optimizes standard 8-bit vectorial nonlinearity over all 255 nonzero output masks and correctly evaluates it on stored S-boxes.
  - Source evidence: `sim.py`; `deposit_sim_regressions.py`; N4 verification report and 150 seed JSONs.
  - Allowed wording: “We optimized the minimum component nonlinearity over all nonzero 8-bit output masks.”
  - Forbidden stronger wording: “The generated S-box is secure” or “has AES-level strength.”
  - Uncertainty: The search is over a constrained map-generated permutation family and does not evaluate a complete cipher.
  - Next check: Independent final code review and citation/definition verification.
  - Decision: keep.
- Claim: FPA improves the search.
  - Source evidence: paired 30-seed contrasts; FPA LE is lower, and all three Holm-adjusted primary p-values are at least 0.057.
  - Allowed wording: “The paired experiment did not establish an FPA benefit; the LE contrast was negative before/after adjustment directionally, with Holm-adjusted p=0.057.”
  - Forbidden stronger wording: “FPA improves,” “equivalent,” or “no effect.”
  - Uncertainty: One seed schedule and constrained generator family; intervals retain uncertainty.
  - Next check: More seeds only if precision targets justify them; an expanded study is separate.
  - Decision: weaken.
- Claim: The 2D map outperforms the 1D map.
  - Source evidence: equal-call 30-seed 1D and 2D-random search outputs; descriptive summaries only.
  - Allowed wording: “The exploratory comparison did not show a clear advantage attributable to coupling.”
  - Forbidden stronger wording: “2D coupling is beneficial.”
  - Uncertainty: Distinct parameter-space dimensionalities; not paired as identical candidate sets.
  - Next check: A controlled ablation using shared parameter/state mapping if the mechanism is central.
  - Decision: weaken.

## Reproducibility

Run summary and reconciliation are generated by `analysis_revision.py` from immutable `run-v3/results/*.json`. Use `python output/sbox-revision/2026-10-19-inventory-01/analysis/analysis_revision.py` to regenerate. Results use fixed bootstrap seed 20261019 and 100,000 paired bootstrap resamples. Exact sign-flip uses a meet-in-the-middle enumeration over all 2^30 sign assignments; Holm correction covers the three primary endpoints within the FPA contrast.
