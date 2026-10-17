# Technical inspection response matrix — revised draft status

| Review item | Work completed in current revision | Evidence / remaining action |
|---|---|---|
| 1. Objective count mismatch | Fitness defined as 3 independent objectives (LE, NLmin, DU); SAC, BIC-SAC deviation, LAP post hoc | `main_Scientific_Reports.tex` Eq. `eq:fit`, Abstract/Methods aligned |
| 2. Runtime is a DU transform | Removed runtime objective and all Cortex-M timing claims | No physical device measurements in workspace |
| 3. LE near 16.8 implausible | Corrected tangent normalization; expected ln(2) check for uncoupled mu=4; results now around 0.67 | 14 regression tests pass; final seed outputs include LE |
| 4. Single-seed headline / FPA contribution | Five paired seeds, explicit normalized ideal-point compromise rule; FPA claim withdrawn | Exact paired tests, Holm adjustment, bootstrap CIs reported. n=5 is low-power; future larger seed set/hypervolume recommended |
| 5. Hardware runtime | Removed unsupported timing/energy claims; paper explicitly simulation-only | Actual Cortex-M benchmarking remains future experiment |
| 6. Code/data availability | Current source, manuscript, results, analyses and figures are on public GitHub `main`; manuscript Data Availability cites that repository | Zenodo v1.0.0 checksums differ from current simulator/results/figure. Do not cite its DOI for this revision; create a new tagged Zenodo version and insert its DOI before resubmission |
| 7. Metrics | Removed entropy as informative outcome; fixed LAP probability-bias scale; removed NPCR/UACI for lack of full cipher | No complete cipher/key schedule found; need author-provided design if retaining these metrics |
| 8. Baselines | Added like-for-like AES 8-bit comparison in prose using the same metrics; removed incompatible/missing static chaotic and GA baselines | AES generator tests: prefix, bijection, NL=112, DU=4, LAP=0.0625. Valid recent chaotic like-for-like baselines still need careful verification/reproduction |
| 9. Citations | Removed orphan items and suspicious TIFS DOI; corrected Adams–Tavares DOI; verified core records via Crossref/arXiv/OpenAlex | Remaining 9 cited refs metadata checked; full text claim verification is partial, especially some technical definitions. Do not claim complete full bibliography audit |
| 10. Statistics | Added five paired seed comparisons; bootstrap 95% CI; exact paired sign-flip tests; Holm correction | All corrected p_Holm=1.0; underpowered; no equivalence or superiority claim |
| 11. Contradictions | Updated generations 20, metrics/objectives, corrected LE values, headline means, SAC precision | Automated SciRep checks pass; confirm final tables/captions from current PDF after last compile |
| 12. gamma=0 | Reparameterized phase genes to valid state values; final selections are gamma 0.0004–0.4923, including near-zero | No 1D-vs-2D comparison yet; paper states coupling benefit not established |

## Citation traceability status
The nine citations retained in the revised text were checked. Metadata checked for the five DOI sources (Adams/Tavares, Carlet, Matsui, Nyberg, Biham/Shamir) and the Yang FPA DOI; full-text checked for Cassal-Quiroga arXiv:1902.06412, Kuznetsov arXiv:2407.03510, and McLaughlin/Clark arXiv:1301.6972. Note: source full text was not retrieved for the DOI chapters/articles; only their metadata and high-level topic were checked. The old Adams DOI was corrected; suspicious TIFS DOI was removed with the unused bibliography entries. This is not a full-text verification of every detailed claim in those sources.

## Upload blockers
1. Public repository DOI not yet available. Do not upload this draft to the journal until deposited and Data Availability statement updated with permanent DOI/link.
2. Citation review is bounded as above; source-specific detailed claims should be checked if the literature discussion expands.
3. Rasterized pages were produced, but visual image inspection is unavailable in this session; a human should review the rendered PDF layout.
4. The response letter and cover letter have not yet been reconciled with this substantive revision.
5. There is a title/page-header overfull box warning; it should be corrected or confirmed benign after visual review.

