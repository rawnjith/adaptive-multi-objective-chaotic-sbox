# N5 Visual QA

All three figure families were rendered to PDF+SVG and rasterized to PNG at 180 DPI. Each raster was visually inspected in this session.

## Figure 1 — objective-comparison
- Renders all three panels (LE, vectorial NLmin, DU) with per-seed jittered points, box summaries, and mean markers with bootstrap 95% CI bars.
- Labels legible; no clipping; x labels rotated for readability.
- DU panel is heavily concentrated at 10 with one 12 outlier per method; boxplot is honest about the discreteness.
- Observation: LE distributions overlap; FPA's mean LE is slightly lower than NSGA-II-only; vectorial NLmin distributions overlap heavily. No superiority read from this figure.

## Figure 2 — paired-hypervolume
- Left panel: paired FPA−NSGA-II differences with zero line; LE deltas are mostly near zero/negative; NLmin deltas are centered near 0 with one +6 outlier; DU delta centered at 0.
- Right panel: per-seed hypervolume by method with mean + bootstrap CI; all four methods cluster in a narrow range 0.710–0.731; intervals overlap.
- Observation: hypervolume separation is small; descriptive only.
- Minor visual note: left-panel y order is LE (bottom), NLmin (middle), DU (top), matching code order; acceptable but could be reordered to the paper's endpoint order.

## Figure 3 — random-permutation-reference
- Three panels (NLmin, DU, LAP) with four method distributions plus a per-seed-mean random-permutation reference.
- Random permutation seed-means are visibly lower in NLmin (~92–93), higher in DU (~11.3) and higher in LAP (~0.138) than optimizer selections.
- Observation: the random-permutation reference is clearly separated from search outputs on NLmin/DU/LAP; no inferential claim from nested reference draws.

## QA result
- Vector exports exist and are non-empty; raster previews inspected.
- No clipping/overlap defects observed at 180 DPI.
- Human final-size check at the target journal column width is still recommended before submission (noted in figure-catalog.md).
