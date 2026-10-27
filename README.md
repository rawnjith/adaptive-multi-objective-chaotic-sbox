# Reproducibility materials: adaptive multi-objective chaotic S-box search

Code, data and analysis for a 30-seed paired simulation study of an 8-bit S-box generator driven by a
modulo-one coupled Logistic map and searched by three-objective NSGA-II with and without periodic Flower
Pollination Algorithm perturbations.

This repository contains **only** the reproducibility materials: simulator, frozen protocol, raw seed-level
outputs, analysis code, derived summaries and figures. It does not contain the manuscript, cover letter or
response documents.

## Scope of the evidence

- Objectives: largest Lyapunov exponent (maximised), standard vectorial S-box nonlinearity over all 255
  nonzero output masks (maximised), differential uniformity (minimised).
- Methods: NSGA-II, NSGA-II+FPA, equal-budget uniform parameter random search, uncoupled 1D Logistic-map
  search, and a descriptive uniform-random-permutation reference.
- Design: 30 fixed seeds x 5 methods = 150 cells, exactly 940 objective evaluations per cell.
- Everything is simulated: no cipher, device timing, energy or side-channel measurement is included.

## Layout

| Path | Contents |
|---|---|
| `sim.py` | current simulator/runner (CLI + objective and metric implementation) |
| `sim_run_v3.py` | recovered snapshot of the source that produced the archived runs |
| `protocol.md` | frozen metric, budget, seed list, comparator and statistical protocol |
| `run_grid.json` | per-cell record: method, seed, hashes, evaluation count, status |
| `results/` | 150 raw result JSONs (one per method x seed) |
| `logs/`, `checkpoints/` | per-cell run logs and provenance checkpoints |
| `analysis/` | verifier, analysis and figure scripts, derived summaries, QA records |
| `figures/`, `supplementary/` | vector exports and the three supplementary figures |
| `tests/` | regression tests for the simulator and analysis code |
| `analysis-stats.md` | human-readable statistics summary |

## Reproduce

Python 3.12.3, NumPy 1.26.4, Matplotlib 3.10.8, SciPy 1.14.1.

```bash
python analysis/verify_run_grid.py      # checks protocol/result hashes, 150 cells, 940 evaluations each
python analysis/analysis_revision.py    # regenerates analysis/results_summary.json and reconciliation.json
python analysis/make_figures.py         # regenerates figures/figure-01..03
python -m pytest tests -q
```

## Provenance note

Run artifacts record the run-time simulator hash
`04775ce3631504457859ed6a8ec0a79e837406785899b4ee62239f31469683cf`. The recovered snapshot
(`sim_run_v3.py`) and the current runner (`sim.py`) are later revisions and are **not** byte-identical to
that hash; integrity of the archived outputs rests on the recorded result hashes, the per-cell checkpoints
and the verifier, which pass for all 150 cells.
