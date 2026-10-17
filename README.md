# Reproducibility release — S-box optimizer study

## Contents
- `sim.py`: map-driven S-box generator, NSGA-II, and optional FPA perturbation.
- `deposit_sim_regressions.py`: regression tests for the corrected scientific implementation and CLI.
- `results/`: final paired JSON outputs for seeds 42, 100, 1234, 2025, 7777 (`depositfinal_fpa_seedN.json`, `depositfinal_nsga_only_seedN.json`).
- `figures/`: manuscript figure PDFs, including the paired ablation plot.
- `plot_verified_ablation.py`: regenerates the NLmin paired figure from `results/`.
- `analyze_paired_runs.py` and `paired_statistics.json`: paired bootstrap and exact sign-flip analysis.
- `analysis-stats.md`: human-readable summary of methods and limitations.

## Environment and configuration
Verified with Python 3.12.3, NumPy 1.26.4, Matplotlib 3.10.8, SciPy 1.14.1. Runs use population 80, 20 generations, per-gene crossover probability 0.5, Gaussian mutation SD 0.01. Search genes are phase variables p_x,p_y in [0,1), mapped to map states x0=sin²(πp_x), y0=sin²(πp_y); mu is in [3.57,4.0), gamma in [0,0.5). FPA is applied every three generations to at most ten members of the first nondominated front. Seeds are paired across variants: 42, 100, 1234, 2025, 7777.

## Reproduce runs
Run from this directory. For each seed, replace `SEED` with one of `42`, `100`, `1234`, `2025`, `7777`:

```sh
python sim.py --seed SEED --output results/depositfinal_fpa_seedSEED.json
python sim.py --seed SEED --no-fpa --output results/depositfinal_nsga_only_seedSEED.json
```

Omit `--figures` for run-only regeneration. Add `--figures` to generate the single-run figures under `figures/`.

## Reproduce analysis and paired figure

```sh
python analyze_paired_runs.py
python plot_verified_ablation.py
```

The scripts default to this bundle's `results/` and `figures/` directories and write `paired_statistics.json` and `figures/fig_08_paired_ablation.pdf`.

## Regression tests

```sh
python -m pytest deposit_sim_regressions.py -q
```

## Scope and limitations
This is a simulation study of an S-box permutation, not a complete cipher evaluation. No key schedule/diffusion construction was specified, so NPCR/UACI are omitted. No MCU runtime, energy, or side-channel measurements are included. Five seed pairs are exploratory and do not establish FPA superiority or equivalence. The selected solutions include near-zero coupling; comparison with a one-dimensional generator remains future work.

## Public archive status
The public GitHub repository is https://github.com/rawnjith/adaptive-multi-objective-chaotic-sbox. The manuscript PDF and TeX on the current `main` branch are the latest working revision. Zenodo v1.0.0 is archived at https://doi.org/10.5281/zenodo.23439389 and points to GitHub commit `05d95b523a39a1a3daf654be88998324bc721f42`; it does not contain the later manuscript/source/data updates present on current `main`. Do not cite DOI 10.5281/zenodo.23439389 as archiving the newer revision. After the final release contents are approved, create a new tagged GitHub release, archive it in Zenodo, then cite the new version DOI in the manuscript and keep the GitHub URL as the repository address.
