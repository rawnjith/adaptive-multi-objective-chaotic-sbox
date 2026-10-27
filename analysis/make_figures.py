from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
SUMMARY = json.loads((HERE / "results_summary.json").read_text(encoding="utf-8"))
METHODS = ["nsga_only", "fpa", "random_search", "1d_random_search"]
LABELS = ["NSGA-II", "NSGA-II + FPA", "2D random search", "1D random search"]
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]

plt.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def interval(ax, x, cell, mean_marker="D"):
    mean = cell["mean"]
    low, high = cell["bootstrap_ci_95"]
    ax.errorbar(x, mean, yerr=[[mean - low], [high - mean]], fmt=mean_marker,
                color="black", markersize=3.5, capsize=2, linewidth=1, zorder=5)


def figure_objectives(out: Path) -> None:
    rng = np.random.default_rng(20261019)
    metrics = [
        ("le", "Largest Lyapunov exponent", "higher is optimizer objective"),
        ("nl_min", "Vectorial NL$_{min}$", "higher is better"),
        ("du", "Differential uniformity", "lower is better"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.8), constrained_layout=True)
    for axis, (metric, title, direction) in zip(axes, metrics):
        values = [SUMMARY["per_seed"][method][metric] for method in METHODS]
        positions = np.arange(len(METHODS))
        boxes = axis.boxplot(values, positions=positions, widths=.56, showfliers=False,
                             patch_artist=True, medianprops={"color": "black", "linewidth": 1.1},
                             whiskerprops={"linewidth": .8}, capprops={"linewidth": .8})
        for box, color in zip(boxes["boxes"], COLORS):
            box.set_facecolor(color)
            box.set_alpha(.22)
        for index, (seed_values, color) in enumerate(zip(values, COLORS)):
            jitter = rng.uniform(-.11, .11, len(seed_values))
            axis.scatter(index + jitter, seed_values, s=12, color=color, alpha=.75,
                         edgecolors="white", linewidths=.25, zorder=3)
            interval(axis, index, SUMMARY["summary"][METHODS[index]][metric])
        axis.set_title(title)
        axis.set_ylabel(metric)
        axis.set_xlabel(direction, fontsize=7)
        axis.set_xticks(positions, LABELS, rotation=25, ha="right")
        axis.grid(axis="y", alpha=.22, linewidth=.5)
    fig.savefig(out / "figure-01-objective-comparison.pdf", bbox_inches="tight")
    fig.savefig(out / "figure-01-objective-comparison.svg", bbox_inches="tight")
    plt.close(fig)


def figure_paired_hv(out: Path) -> None:
    rng = np.random.default_rng(20261020)
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), constrained_layout=True)
    axis = axes[0]
    for index, metric in enumerate(("le", "nl_min", "du")):
        values = SUMMARY["comparisons"]["fpa_vs_nsga_only"][metric]["paired_differences"]
        jitter = rng.uniform(-.12, .12, len(values))
        axis.scatter(values, index + jitter, s=13, color=COLORS[1], alpha=.75,
                     edgecolors="white", linewidths=.25)
        contrast = SUMMARY["comparisons"]["fpa_vs_nsga_only"][metric]
        mean = contrast["mean_difference"]
        low, high = contrast["bootstrap_ci_95"]
        axis.errorbar(mean, index, xerr=[[mean - low], [high - mean]], fmt="D",
                      color="black", markersize=3.5, capsize=2, linewidth=1)
    axis.axvline(0, color="#555555", linewidth=.8, linestyle="--")
    axis.set_yticks(range(3), ["LE", "Vectorial NL$_{min}$", "DU"])
    axis.set_xlabel("FPA − NSGA-II (paired seed)")
    axis.set_title("Paired selected-objective differences")
    axis.grid(axis="x", alpha=.22, linewidth=.5)

    axis = axes[1]
    for index, (method, label, color) in enumerate(zip(METHODS, LABELS, COLORS)):
        values = SUMMARY["per_seed"][method]["hypervolume"]
        axis.scatter(index + rng.uniform(-.11, .11, len(values)), values, s=12,
                     color=color, alpha=.75, edgecolors="white", linewidths=.25)
        interval(axis, index, SUMMARY["summary"][method]["hypervolume"])
    axis.set_xticks(range(4), LABELS, rotation=25, ha="right")
    axis.set_ylabel("3D hypervolume (fixed protocol transform/reference)")
    axis.set_title("Per-seed Pareto-front quality")
    axis.grid(axis="y", alpha=.22, linewidth=.5)
    fig.savefig(out / "figure-02-paired-hypervolume.pdf", bbox_inches="tight")
    fig.savefig(out / "figure-02-paired-hypervolume.svg", bbox_inches="tight")
    plt.close(fig)


def figure_random_permutation(out: Path) -> None:
    grid = json.loads((RUN / "run_grid.json").read_text(encoding="utf-8"))
    by_seed = {}
    for row in grid:
        if row["method"] != "random_permutation_reference":
            continue
        record = json.loads((RUN / row["output"]).read_text(encoding="utf-8"))
        by_seed[row["seed"]] = record["metrics"]
    seeds = SUMMARY["seeds"]
    rng = np.random.default_rng(20261021)
    metrics = [("nl_min", "Vectorial NL$_{min}$"), ("du", "Differential uniformity"),
               ("lap", "LAP probability bias")]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.8), constrained_layout=True)
    colors = COLORS + ["#777777"]
    names = LABELS + ["Random permutations\n(seed means)"]
    for axis, (metric, title) in zip(axes, metrics):
        data = [SUMMARY["per_seed"][method][metric] for method in METHODS]
        # Random permutations are nested within seed; collapse to one mean per seed.
        data.append([float(np.mean([item[metric] for item in by_seed[seed]])) for seed in seeds])
        boxes = axis.boxplot(data, positions=np.arange(5), widths=.54, showfliers=False,
                             patch_artist=True, medianprops={"color": "black", "linewidth": 1.1},
                             whiskerprops={"linewidth": .8}, capprops={"linewidth": .8})
        for box, color in zip(boxes["boxes"], colors):
            box.set_facecolor(color)
            box.set_alpha(.22)
        for index, (seed_values, color) in enumerate(zip(data, colors)):
            axis.scatter(index + rng.uniform(-.11, .11, len(seed_values)), seed_values,
                         s=11, color=color, alpha=.72, edgecolors="white", linewidths=.2)
        axis.set_title(title)
        axis.set_xticks(np.arange(5), names, rotation=25, ha="right")
        axis.grid(axis="y", alpha=.22, linewidth=.5)
    fig.savefig(out / "figure-03-random-permutation-reference.pdf", bbox_inches="tight")
    fig.savefig(out / "figure-03-random-permutation-reference.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    out = HERE.parent / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figure_objectives(out)
    figure_paired_hv(out)
    figure_random_permutation(out)
    print("Rendered 3 figure families:", ", ".join(path.name for path in sorted(out.glob("*.pdf"))))


if __name__ == "__main__":
    main()
