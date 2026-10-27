from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
import sim

PRIMARY = ("le", "nl_min", "du")
SECONDARY = ("sac", "bic_sac", "lap")
SEARCH_METHODS = ("nsga_only", "fpa", "random_search", "1d_random_search")
METRIC_ID = "vectorial_nonlinearity_min_over_255_nonzero_output_masks"


def mean_sd(values: Sequence[float]) -> tuple[float, float]:
    data = np.asarray(values, dtype=float)
    if data.ndim != 1 or data.size < 2 or not np.all(np.isfinite(data)):
        raise ValueError("mean_sd requires at least two finite observations")
    return float(data.mean()), float(data.std(ddof=1))


def percentile_bootstrap_ci(values: Sequence[float], resamples: int = 100_000, seed: int = 20261019) -> list[float]:
    data = np.asarray(values, dtype=float)
    if data.ndim != 1 or data.size < 2 or not np.all(np.isfinite(data)):
        raise ValueError("bootstrap CI requires at least two finite observations")
    if resamples < 1:
        raise ValueError("resamples must be positive")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, data.size, size=(resamples, data.size))
    means = data[indices].mean(axis=1)
    return [float(x) for x in np.quantile(means, [.025, .975])]


def exact_sign_flip_pvalue(differences: Sequence[float], tolerance: float = 1e-12) -> float:
    data = np.asarray(differences, dtype=float)
    if data.ndim != 1 or data.size == 0 or not np.all(np.isfinite(data)):
        raise ValueError("differences must be a non-empty finite vector")
    split = data.size // 2

    def signed_sums(part: np.ndarray) -> np.ndarray:
        sums = np.asarray([0.0])
        for value in part:
            sums = np.concatenate((sums + value, sums - value))
        return sums

    left = signed_sums(data[:split])
    right = np.sort(signed_sums(data[split:]))
    observed = abs(float(data.sum()))
    threshold = observed - tolerance
    if threshold <= 0:
        return 1.0
    count = 0
    for subtotal in left:
        count += int(np.searchsorted(right, -threshold - subtotal, side="right"))
        count += int(right.size - np.searchsorted(right, threshold - subtotal, side="left"))
    return float(count / (2 ** data.size))


def holm_adjust(p_values: Sequence[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must be a finite vector in [0,1]")
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    n = len(p)
    for rank, index in enumerate(order):
        running = max(running, (n - rank) * p[index])
        adjusted[index] = min(1.0, running)
    return [float(x) for x in adjusted]


def exact_sign_flip_distribution(differences: Sequence[float], tolerance: float = 1e-12) -> dict[str, Any]:
    data = np.asarray(differences, dtype=float)
    if data.ndim != 1 or data.size == 0 or not np.all(np.isfinite(data)):
        raise ValueError("differences must be a non-empty finite vector")
    left = [0.0]
    for value in data[:data.size // 2]:
        left = [subtotal + sign * float(value) for subtotal in left for sign in (-1, 1)]
    right = np.sort(np.asarray([0.0]))
    for value in data[data.size // 2:]:
        right = np.concatenate((right + value, right - value))
        right.sort()
    observed = abs(float(data.sum()))
    threshold = observed - tolerance
    if threshold <= 0:
        extreme = 2 ** data.size
    else:
        extreme = 0
        for subtotal in left:
            low = np.searchsorted(right, -threshold - subtotal, side="right")
            high = right.size - np.searchsorted(right, threshold - subtotal, side="left")
            extreme += int(low + high)
    return {"observed_abs_sum": observed, "extreme_assignments": int(extreme),
            "total_assignments": 2 ** data.size, "p_value": float(extreme / (2 ** data.size)),
            "tolerance": tolerance, "algorithm": "meet-in-the-middle exact enumeration"}


def cohen_dz(differences: Sequence[float]) -> float:
    data = np.asarray(differences, dtype=float)
    if data.ndim != 1 or data.size < 2 or not np.all(np.isfinite(data)):
        raise ValueError("Cohen dz requires at least two finite paired differences")
    sd = float(data.std(ddof=1))
    if sd == 0:
        return 0.0 if float(data.mean()) == 0 else math.copysign(math.inf, float(data.mean()))
    return float(data.mean() / sd)


def _paired_effect(differences: np.ndarray, resamples: int, seed: int) -> dict[str, Any]:
    if differences.size == 1:
        return {"kind": "descriptive_only", "n": 1, "mean_difference": float(differences[0]),
                "sample_sd_difference": None, "bootstrap_ci_95": None, "exact_sign_flip_p": None,
                "holm_adjusted_p": None, "cohen_dz": None,
                "paired_differences": [float(differences[0])]}
    mean, sd = mean_sd(differences)
    return {
        "n": int(differences.size),
        "mean_difference": mean,
        "sample_sd_difference": sd,
        "bootstrap_ci_95": percentile_bootstrap_ci(differences, resamples, seed),
        **exact_sign_flip_distribution(differences),
        "exact_sign_flip_p": exact_sign_flip_pvalue(differences),
        "cohen_dz": cohen_dz(differences),
        "paired_differences": [float(x) for x in differences],
    }


def _load_records(run_files: Iterable[dict[str, Any] | str | Path]) -> list[dict[str, Any]]:
    records = []
    for item in run_files:
        data = json.loads(Path(item).read_text(encoding="utf-8")) if isinstance(item, (str, Path)) else item
        records.append(data)
    return records


def _summary(values: Sequence[float], resamples: int, seed: int) -> dict[str, Any]:
    if len(values) == 1:
        return {"mean": float(values[0]), "sample_sd": None, "n": 1, "bootstrap_ci_95": None}
    mean, sd = mean_sd(values)
    return {"mean": mean, "sample_sd": sd, "n": len(values),
            "bootstrap_ci_95": percentile_bootstrap_ci(values, resamples, seed)}


def analyze_runs(run_files: Iterable[dict[str, Any] | str | Path], protocol: dict[str, Any] | None = None,
                 bootstrap_resamples: int = 100_000, bootstrap_seed: int = 20261019) -> dict[str, Any]:
    records = _load_records(run_files)
    if not records:
        raise ValueError("no run records supplied")
    by_method: dict[str, dict[int, dict[str, Any]]] = {}
    permutation_samples: dict[int, list[dict[str, Any]]] = {}
    for record in records:
        method = record.get("variant")
        if record.get("metric_definition") != METRIC_ID:
            raise ValueError("metric definition mismatch")
        if method == "random_permutation_reference":
            seed = int(record["seed"])
            permutation_samples[seed] = record["metrics"]
            continue
        if method not in SEARCH_METHODS:
            raise ValueError(f"unknown analysis method: {method}")
        seed = int(record["seed"])
        selected = record["selected"]
        if "sbox" not in record:
            raise ValueError(f"selected S-box missing for paired seed {seed}")
        computed = sim.vectorial_nonlinearity(record["sbox"])
        if computed != selected["nl_min"]:
            raise ValueError(f"selected vectorial NL mismatch for paired seed {seed}: stored {selected['nl_min']} recomputed {computed}")
        if int(record.get("evaluation_count", -1)) != 940:
            raise ValueError(f"evaluation budget mismatch for paired seed {seed}")
        recomputed_hv = sim.hypervolume_3d([item["normalized_objectives"] for item in record["pareto"]], [1.1, 1.1, 1.1])
        if not math.isclose(recomputed_hv, float(record["hypervolume"]), rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"hypervolume mismatch for {method}/{seed}")
        if method in ("nsga_only", "fpa"):
            initial_population = record.get("initial_population") or []
            if int(record.get("initial_pareto_size", -1)) < 1 or len(initial_population) != 80:
                raise ValueError(f"initial-population diagnostics missing or wrong size for {method}/{seed}")
            if any(not all(math.isfinite(float(item[key])) for key in ("phase_x", "phase_y", "mu", "gamma"))
                   for item in initial_population):
                raise ValueError(f"invalid initial population chromosome for {method}/{seed}")
        if int(record.get("evaluation_count", -1)) != 940:
            raise ValueError(f"evaluation budget mismatch for paired seed {seed}")
        if seed in by_method.setdefault(method, {}):
            raise ValueError(f"duplicate method/paired seed: {method}/{seed}")
        by_method[method][seed] = record

    required_methods = set(SEARCH_METHODS)
    if set(by_method) != required_methods:
        raise ValueError(f"incomplete methods: expected {sorted(required_methods)}, got {sorted(by_method)}")
    seed_sets = {method: set(items) for method, items in by_method.items()}
    common_seeds = set.intersection(*seed_sets.values())
    if not common_seeds or any(seeds != common_seeds for seeds in seed_sets.values()):
        raise ValueError("methods must have the same paired seed set")
    ordered_seeds = sorted(common_seeds)
    for seed in ordered_seeds:
        nsga_initial = by_method["nsga_only"][seed]["initial_population"] or []
        fpa_initial = by_method["fpa"][seed]["initial_population"] or []
        nsga_genes = [(item["phase_x"], item["phase_y"], item["mu"], item["gamma"]) for item in nsga_initial]
        fpa_genes = [(item["phase_x"], item["phase_y"], item["mu"], item["gamma"]) for item in fpa_initial]
        if nsga_genes != fpa_genes:
            raise ValueError(f"paired initial populations differ for seed {seed}")

    per_seed: dict[str, dict[str, list[float]]] = {}
    summaries: dict[str, dict[str, Any]] = {}
    for method in SEARCH_METHODS:
        metrics = (*PRIMARY, *SECONDARY, "hypervolume", "gamma", "coordinate_nl_min", "front_size", "evaluation_count")
        if method in ("nsga_only", "fpa"):
            metrics += ("initial_le", "initial_nl_min", "initial_du", "initial_front_size")
        per_seed[method] = {metric: [] for metric in metrics}
        for seed in ordered_seeds:
            record = by_method[method][seed]
            selected = record["selected"]
            for metric in (*PRIMARY, *SECONDARY):
                per_seed[method][metric].append(float(selected[metric]))
            per_seed[method]["hypervolume"].append(float(record["hypervolume"]))
            per_seed[method]["gamma"].append(float(selected["gamma"]))
            per_seed[method]["coordinate_nl_min"].append(float(sim.coordinate_nonlinearity(record["sbox"])[0]))
            if method in ("nsga_only", "fpa"):
                initial = record["initial_selected"]
                if initial is None or record.get("initial_pareto_size") is None:
                    raise ValueError(f"initial-population diagnostics missing for {method}/{seed}")
                per_seed[method]["initial_le"].append(float(initial["le"]))
                per_seed[method]["initial_nl_min"].append(float(initial["nl_min"]))
                per_seed[method]["initial_du"].append(float(initial["du"]))
                per_seed[method]["initial_front_size"].append(float(record["initial_pareto_size"]))
            per_seed[method]["front_size"].append(float(record["pareto_size"]))
            per_seed[method]["evaluation_count"].append(float(record["evaluation_count"]))
        summaries[method] = {key: _summary(values, bootstrap_resamples, bootstrap_seed)
                             for key, values in per_seed[method].items()}

    contrasts: dict[str, dict[str, Any]] = {}
    for method in ("fpa", "random_search", "1d_random_search"):
        name = f"{method}_vs_nsga_only"
        contrasts[name] = {}
        if method != "fpa":
            for metric in PRIMARY:
                differences = np.asarray(per_seed[method][metric]) - np.asarray(per_seed["nsga_only"][metric])
                if differences.size == 1:
                    contrasts[name][metric] = {"kind": "descriptive_only", "n": 1,
                                               "mean_difference": float(differences[0]),
                                               "sample_sd_difference": None, "cohen_dz": None,
                                               "paired_differences": [float(differences[0])]}
                    continue
                mean, sd = mean_sd(differences)
                contrasts[name][metric] = {"kind": "descriptive_only", "n": int(differences.size),
                                           "mean_difference": mean, "sample_sd_difference": sd,
                                           "bootstrap_ci_95": percentile_bootstrap_ci(differences, bootstrap_resamples, bootstrap_seed + len(contrasts)),
                                           "exact_sign_flip_p": None, "holm_adjusted_p": None,
                                           "cohen_dz": cohen_dz(differences),
                                           "paired_differences": [float(x) for x in differences]}
            continue
        for index, metric in enumerate(PRIMARY):
            differences = np.asarray(per_seed[method][metric]) - np.asarray(per_seed["nsga_only"][metric])
            contrasts[name][metric] = _paired_effect(differences, bootstrap_resamples, bootstrap_seed + index)
        raw_p = [contrasts[name][metric]["exact_sign_flip_p"] for metric in PRIMARY]
        adjusted = holm_adjust(raw_p) if all(value is not None for value in raw_p) else [None] * len(raw_p)
        for metric, p_adj in zip(PRIMARY, adjusted):
            contrasts[name][metric]["holm_adjusted_p"] = p_adj

    permutation_summary = {}
    for seed, metrics in sorted(permutation_samples.items()):
        permutation_summary[str(seed)] = {key: [float(item[key]) for item in metrics]
                                          for key in ("nl_min", "du", "lap", "sac", "bic_sac")}
    return {
        "n_seeds": len(ordered_seeds),
        "seeds": ordered_seeds,
        "metric_definition": METRIC_ID,
        "primary_metrics": list(PRIMARY),
        "secondary_metrics": list(SECONDARY),
        "unit_of_analysis": "independent random seed; FPA/NSGA-II paired on shared initial population",
        "summary": summaries,
        "per_seed": per_seed,
        "comparisons": contrasts,
        "random_permutation_reference": permutation_summary,
        "claim_scope": "simulation-only; no whole-cipher or hardware claim",
        "bootstrap": {"method": "paired percentile bootstrap of seed-level mean", "resamples": bootstrap_resamples,
                      "seed": bootstrap_seed},
        "multiple_comparisons": "Holm adjustment over LE, vectorial NL_min, DU within each prespecified contrast",
        "protocol": protocol,
    }


def reconciliation_rows(summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for method, metrics in summary["summary"].items():
        for metric, cell in metrics.items():
            rows.append({"method": method, "metric": metric, **cell})
    for contrast, metrics in summary["comparisons"].items():
        for metric, cell in metrics.items():
            rows.append({"contrast": contrast, "metric": metric, **cell})
    return rows


def main() -> None:
    run_root = Path(__file__).resolve().parents[1]
    grid = json.loads((run_root / "run_grid.json").read_text(encoding="utf-8"))
    if any(node.get("status") != "complete" for node in grid):
        raise SystemExit("N4 grid is incomplete; refusing to aggregate")
    paths = [run_root / node["output"] for node in grid]
    summary = analyze_runs(paths)
    output_root = Path(__file__).resolve().parent
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "results_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (output_root / "reconciliation.json").write_text(json.dumps(reconciliation_rows(summary), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"n_seeds": summary["n_seeds"], "methods": list(summary["summary"]),
                      "contrasts": list(summary["comparisons"]), "outputs": [str(output_root / "results_summary.json"), str(output_root / "reconciliation.json")]}, indent=2))


if __name__ == "__main__":
    main()
