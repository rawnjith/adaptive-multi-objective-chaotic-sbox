import itertools
import json
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = Path(os.environ.get("SBOX_RESULTS_DIR", HERE / "results"))
SEEDS = [42, 100, 1234, 2025, 7777]
METRICS = ["le", "nl_min", "sac", "bic_sac", "du", "lap"]
BOOTSTRAP_SEED = 20261007
BOOTSTRAP_RESAMPLES = 100_000


def selected(method, seed, metric):
    path = RESULTS / f"depositfinal_{method}_seed{seed}.json"
    return json.loads(path.read_text(encoding="utf-8"))["selected"][metric]


def exact_sign_flip_p(values):
    observed = abs(float(np.mean(values)))
    null = [abs(float(np.mean(values * np.asarray(signs))))
            for signs in itertools.product([-1, 1], repeat=len(values))]
    return sum(value >= observed - 1e-12 for value in null) / len(null)


def holm_adjust(p_values):
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values), dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (len(p_values) - rank) * p_values[index]))
        adjusted[index] = running
    return adjusted


def main():
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    rows = []
    for metric in METRICS:
        fpa = np.asarray([selected("fpa", seed, metric) for seed in SEEDS], dtype=float)
        baseline = np.asarray([selected("nsga_only", seed, metric) for seed in SEEDS], dtype=float)
        delta = fpa - baseline
        bootstrap = np.asarray([
            np.mean(rng.choice(delta, size=len(delta), replace=True))
            for _ in range(BOOTSTRAP_RESAMPLES)
        ])
        rows.append({
            "metric": metric,
            "fpa_mean": float(fpa.mean()),
            "fpa_sd": float(fpa.std(ddof=1)),
            "nsga_only_mean": float(baseline.mean()),
            "nsga_only_sd": float(baseline.std(ddof=1)),
            "paired_delta": float(delta.mean()),
            "bootstrap_ci_95": [float(x) for x in np.quantile(bootstrap, [0.025, 0.975])],
            "exact_sign_flip_p": exact_sign_flip_p(delta),
            "paired_differences": [float(x) for x in delta],
        })
    adjusted = holm_adjust([row["exact_sign_flip_p"] for row in rows])
    for row, p_holm in zip(rows, adjusted):
        row["holm_adjusted_p"] = float(p_holm)
    fronts = {
        method: [json.loads((RESULTS / f"depositfinal_{method}_seed{s}.json").read_text(encoding="utf-8"))["pareto_size"] for s in SEEDS]
        for method in ("fpa", "nsga_only")
    }
    output = {
        "seeds": SEEDS,
        "n_pairs": len(SEEDS),
        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "test": "exact two-sided paired sign-flip; Holm adjustment across six metrics",
        "metrics": rows,
        "pareto_front_sizes": fronts,
    }
    out_path = HERE / "paired_statistics.json"
    out_path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
