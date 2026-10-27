import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR))
sys.path.insert(0, str(ROOT))

import sim

from analysis_revision import analyze_runs, exact_sign_flip_pvalue, holm_adjust, mean_sd, percentile_bootstrap_ci


def _selected_run(seed, method, le, nl, du, hv):
    initial_selected = {"le": le - .1, "nl_min": nl - 2, "du": du + 1}
    initial_population = [
        {"phase_x": index / 80, "phase_y": (index + 1) / 80, "mu": 3.6 + index / 1000, "gamma": index / 200}
        for index in range(80)
    ]
    record = {
        "seed": seed,
        "variant": method,
        "metric_definition": "vectorial_nonlinearity_min_over_255_nonzero_output_masks",
        "evaluation_count": 940,
        "selected": {"le": le, "nl_min": sim.vectorial_nonlinearity(sim.aes_sbox()), "du": du,
                     "sac": .5, "bic_sac": .03, "lap": .1, "gamma": .2},
        "objective_values": [-le, -nl, du],
        "normalized_objectives": [.4, .3, .2],
        "hypervolume": sim.hypervolume_3d([[.4, .3, .2]], [1.1, 1.1, 1.1]),
        "pareto_size": 4,
        "initial_selected": initial_selected,
        "initial_population": initial_population,
        "initial_pareto_size": 80,
        "sbox": list(range(256)),
        "pareto": [{"normalized_objectives": [.4, .3, .2]}],
    }
    return record


def test_mean_sd_and_bootstrap_ci_are_deterministic():
    mean, sd = mean_sd([1, 2, 3])
    assert mean == pytest.approx(2)
    assert sd == pytest.approx(1)
    first = percentile_bootstrap_ci([1, 2, 3, 4], resamples=5000, seed=7)
    second = percentile_bootstrap_ci([1, 2, 3, 4], resamples=5000, seed=7)
    assert first == second
    assert first[0] <= 2.5 <= first[1]


def test_exact_sign_flip_and_holm_match_small_enumeration():
    assert exact_sign_flip_pvalue([1, 2, 3]) == pytest.approx(.25)
    assert holm_adjust([.01, .04, .03]) == pytest.approx([.03, .06, .06])


def test_analyze_runs_keeps_seed_as_unit_and_reports_paired_contrasts():
    runs = []
    for seed, baseline, treatment in [(11, .4, .5), (12, .5, .45), (13, .6, .7), (14, .7, .65)]:
        for method, le, nl, du, hv in [
            ("nsga_only", baseline, 90 + seed % 2, 10, .3),
            ("fpa", treatment, 92, 10, .35),
            ("random_search", .35, 88, 12, .2),
            ("1d_random_search", .3, 86, 14, .15),
        ]:
            record = _selected_run(seed, method, le, nl, du, hv)
            sbox = sim.aes_sbox().tolist()
            record["sbox"] = sbox
            record["selected"]["nl_min"] = sim.vectorial_nonlinearity(sbox)
            runs.append(record)
    summary = analyze_runs(runs, bootstrap_resamples=1000, bootstrap_seed=19)
    assert summary["n_seeds"] == 4
    assert summary["per_seed"]["nsga_only"]["le"] == [.4, .5, .6, .7]
    contrast = summary["comparisons"]["fpa_vs_nsga_only"]["le"]
    assert contrast["n"] == 4
    assert contrast["mean_difference"] == pytest.approx(.025)
    assert "holm_adjusted_p" in contrast
    random_contrast = summary["comparisons"]["random_search_vs_nsga_only"]["le"]
    assert random_contrast["kind"] == "descriptive_only"
    assert len(random_contrast["bootstrap_ci_95"]) == 2
    assert random_contrast["exact_sign_flip_p"] is None
    assert len(contrast["bootstrap_ci_95"]) == 2
    assert summary["summary"]["nsga_only"]["hypervolume"]["mean"] == pytest.approx(sim.hypervolume_3d([[.4, .3, .2]], [1.1, 1.1, 1.1]))
    assert "initial_le" in summary["summary"]["nsga_only"]
    assert "initial_le" not in summary["summary"]["random_search"]
    assert summary["claim_scope"] == "simulation-only; no whole-cipher or hardware claim"


def test_analysis_rejects_pareto_hypervolume_mismatch():
    records = []
    for method in ("nsga_only", "fpa", "random_search", "1d_random_search"):
        record = _selected_run(1, method, .5, 90, 10, .6)
        record["sbox"] = sim.aes_sbox().tolist()
        record["selected"]["nl_min"] = sim.vectorial_nonlinearity(record["sbox"])
        record["pareto"][0]["normalized_objectives"] = [.2, .3, .4]
        record["hypervolume"] = .01
        records.append(record)
    with pytest.raises(ValueError, match="hypervolume mismatch"):
        analyze_runs(records, bootstrap_resamples=100, bootstrap_seed=1)


def test_analysis_rejects_unpaired_initial_populations():
    records = []
    for method in ("nsga_only", "fpa", "random_search", "1d_random_search"):
        record = _selected_run(1, method, .5, 90, 10, .6)
        record["sbox"] = sim.aes_sbox().tolist()
        record["selected"]["nl_min"] = sim.vectorial_nonlinearity(record["sbox"])
        record["initial_population"] = [{"phase_x": index / 80, "phase_y": (index + 1) / 80, "mu": 3.6 + index / 1000, "gamma": index / 200} for index in range(80)]
        record["initial_pareto_size"] = 80
        records.append(record)
    records[1]["initial_population"][0]["phase_x"] = .9
    with pytest.raises(ValueError, match="initial populations differ"):
        analyze_runs(records, bootstrap_resamples=100, bootstrap_seed=1)


def test_analysis_rejects_missing_initial_population_diagnostics():
    record = _selected_run(1, "nsga_only", .4, 90, 10, .3)
    record["sbox"] = sim.aes_sbox().tolist()
    record["selected"]["nl_min"] = sim.vectorial_nonlinearity(record["sbox"])
    record["initial_selected"] = None
    record["initial_population"] = None
    full = [record]
    for method in ("fpa", "random_search", "1d_random_search"):
        other = _selected_run(1, method, .4, 90, 10, .3)
        other["sbox"] = sim.aes_sbox().tolist()
        other["selected"]["nl_min"] = sim.vectorial_nonlinearity(other["sbox"])
        full.append(other)
    with pytest.raises(ValueError, match="initial-population diagnostics missing"):
        analyze_runs(full, bootstrap_resamples=100, bootstrap_seed=1)


def test_smoke_analysis_with_single_seed_marks_uncertainty_unavailable():
    records = []
    for method in ("nsga_only", "fpa", "random_search", "1d_random_search"):
        record = _selected_run(1, method, .5, 90, 10, .6)
        record["sbox"] = sim.aes_sbox().tolist()
        record["selected"]["nl_min"] = sim.vectorial_nonlinearity(record["sbox"])
        records.append(record)
    summary = analyze_runs(records, bootstrap_resamples=50, bootstrap_seed=1)
    assert summary["n_seeds"] == 1
    assert summary["summary"]["nsga_only"]["le"]["sample_sd"] is None
    assert summary["summary"]["nsga_only"]["le"]["bootstrap_ci_95"] is None


def test_analysis_rejects_missing_paired_seed_or_bad_metric_provenance():
    runs = [_selected_run(1, "nsga_only", .4, 90, 10, .3),
            _selected_run(2, "fpa", .5, 92, 10, .35)]
    with pytest.raises(ValueError, match="paired seed"):
        analyze_runs(runs, bootstrap_resamples=100, bootstrap_seed=1)
    bad = _selected_run(1, "nsga_only", .4, 90, 10, .3)
    bad["metric_definition"] = "coordinate_only"
    with pytest.raises(ValueError, match="metric definition"):
        analyze_runs([bad], bootstrap_resamples=100, bootstrap_seed=1)
