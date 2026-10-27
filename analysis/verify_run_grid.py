from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
import sim

RUN_DIR = Path(__file__).resolve().parents[1]
EXPECTED_PROTOCOL = "662824ab0364f19f9f7e61ca7ddd69ad568780dbd38e5305786b4625a4e8d551"
EXPECTED_SIMULATOR = "04775ce3631504457859ed6a8ec0a79e837406785899b4ee62239f31469683cf"  # run-time code hash; current file hash is reported separately
EXPECTED_SEEDS = [104729, 130363, 155921, 196613, 228017, 262147, 300007, 350377, 400009, 450001,
                  500009, 550007, 600011, 650011, 700001, 750019, 800011, 850009, 900001, 950009,
                  1000003, 1050001, 1100009, 1150003, 1200007, 1250009, 1300003, 1350001, 1400009, 1450003]
METHODS = ["nsga_only", "fpa", "random_search", "1d_random_search", "random_permutation_reference"]
METRIC_ID = "vectorial_nonlinearity_min_over_255_nonzero_output_masks"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify() -> dict:
    assert sha256(RUN_DIR / "protocol.md") == EXPECTED_PROTOCOL, "protocol hash mismatch"
    grid = json.loads((RUN_DIR / "run_grid.json").read_text(encoding="utf-8"))
    current_simulator = sha256(REPO_ROOT / "sim.py")
    run_hashes = {node.get("simulator_hash") for node in grid}
    assert run_hashes == {EXPECTED_SIMULATOR}, f"run grid simulator hashes disagree: {run_hashes}"
    expected = {(method, seed) for method in METHODS for seed in EXPECTED_SEEDS}
    actual = {(node["method"], node["seed"]) for node in grid}
    assert len(grid) == 150 and actual == expected, "run-grid cells do not match protocol"
    assert all(node["status"] == "complete" for node in grid), "grid contains incomplete nodes"

    results: dict[tuple[str, int], dict] = {}
    paired_initial = {"nsga_only": {}, "fpa": {}}
    for node in grid:
        method, seed = node["method"], node["seed"]
        result_path = RUN_DIR / node["output"]
        log_path = RUN_DIR / node["log"]
        checkpoint_path = RUN_DIR / node["checkpoint"]
        assert result_path.is_file() and log_path.is_file() and checkpoint_path.is_file(), f"missing artifact: {method}/{seed}"
        assert sha256(result_path) == node["outputs_hash"]["result"], f"grid output hash mismatch: {method}/{seed}"
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        assert checkpoint["status"] == "complete" and checkpoint["exit_status"] == 0, f"failed checkpoint: {method}/{seed}"
        assert checkpoint["evaluation_count"] == 940, f"checkpoint call count: {method}/{seed}"
        assert checkpoint["input_sha256"]["protocol"] == EXPECTED_PROTOCOL, f"checkpoint protocol: {method}/{seed}"
        assert checkpoint["input_sha256"]["simulator"] == EXPECTED_SIMULATOR, f"checkpoint source hash: {method}/{seed}"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        assert result["evaluation_count"] == 940, f"result call count: {method}/{seed}"
        assert result["run_protocol_hash"] == EXPECTED_PROTOCOL, f"result protocol: {method}/{seed}"
        assert result["simulator_hash"] == EXPECTED_SIMULATOR, f"result source hash: {method}/{seed}"
        assert result["seed"] == seed and result["variant"] == method, f"result key mismatch: {method}/{seed}"
        assert result["metric_definition"] == METRIC_ID, f"metric definition: {method}/{seed}"
        results[(method, seed)] = result

        if method == "random_permutation_reference":
            assert result["permutation_count"] == 940 and len(result["metrics"]) == 940, f"permutation count: {seed}"
            for i, metric in enumerate(result["metrics"]):
                sbox = metric["sbox"]
                assert sorted(sbox) == list(range(256)), f"non-bijection random permutation: seed={seed} sample={i}"
                assert math.isfinite(metric["nl_min"]) and math.isfinite(metric["du"]) and math.isfinite(metric["lap"])
            continue

        sbox = result["sbox"]
        assert sorted(sbox) == list(range(256)), f"non-bijection selected S-box: {method}/{seed}"
        selected = result["selected"]
        assert all(math.isfinite(float(selected[k])) for k in ("le", "nl_min", "du", "lap", "sac", "bic_sac"))
        assert 0 <= selected["nl_min"] <= 128 and 2 <= selected["du"] <= 256
        computed_nl = sim.vectorial_nonlinearity(sbox)
        assert computed_nl == selected["nl_min"], f"selected NL mismatch: {method}/{seed}"
        recomputed_hv = sim.hypervolume_3d([item["normalized_objectives"] for item in result["pareto"]], [1.1, 1.1, 1.1])
        assert math.isclose(recomputed_hv, result["hypervolume"], rel_tol=1e-12, abs_tol=1e-12), f"HV mismatch: {method}/{seed}"
        if method in paired_initial:
            paired_initial[method][seed] = [(item["phase_x"], item["phase_y"], item["mu"], item["gamma"])
                                            for item in result["initial_population"]]

    for seed in EXPECTED_SEEDS:
        assert paired_initial["nsga_only"][seed] == paired_initial["fpa"][seed], f"paired initial populations differ: {seed}"
    return {"nodes": len(grid), "methods": len(METHODS), "seeds": len(EXPECTED_SEEDS),
            "completed": sum(node["status"] == "complete" for node in grid),
            "evaluation_budget": 940, "protocol_hash": EXPECTED_PROTOCOL,
            "run_time_simulator_hash": EXPECTED_SIMULATOR,
            "current_simulator_hash": current_simulator,
            "provenance_note": "run artifacts record the immutable run-time source hash; current source hash is a later analysis/runner-hardening revision and is not asserted equal",
            "status": "PASS"}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
