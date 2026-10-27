"""NSGA-II-FPA / 2D-Logistic S-box simulator — honest, reproducible.

Seed = 20260814.
Outputs:
  figures/*.pdf           publication-grade figures
  results.json            full Pareto front + selected S-box + metrics
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import time
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


SEED = 20260814
np.random.seed(SEED)
# Allow runtime seed override via tools/sim.py --seed=N (handled in __main__)
FIG_DIR = Path("figures")
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.grid": False,
    "lines.linewidth": 1.2, "lines.markersize": 3.5,
    "savefig.dpi": 400, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
WONG = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7"]


# ---------------------------- 2D Logistic map ----------------------------
def step_2d_logistic(x, y, mu, gamma):
    """Apply one simultaneous modulo-one step of the coupled map."""
    x_next = mu * x * (1 - x) + gamma * y
    y_next = mu * y * (1 - y) + gamma * x
    return x_next - np.floor(x_next), y_next - np.floor(y_next)


def initial_state(x_gene, y_gene):
    """Map phase genes in [0,1) to valid unit-square Logistic states."""
    return float(np.sin(np.pi * x_gene) ** 2), float(np.sin(np.pi * y_gene) ** 2)


def wrap_parameters(chromosome):
    """Wrap each optimizer gene periodically within its search interval."""
    x_gene, y_gene, mu, gamma = chromosome
    wrapped = [
        float(x_gene % 1.0),
        float(y_gene % 1.0),
        3.57 + ((mu - 3.57) % 0.43),
        gamma % 0.5,
    ]
    return wrapped


def iterate_2d_logistic(x0, y0, mu, gamma, n, transient=200):
    x = np.empty(n)
    y = np.empty(n)
    x_state, y_state = x0, y0
    for i in range(n + transient):
        x_state, y_state = step_2d_logistic(x_state, y_state, mu, gamma)
        if i >= transient:
            x[i - transient], y[i - transient] = x_state, y_state
    return x, y


def lyapunov_1d_logistic(x0, mu, n=1500, transient=200):
    """Estimate the largest Lyapunov exponent for the 1D logistic recurrence."""
    x = float(x0)
    total = 0.0
    for index in range(n + transient):
        derivative = abs(float(mu) * (1.0 - 2.0 * x))
        if derivative <= np.finfo(float).tiny:
            return float("-inf")
        x = (mu * x * (1.0 - x)) % 1.0
        if index >= transient:
            total += np.log(derivative)
    return float(total / n)


def lyapunov_2d_logistic(x0, y0, mu, gamma, n=1500, transient=200):
    """Estimate the largest Lyapunov exponent using a unit tangent vector."""
    le_sum = 0.0
    x, y = x0, y0
    tangent = np.array([1.0, 0.0], dtype=float)
    for i in range(n + transient):
        x_next, y_next = step_2d_logistic(x, y, mu, gamma)
        jacobian = np.array([
            [mu * (1 - 2 * x), gamma],
            [gamma, mu * (1 - 2 * y)],
        ])
        tangent = jacobian @ tangent
        norm = np.linalg.norm(tangent)
        if not np.isfinite(norm) or norm <= np.finfo(float).tiny:
            return float("-inf")
        tangent /= norm
        if i >= transient:
            le_sum += np.log(norm)
        x, y = x_next, y_next
    return le_sum / n


# ---------------------------- map -> S-box -------------------------------
def map_to_sbox(x0, y0, mu, gamma, transient=500, target_len=256):
    """Construct an 8-bit S-box by using the 2D Logistic map as a Fisher-Yates
    shuffle driver (Jakimoski-Kocarev 2001, Ozkaynak 2010 variant).

    Procedure:
      1. Burn 500 transient iterates.
      2. Collect 256 chaotic key values k_0..k_{255}.
      3. For i = 255 down to 1:
            j = floor(k_i * (i + 1))   # uniform in [0, i]
            swap(perm[i], perm[j])
      This produces a deterministic permutation of length 256 with no collisions.
      The Fisher-Yates step uses chaotic values to mix — the resulting S-box
      retains strong avalanche because small input changes perturb many keys.
    """
    perm = np.arange(target_len, dtype=np.int64)
    x, y = x0, y0
    for _ in range(transient):
        x, y = step_2d_logistic(x, y, mu, gamma)
    # Fisher-Yates driven by successive x-coordinates of simultaneous map states.
    for i in range(target_len - 1, 0, -1):
        x, y = step_2d_logistic(x, y, mu, gamma)
        key = x
        j = int(np.floor(key * (i + 1)))
        if j < 0:
            j = 0
        elif j > i:
            j = i
        perm[i], perm[j] = perm[j], perm[i]
    return perm


def map_to_sbox_1d(phase_x, mu, transient=500, target_len=256):
    """Construct a permutation from a 1D map using the shared phase-gene transform."""
    perm = np.arange(target_len, dtype=np.int64)
    x = float(np.sin(np.pi * phase_x) ** 2)
    for _ in range(transient):
        x = (mu * x * (1 - x)) % 1.0
    for index in range(target_len - 1, 0, -1):
        x = (mu * x * (1 - x)) % 1.0
        swap_index = int(np.floor(x * (index + 1)))
        perm[index], perm[swap_index] = perm[swap_index], perm[index]
    return perm


def aes_sbox():
    """Generate the standard AES S-box from GF(2^8) inversion and affine map."""
    def gf_mul(a, b):
        result = 0
        for _ in range(8):
            if b & 1:
                result ^= a
            high = a & 0x80
            a = (a << 1) & 0xFF
            if high:
                a ^= 0x1B
            b >>= 1
        return result

    def gf_pow(a, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = gf_mul(result, a)
            a = gf_mul(a, a)
            exponent >>= 1
        return result

    def rotl8(value, count):
        return ((value << count) | (value >> (8 - count))) & 0xFF

    table = np.empty(256, dtype=np.int64)
    for value in range(256):
        inverse = 0 if value == 0 else gf_pow(value, 254)
        table[value] = (inverse ^ rotl8(inverse, 1) ^ rotl8(inverse, 2)
                        ^ rotl8(inverse, 3) ^ rotl8(inverse, 4) ^ 0x63)
    return table


# ----------------------- Pre-computed parity tables ----------------------
_X = np.arange(256, dtype=np.int64)
_BITS_X = np.stack([(_X >> k) & 1 for k in range(8)], axis=0).astype(np.int64)
_X_PM = 1 - 2 * _BITS_X  # 8x256, +1/-1 per bit


# ---------------------------- Metrics ------------------------------------
def walsh_hadamard(f):
    n = int(round(np.log2(len(f))))
    h = f.astype(np.int64).copy()
    for i in range(n):
        step = 1 << i
        for j in range(0, len(h), 2 * step):
            for k in range(step):
                u = h[j+k]; v = h[j+k+step]
                h[j+k] = u + v
                h[j+k+step] = u - v
    return h


def _validated_sbox(sbox):
    values = np.asarray(sbox)
    if values.shape != (256,) or not np.issubdtype(values.dtype, np.integer):
        raise ValueError("an 8-bit S-box must contain 256 integer outputs")
    if np.any(values < 0) or np.any(values > 255) or np.unique(values).size != 256:
        raise ValueError("an 8-bit S-box must be a permutation of 0..255")
    return values.astype(np.int64, copy=False)


def _component_spectra(sbox, masks):
    components = np.zeros((len(masks), 256), dtype=np.int64)
    for bit in range(8):
        selected = np.array([(mask >> bit) & 1 for mask in masks], dtype=bool)
        if np.any(selected):
            components[selected] ^= ((sbox >> bit) & 1)
    spectra = 1 - 2 * components
    width = 1
    while width < 256:
        blocks = spectra.reshape(len(masks), -1, 2, width)
        left = blocks[:, :, 0, :].copy()
        right = blocks[:, :, 1, :].copy()
        blocks[:, :, 0, :] = left + right
        blocks[:, :, 1, :] = left - right
        width *= 2
    return spectra


def component_nonlinearity(sbox, output_mask):
    """Return the Boolean nonlinearity of the component output_mask · S(x)."""
    values = _validated_sbox(sbox)
    if not isinstance(output_mask, (int, np.integer)) or not 1 <= output_mask <= 255:
        raise ValueError("output_mask must be an integer in 1..255")
    spectrum = np.abs(_component_spectra(values, [int(output_mask)]))
    return 128 - int(spectrum.max()) // 2


def vectorial_nonlinearity(sbox):
    """Return minimum component nonlinearity over all 255 nonzero output masks."""
    values = _validated_sbox(sbox)
    spectra = np.abs(_component_spectra(values, range(1, 256)))
    component_maxima = spectra.max(axis=1)
    return int(np.min(128 - component_maxima // 2))


def random_permutation_sample(seed, count):
    """Return seeded uniform random permutations as a descriptive reference."""
    rng = np.random.default_rng(seed)
    return [rng.permutation(256) for _ in range(count)]


def random_parameter_candidates(seed, budget):
    """Return a seeded uniform sample of valid 2D-map generator parameters."""
    rng = np.random.default_rng(seed)
    candidates = []
    for _ in range(budget):
        candidates.append((
            float(rng.uniform(0.0, 1.0)),
            float(rng.uniform(0.0, 1.0)),
            float(rng.uniform(3.57, 4.0)),
            float(rng.uniform(0.0, 0.5)),
        ))
    return candidates


def normalize_objectives(candidate):
    """Map objectives monotonically into fixed [0,1) minimization coordinates."""
    le = float(candidate["le"])
    nl = float(candidate["nl_min"])
    du = float(candidate["du"])
    if not np.isfinite(le) or not np.isfinite(nl) or not np.isfinite(du):
        raise ValueError("objective values must be finite")
    return (0.5 - np.arctan(le) / np.pi, (128.0 - nl) / 128.0, du / 256.0)


def exact_sign_flip_pvalue(differences, tolerance=1e-12):
    """Compute an exact two-sided paired sign-flip p-value by meet-in-the-middle."""
    values = np.asarray(differences, dtype=float)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("differences must be a non-empty finite one-dimensional array")
    split = values.size // 2
    left = values[:split]
    right = values[split:]

    def signed_sums(part):
        sums = np.array([0.0])
        for value in part:
            sums = np.concatenate((sums + value, sums - value))
        return sums

    left_sums = signed_sums(left)
    right_sums = np.sort(signed_sums(right))
    observed = abs(float(values.sum()))
    threshold = observed - tolerance
    if threshold <= 0:
        return 1.0
    total = 0
    for left_sum in left_sums:
        lower_count = np.searchsorted(right_sums, -threshold - left_sum, side="right")
        upper_count = right_sums.size - np.searchsorted(right_sums, threshold - left_sum, side="left")
        total += lower_count + upper_count
    return float(total / (2 ** values.size))


def hypervolume_3d(front, reference_point):
    """Compute exact minimization hypervolume for a 3D front and fixed reference."""
    points = np.asarray(front, dtype=float)
    reference = np.asarray(reference_point, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or reference.shape != (3,):
        raise ValueError("front and reference_point must be 3-dimensional")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(reference)):
        raise ValueError("front and reference_point must contain finite coordinates")
    if np.any(points >= reference):
        raise ValueError("reference_point must be strictly worse than every front point")
    nondominated = []
    for point in points:
        if not any(np.all(other <= point) and np.any(other < point) for other in points):
            nondominated.append(point)
    xs = sorted({float(p[0]) for p in nondominated} | {float(reference[0])})
    volume = 0.0
    for left, right in zip(xs[:-1], xs[1:]):
        active = [p[1:] for p in nondominated if p[0] <= left]
        if active:
            volume += (right - left) * _hypervolume_2d(active, reference[1:])
    return float(volume)


def _hypervolume_2d(points, reference):
    points = np.asarray(points, dtype=float)
    ys = sorted({float(p[0]) for p in points} | {float(reference[0])})
    area = 0.0
    for low, high in zip(ys[:-1], ys[1:]):
        active_z = [p[1] for p in points if p[0] <= low]
        if active_z:
            area += (high - low) * (reference[1] - min(active_z))
    return float(area)


def coordinate_nonlinearity(sbox):
    """Return coordinate-function NL minimum and mean as diagnostic values."""
    values = _validated_sbox(sbox)
    nls = [component_nonlinearity(values, 1 << bit) for bit in range(8)]
    return int(min(nls)), float(np.mean(nls))


def nonlinearity(sbox):
    """Return minimum and mean coordinate-function nonlinearity."""
    return coordinate_nonlinearity(sbox)


def sac_matrix(sbox):
    A = np.zeros((8, 8))
    sb = sbox[_X]
    for i in range(8):
        diff = sb ^ sbox[_X ^ (1 << i)]
        for j in range(8):
            A[i, j] = ((diff >> j) & 1).mean()
    return A


def sac(sbox):
    return float(sac_matrix(sbox).mean())


def bic_sac(sbox):
    sb = sbox[_X]
    vals = []
    for i in range(8):
        diff = sb ^ sbox[_X ^ (1 << i)]
        for j in range(8):
            for k in range(j+1, 8):
                v = ((diff >> j) & 1) ^ ((diff >> k) & 1)
                vals.append(v.mean())
    arr = np.array(vals)
    return float(np.mean(np.abs(arr - 0.5)))


def differential_uniformity(sbox):
    sb = sbox[_X]
    ddt = np.zeros((256, 256), dtype=np.int32)
    for dx in range(256):
        dy = sb ^ sbox[_X ^ dx]
        binc = np.bincount(dy, minlength=256)
        ddt[dx] = binc
    ddt[0, 0] = 0
    return int(ddt.max())


def lap(sbox):
    """Maximum probability bias |Pr[a.x=b.S(x)] - 1/2| over nonzero masks."""
    sx = sbox[_X]
    bits_sx = np.stack([(sx >> k) & 1 for k in range(8)], axis=0)
    sx_pm = 1 - 2 * bits_sx  # 8x256, +1/-1 per bit of S(x)
    n = 256
    pa = np.empty((255, n), dtype=np.int64)
    pb = np.empty((255, n), dtype=np.int64)
    for i, mask in enumerate(range(1, 256)):
        v = np.ones(n, dtype=np.int64)
        for k in range(8):
            if (mask >> k) & 1:
                v *= _X_PM[k]
        pa[i] = v
        v2 = np.ones(n, dtype=np.int64)
        for k in range(8):
            if (mask >> k) & 1:
                v2 *= sx_pm[k]
        pb[i] = v2
    M = pa @ pb.T  # 255x255
    absolute_correlation = np.abs(M) / n
    return float(absolute_correlation.max() / 2.0)


def shannon_entropy_bits(sbox):
    counts = np.bincount(sbox, minlength=256)
    p = counts / counts.sum()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def npcr_uaci(sbox, image_seed=2026, size=64):
    """NPCR/UACI for a chaotic image encryption pipeline.

    Standard procedure: encrypt the same plaintext under two cipher keystreams
    that differ by a 1-bit key perturbation. The keystream is generated by the
    2D Logistic map (different initial condition -> different orbit) and the
    S-box is applied per byte to produce the ciphertext. We report NPCR and
    UACI averaged across the full image.
    """
    rng = np.random.default_rng(image_seed)
    P = rng.integers(0, 256, size=(size, size), dtype=np.uint8)
    n_bytes = size * size

    def keystream(x0, y0, mu, gamma, n):
        ks = np.empty(n, dtype=np.int64)
        x, y = x0, y0
        for _ in range(500):
            x = mu * x * (1 - x) + gamma * y
            y = mu * y * (1 - y) + gamma * x
            x -= np.floor(x); y -= np.floor(y)
        for i in range(n):
            x = mu * x * (1 - x) + gamma * y
            x -= np.floor(x)
            ks[i] = int(np.floor(x * 256)) & 0xFF
            x_new = mu * x * (1 - x) + gamma * y
            y_new = mu * y * (1 - y) + gamma * x
            x, y = x_new - np.floor(x_new), y_new - np.floor(y_new)
        return ks

    K1 = keystream(0.3, 0.7, 3.91, 0.21, n_bytes)
    # 1-bit perturbation: change initial condition x0 by 1 ULP (1e-9)
    K2 = keystream(0.3 + 1e-9, 0.7, 3.91, 0.21, n_bytes)
    # Round 1: S-box(plaintext) XOR keystream
    C1 = sbox[P.reshape(-1)].reshape(size, size).astype(np.uint8) ^ K1.reshape(size, size).astype(np.uint8)
    C2 = sbox[P.reshape(-1)].reshape(size, size).astype(np.uint8) ^ K2.reshape(size, size).astype(np.uint8)
    diff = (C1 != C2).astype(np.float32)
    npcr = diff.mean() * 100
    uaci = (np.abs(C1.astype(np.float32) - C2.astype(np.float32)) / 255).mean() * 100
    return float(npcr), float(uaci)


# ---------------------------- NSGA-II ------------------------------------
def objective_keys():
    """The independent objectives actually optimized by the algorithm."""
    return ("le", "nl_min", "du")


def objective_vector(candidate):
    """Map one candidate to the ordered minimization vector used by NSGA-II."""
    return (-float(candidate["le"]), -float(candidate["nl_min"]), float(candidate["du"]))


def select_knee(front):
    """Select the normalized Euclidean-distance-to-ideal compromise point."""
    if not front:
        raise ValueError("Cannot select a compromise from an empty Pareto front")
    values = -np.array([objective_vector(candidate) for candidate in front], dtype=float)
    low = values.min(axis=0)
    span = values.max(axis=0) - low
    normalized = np.divide(values - low, span, out=np.zeros_like(values), where=span > 0)
    distances = np.linalg.norm(1.0 - normalized, axis=1)
    return int(np.argmin(distances))


def dominates(a, b):
    """Return whether candidate a is no worse on every objective and better on one."""
    a_values = objective_vector(a)
    b_values = objective_vector(b)
    return all(av <= bv for av, bv in zip(a_values, b_values)) and any(
        av < bv for av, bv in zip(a_values, b_values)
    )


def fast_non_dominated_sort(pop):
    n_pop = len(pop)
    S = [set() for _ in range(n_pop)]
    n_dom = [0] * n_pop
    fronts = [[]]
    for p in range(n_pop):
        for q in range(n_pop):
            if p == q:
                continue
            if dominates(pop[p], pop[q]):
                S[p].add(q)
            elif dominates(pop[q], pop[p]):
                n_dom[p] += 1
        if n_dom[p] == 0:
            fronts[0].append(p)
    i = 0
    while fronts[i]:
        Q = []
        for p in fronts[i]:
            for q in S[p]:
                n_dom[q] -= 1
                if n_dom[q] == 0:
                    Q.append(q)
        i += 1
        fronts.append(Q)
    return fronts


def tournament_winner(a, b, ranks, crowding):
    """Return the rank/crowding winner of one binary tournament."""
    if ranks[a] < ranks[b]:
        return a
    if ranks[b] < ranks[a]:
        return b
    if crowding[a] > crowding[b]:
        return a
    if crowding[b] > crowding[a]:
        return b
    return min(a, b)


def crowding_distance(pop, front):
    n = len(front)
    if n == 0:
        return {}
    dist = {i: 0.0 for i in front}
    values = np.asarray([objective_vector(pop[index]) for index in front], dtype=float)
    for column_index in range(values.shape[1]):
        order = np.argsort(values[:, column_index])
        sf = [front[index] for index in order]
        dist[sf[0]] = float("inf")
        dist[sf[-1]] = float("inf")
        vmin = values[order[0], column_index]
        vmax = values[order[-1], column_index]
        if vmax == vmin:
            continue
        for position in range(1, n - 1):
            before = values[order[position - 1], column_index]
            after = values[order[position + 1], column_index]
            dist[sf[position]] += (after - before) / (vmax - vmin)
    return dist


def _evaluate_1d_candidate(phase_x, mu):
    x0 = float(np.sin(np.pi * phase_x) ** 2)
    sbox = map_to_sbox_1d(phase_x, mu)
    return {
        "phase_x": float(phase_x), "phase_y": 0.0, "x0": x0, "y0": 0.0,
        "mu": float(mu), "gamma": 0.0,
        "le": lyapunov_1d_logistic(x0, mu, n=1500, transient=200),
        "nl_min": vectorial_nonlinearity(sbox),
        "nl_avg": coordinate_nonlinearity(sbox)[1],
        "sac": sac(sbox), "bic_sac": bic_sac(sbox),
        "du": differential_uniformity(sbox), "lap": lap(sbox),
        "entropy": shannon_entropy_bits(sbox), "sbox": sbox.tolist(),
    }


def _canonical_config_hash(config):
    return hashlib.sha256(json.dumps(config, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _result_record(candidate, variant, seed, config_hash):
    return {
        **candidate,
        "metric_definition": "vectorial_nonlinearity_min_over_255_nonzero_output_masks",
        "variant": variant,
        "seed": int(seed),
        "config_hash": config_hash,
        "objective_values": list(objective_vector(candidate)),
        "normalized_objectives": list(normalize_objectives(candidate)),
    }


def _build_run_result(evaluations, variant, seed, config, evaluation_count):
    if not evaluations:
        raise ValueError("at least one candidate evaluation is required")
    if evaluation_count != len(evaluations):
        raise ValueError("evaluation_count must equal the number of evaluated candidates")
    if evaluation_count != int(config.get("evaluation_budget", evaluation_count)):
        raise ValueError("run did not consume the configured exact evaluation budget")
    fronts = fast_non_dominated_sort(evaluations)
    front = [evaluations[index] for index in fronts[0]]
    selected = front[select_knee(front)]
    config_hash = _canonical_config_hash(config)
    records = [_result_record(candidate, variant, seed, config_hash) for candidate in evaluations]
    front_records = [_result_record(candidate, variant, seed, config_hash) for candidate in front]
    return {
        "provenance": {"simulated": True, "seed": int(seed), "fpa_enabled": variant == "fpa",
                       "notes": "Algorithm-generated simulation; not hardware or full-cipher measurement."},
        "metric_definition": "vectorial_nonlinearity_min_over_255_nonzero_output_masks",
        "evaluation_count": int(evaluation_count), "seed": int(seed), "variant": variant,
        "config_hash": config_hash, "objective_values": list(objective_vector(selected)),
        "normalized_objectives": list(normalize_objectives(selected)),
        "initial_population": config.get("initial_population_records"),
        "hypervolume": hypervolume_3d(
            [normalize_objectives(candidate) for candidate in front], [1.1, 1.1, 1.1]
        ),
        "selected": {key: value for key, value in selected.items() if key != "sbox"},
        "sbox": selected["sbox"], "pareto_size": len(front),
        "pareto": [{key: value for key, value in record.items() if key != "sbox"} for record in front_records],
        "initial_selected": config.get("initial_selected_record"),
        "initial_pareto_size": config.get("initial_pareto_size"),
        "evaluations": records,
        "initial_population": config.get("initial_population_records", []),
        "initial_selected": config.get("initial_selected_record"),
        "initial_pareto_size": config.get("initial_pareto_size"),
        "pop_size": config.get("pop_size"),
        "n_gen": config.get("n_gen"), "fpa_interval": config.get("fpa_interval"),
        "evaluation_budget": config.get("evaluation_budget"),
    }


def run_optimizer(seed, fpa_enabled, pop_size, n_gen, fpa_interval, evaluation_budget, initial_population=None):
    """Run NSGA-II with an exact objective-call budget; FPA calls consume budget."""
    if pop_size < 2 or n_gen < 1 or evaluation_budget < pop_size:
        raise ValueError("require pop_size >= 2, n_gen >= 1, and budget >= population size")
    np.random.seed(int(seed))
    variant = "fpa" if fpa_enabled else "nsga_only"
    config = {"seed": int(seed), "variant": variant, "pop_size": int(pop_size),
              "n_gen": int(n_gen), "fpa_interval": int(fpa_interval),
              "evaluation_budget": int(evaluation_budget),
              "metric_definition": "vectorial_nonlinearity_min_over_255_nonzero_output_masks",
              "objective_keys": list(objective_keys())}
    evaluated = []
    initial_candidates = []
    if initial_population is not None and len(initial_population) != pop_size:
        raise ValueError("initial_population length must match pop_size")
    for index in range(min(pop_size, evaluation_budget)):
        if initial_population is None:
            chromosome = [float(np.random.uniform(0, 1)), float(np.random.uniform(0, 1)),
                          float(np.random.uniform(3.57, 4)), float(np.random.uniform(0, .5))]
        else:
            chromosome = list(initial_population[index])
        candidate = evaluate_chromosome(chromosome)
        candidate["generation"] = 0
        initial_candidates.append(candidate)
    evaluated.extend(initial_candidates)
    if not evaluated:
        raise ValueError("evaluation_budget must be positive")
    calls = len(evaluated)
    population = evaluated[:]
    for generation in range(n_gen):
        if calls >= evaluation_budget:
            break
        fronts = fast_non_dominated_sort(population)
        ranks = [0] * len(population)
        crowding = {}
        for rank, front in enumerate(fronts):
            for index in front:
                ranks[index] = rank
            crowding.update(crowding_distance(population, front))
        if fpa_enabled and fpa_interval > 0 and (generation + 1) % fpa_interval == 0:
            candidates = sorted(fronts[0], key=lambda index: -crowding.get(index, 0.0))[:10]
            for index in candidates:
                if calls >= evaluation_budget:
                    break
                updated = evaluate_chromosome(fpa_mutate(population, index))
                updated["generation"] = generation
                calls += 1
                population[index] = updated
                evaluated.append(updated)
        offspring = []
        target_offspring = min(pop_size // 2, evaluation_budget - calls)
        while len(offspring) < target_offspring:
            a1, a2 = np.random.choice(len(population), 2, replace=len(population) < 2)
            b1, b2 = np.random.choice(len(population), 2, replace=len(population) < 2)
            parent_a = population[tournament_winner(int(a1), int(a2), ranks, crowding)]
            parent_b = population[tournament_winner(int(b1), int(b2), ranks, crowding)]
            genes_a = [parent_a[key] for key in ("phase_x", "phase_y", "mu", "gamma")]
            genes_b = [parent_b[key] for key in ("phase_x", "phase_y", "mu", "gamma")]
            mask = np.random.random(4) > .5
            child = [genes_a[i] if mask[i] else genes_b[i] for i in range(4)]
            for gene in range(4):
                child[gene] += np.random.normal(0, .01)
            evaluated_child = evaluate_chromosome(wrap_parameters(child))
            evaluated_child["generation"] = generation + 1
            calls += 1
            offspring.append(evaluated_child)
            evaluated.append(evaluated_child)
        combined = population + offspring
        fronts = fast_non_dominated_sort(combined)
        survivors = []
        for front in fronts:
            if len(survivors) + len(front) <= pop_size:
                survivors.extend(combined[index] for index in front)
            else:
                distances = crowding_distance(combined, front)
                order = sorted(front, key=lambda index: -distances.get(index, 0.0))
                survivors.extend(combined[index] for index in order[:pop_size-len(survivors)])
                break
        population = survivors
    while calls < evaluation_budget:
        fronts = fast_non_dominated_sort(population)
        ranks = [0] * len(population)
        crowding = {}
        for rank, front in enumerate(fronts):
            for index in front:
                ranks[index] = rank
            crowding.update(crowding_distance(population, front))
        a1, a2 = np.random.choice(len(population), 2, replace=len(population) < 2)
        b1, b2 = np.random.choice(len(population), 2, replace=len(population) < 2)
        parent_a = population[tournament_winner(int(a1), int(a2), ranks, crowding)]
        parent_b = population[tournament_winner(int(b1), int(b2), ranks, crowding)]
        genes_a = [parent_a[key] for key in ("phase_x", "phase_y", "mu", "gamma")]
        genes_b = [parent_b[key] for key in ("phase_x", "phase_y", "mu", "gamma")]
        mask = np.random.random(4) > .5
        child = [genes_a[i] if mask[i] else genes_b[i] for i in range(4)]
        for gene in range(4):
            child[gene] += np.random.normal(0, .01)
        candidate = evaluate_chromosome(wrap_parameters(child))
        candidate["generation"] = n_gen
        calls += 1
        evaluated.append(candidate)
        combined = population + [candidate]
        fronts = fast_non_dominated_sort(combined)
        survivors = []
        for front in fronts:
            if len(survivors) + len(front) <= pop_size:
                survivors.extend(combined[index] for index in front)
            else:
                distances = crowding_distance(combined, front)
                order = sorted(front, key=lambda index: -distances.get(index, 0.0))
                survivors.extend(combined[index] for index in order[:pop_size-len(survivors)])
                break
        population = survivors
    config["initial_population_records"] = [
        _result_record(candidate, variant, seed, _canonical_config_hash(config))
        for candidate in initial_candidates
    ]
    result = _build_run_result(evaluated, variant, seed, config, calls)
    initial_front = fast_non_dominated_sort(initial_candidates)[0]
    initial_pareto = [initial_candidates[index] for index in initial_front]
    initial_selected = initial_pareto[select_knee(initial_pareto)]
    config["initial_selected_record"] = _result_record(initial_selected, variant, seed, result["config_hash"])
    config["initial_pareto_size"] = len(initial_pareto)
    result["initial_selected"] = config["initial_selected_record"]
    result["initial_pareto_size"] = len(initial_pareto)
    result["initial_population"] = [
        _result_record(candidate, variant, seed, result["config_hash"])
        for candidate in initial_candidates
    ]
    result["initial_selected"] = _result_record(initial_selected, variant, seed, result["config_hash"])
    result["initial_pareto_size"] = len(initial_pareto)
    config["initial_population_records"] = result["initial_population"]
    config["initial_selected_record"] = result["initial_selected"]
    config["initial_pareto_size"] = result["initial_pareto_size"]
    return result


def _write_json_atomic(path, payload):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    temporary.replace(destination)


def execute_run_node(run_root, method, seed, protocol_hash, evaluation_budget=940, retry=0):
    """Run one attempt of a frozen-protocol method/seed cell."""
    if retry not in (0, 1):
        raise ValueError("retry must be zero or one")
    if protocol_hash != _file_sha256(Path(run_root) / "protocol.md"):
        raise ValueError("protocol hash does not match frozen protocol file")
    run_root = Path(run_root)
    run_root.mkdir(parents=True, exist_ok=True)
    protocol_copy = run_root / "protocol.md"
    if not protocol_copy.exists():
        raise FileNotFoundError(f"frozen protocol copy is missing: {protocol_copy}")
    suffix = "" if retry == 0 else "_retry1"
    result_path = run_root / "results" / f"seed_{method}_{seed}{suffix}.json"
    log_path = run_root / "logs" / f"seed_{method}_{seed}{suffix}.jsonl"
    node_path = run_root / "checkpoints" / f"N4_{method}_{seed}{suffix}.json"
    run_root.joinpath("results").mkdir(parents=True, exist_ok=True)
    run_root.joinpath("logs").mkdir(parents=True, exist_ok=True)
    run_root.joinpath("checkpoints").mkdir(parents=True, exist_ok=True)
    if result_path.exists() or node_path.exists():
        raise FileExistsError(f"run node already exists: {result_path}")
    source_hash = _file_sha256(Path(__file__))
    config = {"method": method, "seed": int(seed), "evaluation_budget": int(evaluation_budget),
              "protocol_hash": protocol_hash, "simulator_hash": source_hash}
    config_hash = _canonical_config_hash(config)
    node = {"node_id": f"N4:{method}:{seed}", "status": "running", "inputs": config,
            "input_sha256": {"protocol": protocol_hash, "simulator": source_hash},
            "outputs": [str(result_path), str(log_path)], "output_sha256": None,
            "started_at": time.time(), "finished_at": None, "command": "execute_run_node",
            "evaluation_count": 0, "exit_status": None, "retry": int(retry)}
    _write_json_atomic(node_path, node)
    with log_path.open("a", encoding="utf-8") as log:
        log.write(json.dumps({"event": "start", "node": node["node_id"], "ts": node["started_at"], "config": config}) + "\n")
    try:
        if method in ("nsga_only", "fpa"):
            if evaluation_budget < 80:
                raise ValueError("NSGA-II evaluation budget must be at least the population size")
            initial_rng = np.random.default_rng(seed)
            initial_population = [
                [float(initial_rng.uniform(0, 1)), float(initial_rng.uniform(0, 1)),
                 float(initial_rng.uniform(3.57, 4)), float(initial_rng.uniform(0, .5))]
                for _ in range(80)
            ]
            result = run_optimizer(seed, method == "fpa", 80, 20, 3, evaluation_budget,
                                   initial_population=initial_population)
        elif method == "random_search":
            result = run_random_search(seed, evaluation_budget)
        elif method == "1d_random_search":
            result = run_1d_random_search(seed, evaluation_budget)
        elif method == "random_permutation_reference":
            result = run_random_permutation_reference(seed, evaluation_budget)
        else:
            raise ValueError(f"unknown run method: {method}")
        if result.get("evaluation_count") != evaluation_budget:
            raise ValueError(f"expected {evaluation_budget} objective calls, got {result.get('evaluation_count')}")
        result.update({"run_protocol_hash": protocol_hash, "run_config_hash": config_hash,
                       "simulator_hash": source_hash, "source_revision": source_hash})
        _write_json_atomic(result_path, result)
        validate_run_artifact(result_path, protocol_hash, source_hash, evaluation_budget, int(seed),
                              "random_permutation_reference" if method == "random_permutation_reference" else method)
        node.update({"status": "complete", "finished_at": time.time(), "exit_status": 0,
                     "evaluation_count": result.get("evaluation_count", result.get("permutation_count")),
                     "output_sha256": {"result": _file_sha256(result_path), "log": _file_sha256(log_path)}})
        with log_path.open("a", encoding="utf-8") as log:
            log.write(json.dumps({"event": "complete", "ts": node["finished_at"],
                                  "evaluation_count": node["evaluation_count"],
                                  "elapsed_seconds": node["finished_at"]-node["started_at"]}) + "\n")
        _write_json_atomic(node_path, node)
        return node
    except Exception as error:
        if result_path.exists():
            result_path.unlink()
        node.update({"status": "failed", "finished_at": time.time(), "exit_status": 1,
                     "error": f"{type(error).__name__}: {error}"})
        with log_path.open("a", encoding="utf-8") as log:
            log.write(json.dumps({"event": "failed", "ts": node["finished_at"], "error": node["error"]}) + "\n")
        _write_json_atomic(node_path, node)
        raise


def validate_run_artifact(path, protocol_hash, simulator_hash, evaluation_budget, seed, method):
    """Validate provenance and metric ranges for one completed run artifact."""
    result = json.loads(Path(path).read_text(encoding="utf-8"))
    required = {"run_protocol_hash", "run_config_hash", "simulator_hash", "evaluation_count",
                "metric_definition", "seed", "variant"}
    if method != "random_permutation_reference":
        required |= {"selected", "objective_values"}
    else:
        required.add("metrics")
    missing = required - set(result)
    if missing:
        raise ValueError(f"run artifact missing fields: {sorted(missing)}")
    if result["run_protocol_hash"] != protocol_hash or result["simulator_hash"] != simulator_hash:
        raise ValueError("run artifact provenance hash mismatch")
    if result["evaluation_count"] != evaluation_budget or result["seed"] != seed:
        raise ValueError("run artifact evaluation count or seed mismatch")
    if result["variant"] != method or result["metric_definition"] != "vectorial_nonlinearity_min_over_255_nonzero_output_masks":
        raise ValueError("run artifact method or metric definition mismatch")
    if method == "random_permutation_reference":
        if result.get("permutation_count") != evaluation_budget or len(result["metrics"]) != evaluation_budget:
            raise ValueError("random-permutation sample count mismatch")
        for metrics in result["metrics"]:
            if not 0 <= metrics["nl_min"] <= 128 or not 2 <= metrics["du"] <= 256:
                raise ValueError("random-permutation artifact contains out-of-domain metrics")
    else:
        if not np.isfinite(result["selected"]["le"]) or not 0 <= result["selected"]["nl_min"] <= 128 or not 2 <= result["selected"]["du"] <= 256:
            raise ValueError("run artifact contains out-of-domain selected metrics")
    if result.get("hypervolume") is not None and (not np.isfinite(result["hypervolume"]) or result["hypervolume"] < 0):
        raise ValueError("run artifact hypervolume is invalid")
    return result


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_random_search(seed, budget):
    """Uniform random search over valid 2D-map parameters with exact call budget."""
    if budget < 1:
        raise ValueError("budget must be positive")
    rng = np.random.default_rng(seed)
    config = {"seed": int(seed), "variant": "random_search", "evaluation_budget": int(budget),
              "metric_definition": "vectorial_nonlinearity_min_over_255_nonzero_output_masks",
              "objective_keys": list(objective_keys())}
    evaluations = []
    for _ in range(budget):
        candidate = [float(rng.uniform(0, 1)), float(rng.uniform(0, 1)),
                     float(rng.uniform(3.57, 4)), float(rng.uniform(0, .5))]
        evaluations.append(evaluate_chromosome(candidate))
    return _build_run_result(evaluations, "random_search", seed, config, len(evaluations))


def run_1d_random_search(seed, budget):
    """Uniform parameter search for the x-only 1D map with exact call budget."""
    if budget < 1:
        raise ValueError("budget must be positive")
    rng = np.random.default_rng(seed)
    config = {"seed": int(seed), "variant": "1d_random_search", "evaluation_budget": int(budget),
              "metric_definition": "vectorial_nonlinearity_min_over_255_nonzero_output_masks",
              "objective_keys": list(objective_keys()), "map_dimension": 1}
    evaluations = []
    for _ in range(budget):
        candidate = _evaluate_1d_candidate(float(rng.uniform(0, 1)), float(rng.uniform(3.57, 4)))
        evaluations.append(candidate)
    return _build_run_result(evaluations, "1d_random_search", seed, config, len(evaluations))


def run_random_permutation_reference(seed, count):
    """Generate a descriptive seed-stable metric distribution over uniform permutations."""
    if count < 1:
        raise ValueError("count must be positive")
    records = []
    for index, sbox in enumerate(random_permutation_sample(seed, count)):
        metrics = {"nl_min": vectorial_nonlinearity(sbox), "nl_avg": coordinate_nonlinearity(sbox)[1],
                   "du": differential_uniformity(sbox), "sac": sac(sbox), "bic_sac": bic_sac(sbox),
                   "lap": lap(sbox), "sbox": sbox.tolist()}
        records.append({"seed": int(seed), "sample_index": index, **metrics})
    return {"seed": int(seed), "permutation_count": int(count), "evaluation_count": int(count),
            "hypervolume": None, "variant": "random_permutation_reference",
            "metric_definition": "vectorial_nonlinearity_min_over_255_nonzero_output_masks", "metrics": records}


def evaluate_chromosome(chrom, include_cipher_metrics=False):
    x_gene, y_gene, mu, gamma = chrom
    x0, y0 = initial_state(x_gene, y_gene)
    sbox = map_to_sbox(x0, y0, mu, gamma)
    le = lyapunov_2d_logistic(x0, y0, mu, gamma, n=1500, transient=200)
    nl_min = vectorial_nonlinearity(sbox)
    _, nl_avg = coordinate_nonlinearity(sbox)
    sac_v = sac(sbox)
    bic_s = bic_sac(sbox)
    du = differential_uniformity(sbox)
    lap_v = lap(sbox)
    h = shannon_entropy_bits(sbox)
    if include_cipher_metrics:
        npcr, uaci = npcr_uaci(sbox)
    result = dict(
        phase_x=float(x_gene), phase_y=float(y_gene),
        x0=float(x0), y0=float(y0), mu=float(mu), gamma=float(gamma),
        le=float(le), nl_min=int(nl_min), nl_avg=float(nl_avg),
        sac=float(sac_v), bic_sac=float(bic_s),
        du=int(du), lap=float(lap_v), entropy=float(h),
        sbox=sbox.tolist())
    if include_cipher_metrics:
        result.update(npcr=float(npcr), uaci=float(uaci))
    return result


def levy_flight(Lambda=1.5, size=4):
    sigma_u = (math.gamma(1 + Lambda) * np.sin(np.pi * Lambda / 2) /
               (math.gamma((1 + Lambda) / 2) * Lambda * 2 ** ((Lambda - 1) / 2))) ** (1 / Lambda)
    u = np.random.normal(0, sigma_u, size=size)
    v = np.random.normal(0, 1, size=size)
    step = u / (np.abs(v) ** (1 / Lambda))
    return 0.01 * step


def fpa_mutate(pop, idx, p_switch=0.8):
    new = np.array([pop[idx]["phase_x"], pop[idx]["phase_y"], pop[idx]["mu"], pop[idx]["gamma"]])
    if np.random.random() < p_switch:
        others = [j for j in range(len(pop)) if j != idx]
        if others:
            other = np.random.choice(others)
            target = np.array([pop[other]["phase_x"], pop[other]["phase_y"], pop[other]["mu"], pop[other]["gamma"]])
            new = new + levy_flight(size=4) * (new - target)
    else:
        new = new + np.random.normal(0, 0.02, size=4)
    return wrap_parameters(new.tolist())


# ---------------------------- Main run -----------------------------------
def main(seed=SEED, fpa_enabled=True, pop_size=80, n_gen=20, fpa_interval=3, output_path="results.json", generate_figures=True, evaluation_budget=None):
    """Run one reproducible optimizer trial and optionally render its figures."""
    if evaluation_budget is None:
        evaluation_budget = pop_size + n_gen * (pop_size // 2)
        if fpa_enabled and fpa_interval > 0:
            evaluation_budget += (n_gen // fpa_interval) * min(10, pop_size)
    run = run_optimizer(seed, fpa_enabled, pop_size, n_gen, fpa_interval, evaluation_budget)
    Path(output_path).write_text(json.dumps(run, indent=2, default=str), encoding="utf-8")
    print(f"Wrote {output_path} (variant={run['variant']}, evaluations={run['evaluation_count']})")
    if generate_figures:
        warnings.warn("figure regeneration is deferred until the result bundle is frozen", RuntimeWarning)
    return run

def cli_main():
    parser = argparse.ArgumentParser(description="Run the S-box optimizer.")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--no-fpa", action="store_true")
    parser.add_argument("--pop-size", type=int, default=80)
    parser.add_argument("--generations", type=int, default=20)
    parser.add_argument("--fpa-interval", type=int, default=3)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()
    if args.output is None:
        suffix = "nsga_only" if args.no_fpa else "fpa"
        args.output = f"results_{suffix}_seed{args.seed}.json"
    main(seed=args.seed, fpa_enabled=not args.no_fpa, pop_size=args.pop_size,
         n_gen=args.generations, fpa_interval=args.fpa_interval,
         output_path=args.output)


if __name__ == "__main__":
    cli_main()
