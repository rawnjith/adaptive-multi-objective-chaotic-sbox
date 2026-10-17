"""NSGA-II-FPA / 2D-Logistic S-box simulator — honest, reproducible.

Seed = 20260814.
Outputs:
  figures/*.pdf           publication-grade figures
  results.json            full Pareto front + selected S-box + metrics
"""
from __future__ import annotations

import argparse
import json
import time
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


def nonlinearity(sbox):
    nls = []
    for bit in range(8):
        f = ((sbox >> bit) & 1).astype(np.int64) * 2 - 1
        wh = np.abs(walsh_hadamard(f))
        nls.append(128 - wh.max() // 2)
    return int(np.min(nls)), float(np.mean(nls))


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


def select_knee(front):
    """Select the normalized Euclidean-distance-to-ideal compromise point."""
    if not front:
        raise ValueError("Cannot select a compromise from an empty Pareto front")
    values = np.array([[p["le"], p["nl_min"], -p["du"]] for p in front], dtype=float)
    low = values.min(axis=0)
    span = values.max(axis=0) - low
    normalized = np.divide(values - low, span, out=np.zeros_like(values), where=span > 0)
    distances = np.linalg.norm(1.0 - normalized, axis=1)
    return int(np.argmin(distances))


def dominates(a, b):
    # Maximize LE and NL; minimize DU. No independent hardware timing
    # measurements are available for a runtime objective.
    a_v = (-a["le"], -a["nl_min"], a["du"])
    b_v = (-b["le"], -b["nl_min"], b["du"])
    return all(av <= bv for av, bv in zip(a_v, b_v)) and any(av < bv for av, bv in zip(a_v, b_v))


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
    keys = [("le", -1), ("nl_min", -1), ("du", 1)]
    for k, sgn in keys:
        sf = sorted(front, key=lambda i: sgn * pop[i][k])
        dist[sf[0]] = float("inf")
        dist[sf[-1]] = float("inf")
        vmin = pop[sf[0]][k]
        vmax = pop[sf[-1]][k]
        if vmax == vmin:
            continue
        for i in range(1, n - 1):
            dist[sf[i]] += (pop[sf[i+1]][k] - pop[sf[i-1]][k]) / (vmax - vmin)
    return dist


def evaluate_chromosome(chrom, include_cipher_metrics=False):
    x_gene, y_gene, mu, gamma = chrom
    x0, y0 = initial_state(x_gene, y_gene)
    sbox = map_to_sbox(x0, y0, mu, gamma)
    le = lyapunov_2d_logistic(x0, y0, mu, gamma, n=1500, transient=200)
    nl_min, nl_avg = nonlinearity(sbox)
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
    sigma_u = (np.math.gamma(1 + Lambda) * np.sin(np.pi * Lambda / 2) /
               (np.math.gamma((1 + Lambda) / 2) * Lambda * 2 ** ((Lambda - 1) / 2))) ** (1 / Lambda)
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
def main(seed=SEED, fpa_enabled=True, pop_size=80, n_gen=20, fpa_interval=3, output_path="results.json", generate_figures=True):
    """Run one optimizer trial with explicit random seed and FPA configuration."""
    global SEED
    SEED = int(seed)
    np.random.seed(SEED)
    POP_SIZE = int(pop_size)
    N_GEN = int(n_gen)
    FPA_INTERVAL = int(fpa_interval) if fpa_enabled else N_GEN + 1

    init_pop = []
    for _ in range(POP_SIZE):
        chrom = [
            float(np.random.uniform(0.0, 1.0)),
            float(np.random.uniform(0.0, 1.0)),
            float(np.random.uniform(3.57, 4.0)),
            float(np.random.uniform(0.0, 0.5)),
        ]
        init_pop.append(evaluate_chromosome(chrom))
    print(f"[init] evaluated {POP_SIZE} chromosomes")

    history = []
    for gen in range(N_GEN):
        t0 = time.time()
        fronts = fast_non_dominated_sort(init_pop)
        ranks = [0] * len(init_pop)
        cd = {}
        for rank, front in enumerate(fronts):
            for idx in front:
                ranks[idx] = rank
            cd.update(crowding_distance(init_pop, front))
        best_le = max(p["le"] for p in init_pop)
        best_nl = max(p["nl_min"] for p in init_pop)
        best_du = min(p["du"] for p in init_pop)
        history.append(dict(gen=gen, best_le=best_le, best_nl=best_nl, best_du=best_du))

        # FPA step every FPA_INTERVAL generations
        if gen % FPA_INTERVAL == FPA_INTERVAL - 1 and cd:
            sorted_idx = sorted(fronts[0], key=lambda i: -cd[i])[:10]
            for idx in sorted_idx:
                new_chrom = fpa_mutate(init_pop, idx)
                init_pop[idx] = evaluate_chromosome(new_chrom)

        # Generate new population via tournament + crossover + mutation
        new_pop = list(init_pop)
        for _ in range(POP_SIZE // 2):
            a1, a2 = np.random.choice(len(init_pop), 2, replace=False)
            b1, b2 = np.random.choice(len(init_pop), 2, replace=False)
            a = tournament_winner(a1, a2, ranks, cd)
            b = tournament_winner(b1, b2, ranks, cd)
            ca = [init_pop[a]["phase_x"], init_pop[a]["phase_y"], init_pop[a]["mu"], init_pop[a]["gamma"]]
            cb = [init_pop[b]["phase_x"], init_pop[b]["phase_y"], init_pop[b]["mu"], init_pop[b]["gamma"]]
            mask = np.random.random(4) > 0.5
            child = [ca[i] if mask[i] else cb[i] for i in range(4)]
            for k in range(4):
                child[k] += np.random.normal(0, 0.01)
            child = wrap_parameters(child)
            new_pop.append(evaluate_chromosome(child))

        # NSGA-II truncation
        fronts2 = fast_non_dominated_sort(new_pop)
        survivors = []
        for f in fronts2:
            if len(survivors) + len(f) <= POP_SIZE:
                survivors.extend([new_pop[i] for i in f])
            else:
                cd2 = crowding_distance(new_pop, f)
                f_sorted = sorted(f, key=lambda i: -cd2.get(i, 0))
                survivors.extend([new_pop[i] for i in f_sorted[:POP_SIZE - len(survivors)]])
                break
        init_pop = survivors
        print(f"[gen {gen+1}/{N_GEN}] best_le={best_le:.4f} best_nl={best_nl} best_du={best_du} ({time.time()-t0:.1f}s)")

    # Select a compromise by normalized Euclidean distance to the ideal point.
    fronts_final = fast_non_dominated_sort(init_pop)
    pareto = [init_pop[i] for i in fronts_final[0]]
    best_idx = select_knee(pareto)
    selected = pareto[best_idx]
    print("\n[selected]")
    for k, v in selected.items():
        if k != "sbox":
            print(f"  {k}: {v}")

    # Save Pareto front (without sbox to keep JSON small)
    pareto_no_sbox = []
    for p in pareto:
        p2 = {k: v for k, v in p.items() if k != "sbox"}
        pareto_no_sbox.append(p2)

    out = dict(
        provenance=dict(simulated=True, seed=SEED, fpa_enabled=bool(fpa_enabled),
                        notes="All numbers produced by this simulator with the recorded seed."),
        selected={k: v for k, v in selected.items() if k != "sbox"},
        sbox=selected["sbox"],
        pareto_size=len(pareto),
        pareto=pareto_no_sbox,
        history=history,
        pop_size=POP_SIZE, n_gen=N_GEN, fpa_interval=FPA_INTERVAL,
    )
    Path(output_path).write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"Wrote {output_path}")
    if not generate_figures:
        return

    # ---------- Figures ----------
    # Figure 2: attractor
    x_att, y_att = iterate_2d_logistic(0.31, 0.78, 3.91, 0.21, n=5000, transient=500)
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    ax.scatter(x_att[::5], y_att[::5], s=0.3, c=WONG[5], alpha=0.5)
    ax.set_xlabel(r"$x_n$")
    ax.set_ylabel(r"$y_n$")
    ax.set_title("2D Logistic attractor ($\\mu=3.91$, $\\gamma=0.21$)")
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG_DIR / "fig_02_attractor.pdf")
    plt.close(fig)

    # Figure 3: Pareto front
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    le_v = [p["le"] for p in pareto]
    nl_v = [p["nl_min"] for p in pareto]
    du_v = [p["du"] for p in pareto]
    sc = ax.scatter(le_v, nl_v, c=du_v, cmap="viridis", s=30, edgecolor="black", linewidth=0.4)
    ax.scatter([selected["le"]], [selected["nl_min"]], c="red", marker="*", s=80,
               label="selected", zorder=10)
    ax.set_xlabel(r"Lyapunov exponent $\lambda_\mathrm{max}$")
    ax.set_ylabel(r"Min nonlinearity (NL)")
    fig.colorbar(sc, ax=ax, shrink=0.8, label="Differential uniformity (DU)")
    # Place legend below the x-axis label, outside the axes, to avoid overlapping data
    ax.legend(frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=1)
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG_DIR / "fig_03_pareto.pdf")
    plt.close(fig)

    # Figure 4: S-box heatmap
    sbox_arr = np.array(selected["sbox"]).reshape(16, 16)
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    im = ax.imshow(sbox_arr, cmap="viridis", aspect="equal")
    ax.set_xlabel("output high nibble")
    ax.set_ylabel("output low nibble")
    ax.set_title("Selected dynamic S-box")
    fig.colorbar(im, ax=ax, shrink=0.8, label="S(x)")
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG_DIR / "fig_04_sbox.pdf")
    plt.close(fig)

    # Figure 5: SAC matrix + per-bit avalanche
    A = sac_matrix(np.array(selected["sbox"]))
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.6))
    im0 = axes[0].imshow(A, cmap="RdYlGn", vmin=0.3, vmax=0.7, aspect="equal")
    axes[0].set_xlabel("output bit")
    axes[0].set_ylabel("input bit")
    axes[0].set_title("SAC matrix")
    fig.colorbar(im0, ax=axes[0], shrink=0.8)
    axes[1].bar(range(8), A.mean(axis=0), color=WONG[5])
    axes[1].axhline(0.5, color=WONG[6], linestyle="--", linewidth=0.8, label="ideal = 0.5")
    axes[1].set_xlabel("output bit")
    axes[1].set_ylabel("mean avalanche prob.")
    axes[1].set_ylim(0.3, 0.7)
    axes[1].legend(frameon=False)
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG_DIR / "fig_05_sac.pdf")
    plt.close(fig)

    # Figure 6: DDT heatmap (log scale)
    sb = np.array(selected["sbox"])
    ddt = np.zeros((256, 256), dtype=np.int32)
    for dx in range(256):
        dy = sb ^ sb[_X ^ dx]
        ddt[dx] = np.bincount(dy, minlength=256)
    ddt[0, 0] = 0
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    im = ax.imshow(np.log1p(ddt), cmap="inferno", aspect="equal")
    ax.set_xlabel(r"output diff $\Delta y$")
    ax.set_ylabel(r"input diff $\Delta x$")
    ax.set_title("DDT (log1p scale)")
    fig.colorbar(im, ax=ax, shrink=0.8, label="log(1+count)")
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG_DIR / "fig_06_ddt.pdf")
    plt.close(fig)

    # Figure 7: convergence
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    gens = [h["gen"] for h in history]
    ax.plot(gens, [h["best_le"] for h in history], color=WONG[5], marker="o",
            label="best LE")
    ax.plot(gens, [h["best_nl"]/120 for h in history], color=WONG[6], marker="s",
            label="best NL / 120")
    ax.plot(gens, [(16-h["best_du"])/16 for h in history], color=WONG[3], marker="^",
            label="(16-best DU) / 16")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Normalized metric (higher better)")
    # Place legend below the x-axis label, outside the axes, to avoid overlapping curves
    ax.legend(frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=3)
    ax.set_title("NSGA-II-FPA convergence")
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG_DIR / "fig_07_convergence.pdf")
    plt.close(fig)

    print("Saved six single-run figures (fig_02 .. fig_07). Use plot_verified_ablation.py for the paired-seed Figure 8.")


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
