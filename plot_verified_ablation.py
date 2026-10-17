import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA_DIR = Path(os.environ.get("SBOX_RESULTS_DIR", HERE / "results"))
OUT = Path(os.environ.get("SBOX_FIGURE_OUT", HERE / "figures" / "fig_08_paired_ablation.pdf"))
SEEDS = [42, 100, 1234, 2025, 7777]
METRIC = "nl_min"

fpa = np.array([json.loads((DATA_DIR / f"depositfinal_fpa_seed{s}.json").read_text())["selected"][METRIC] for s in SEEDS], float)
nsga = np.array([json.loads((DATA_DIR / f"depositfinal_nsga_only_seed{s}.json").read_text())["selected"][METRIC] for s in SEEDS], float)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "pdf.fonttype": 42,
})
fig, ax = plt.subplots(figsize=(4.8, 3.2), constrained_layout=True)
for i, seed in enumerate(SEEDS):
    ax.plot([0, 1], [nsga[i], fpa[i]], color="#8c8c8c", linewidth=0.8, zorder=1)
ax.scatter(np.zeros(len(SEEDS)), nsga, color="#2a78d6", edgecolor="white", linewidth=0.5,
           s=42, label="NSGA-II only", zorder=2)
ax.scatter(np.ones(len(SEEDS)), fpa, color="#1baf7a", edgecolor="white", linewidth=0.5,
           s=42, label="NSGA-II + FPA", zorder=2)
ax.set_xticks([0, 1], ["NSGA-II only", "NSGA-II + FPA"])
ax.set_ylabel(r"Selected minimum nonlinearity, $NL_{\min}$")
ax.set_xlim(-0.25, 1.25)
ax.set_ylim(min(nsga.min(), fpa.min()) - 1.2, max(nsga.max(), fpa.max()) + 1.2)
ax.grid(axis="y", color="#dddddd", linewidth=0.5)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
fig.savefig(OUT, bbox_inches="tight")
plt.close(fig)
print(OUT)
