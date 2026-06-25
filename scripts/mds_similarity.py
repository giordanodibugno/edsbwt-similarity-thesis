#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.manifold import MDS

if len(sys.argv) != 2:
    print("Uso: python3 mds_graph.py similarity_matrix.csv")
    sys.exit(1)

INPUT = sys.argv[1]
OUTPUT = INPUT.replace(".csv", "_mds_graph.png")

# =====================
# Lettura matrice
# =====================
df = pd.read_csv(INPUT, index_col=0)
df = df.apply(pd.to_numeric, errors="coerce")

# Assicura stesso ordine righe/colonne
common = df.index.intersection(df.columns)
df = df.loc[common, common]

names = list(df.index)
S = df.values.astype(float)
n = len(names)

# =====================
# Similarità -> distanza
# =====================
D = 1.0 - S
D[np.isnan(D)] = 1.0
np.fill_diagonal(D, 0)

# =====================
# MDS
# =====================
mds = MDS(
    n_components=2,
    dissimilarity="precomputed",
    random_state=42,
    normalized_stress="auto"
)

coords = mds.fit_transform(D)

# Aumenta leggermente la separazione visiva
coords *= 10

# =====================
# Colori per anno
# =====================
year_colors = {
    "19": "tab:blue",
    "20": "tab:green",
    "21": "tab:orange",
    "22": "tab:red",
    "23": "tab:purple"
}

# =====================
# Disegno
# =====================
fig, ax = plt.subplots(figsize=(16, 10))

# Disegna TUTTI gli archi
for i in range(n):
    for j in range(i + 1, n):
        x1, y1 = coords[i]
        x2, y2 = coords[j]

        sim = S[i, j]

        ax.plot(
            [x1, x2],
            [y1, y2],
            color="black",
            alpha=0.03 + 0.25 * sim,
            linewidth=0.10 + 0.90 * sim,
            zorder=1
        )

# Disegna nodi ed etichette
for i, name in enumerate(names):
    x, y = coords[i]

    year = name[:2]
    color = year_colors.get(year, "lightgray")

    ax.scatter(
        x,
        y,
        s=700,
        c=color,
        edgecolors="black",
        linewidths=1,
        zorder=3
    )

    ax.text(
        x,
        y,
        name,
        ha="center",
        va="center",
        fontsize=9,
        fontweight="bold",
        color="white",
        zorder=4
    )

# Legenda
legend_elements = [
    Line2D(
        [0],
        [0],
        marker='o',
        color='w',
        label=year,
        markerfacecolor=color,
        markeredgecolor='black',
        markersize=12
    )
    for year, color in year_colors.items()
]

ax.legend(
    handles=legend_elements,
    title="Anno",
    loc="upper right"
)

ax.set_title("MDS della matrice di similarità")
ax.axis("off")

plt.tight_layout()
plt.savefig(OUTPUT, dpi=300, bbox_inches="tight")
plt.show()

print(f"Grafico salvato in: {OUTPUT}")
