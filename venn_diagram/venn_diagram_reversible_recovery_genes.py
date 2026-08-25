# ------------------------------------------------------------------------------
# Script Name: venn_diagram_reversible_recovery_genes.py
# Project: Bulk RNA-seq Analysis of ANDY CD/NFD and young/old WT mice
# Study: NFD recovery
# Description: Venn diagram comparing NFD vs CD DEGs (reversible genes) and
#              Rec vs CD DEGs (recovery-specific genes).
# Dependencies: pandas, matplotlib, matplotlib-venn
# Last Updated: 2026-08-25
# ------------------------------------------------------------------------------

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib_venn import venn2

# NFD vs CD DEG results
nfd_cd_degs_file = "path/to/nfd_cd_degs.csv"  # NFD vs CD DEG results (with Symbol column)
nfd_cd_degs = pd.read_csv(nfd_cd_degs_file)

# Rec vs CD DEG results
rec_cd_degs_file = "path/to/rec_cd_degs.csv"  # Rec vs CD DEG results (with Symbol column)
rec_cd_degs = pd.read_csv(rec_cd_degs_file)

rcParams['font.weight'] = 'bold'

s1 = set(nfd_cd_degs["Symbol"].dropna().astype(str))
s2 = set(rec_cd_degs["Symbol"].dropna().astype(str))

fig, ax = plt.subplots(figsize=(6, 6))
out = venn2([s1, s2],
            set_labels=('Reversible genes', 'Recovery-specific genes'),
            set_colors=("#a00000", "#eab77c"),
            alpha=0.80, ax=ax)

# Fonts (guard against None)
for t in out.set_labels:
    if t: t.set_fontsize(18)
for t in out.subset_labels:
    if t: t.set_fontsize(21)

# --- move the "Rec vs CD" set label (index 1) ---
lbl = out.set_labels[1]
x, y = lbl.get_position()
lbl.set_position((x + 0.60, y - 0.0))
lbl.set_ha('center')

plt.show()
