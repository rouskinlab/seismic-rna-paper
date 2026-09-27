from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import font_manager, patches, pyplot as plt

from seismicrna.core.rna import from_ct


# Configure matplotlib to write text as text in SVG files, not as paths
plt.rcParams['svg.fonttype'] = 'none'

# Add Helvetica Neue and Helvetica Neue Light to matplotlib's font list
font_dirs = []  # Add custom font directories if needed
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

# Set font family to Helvetica Neue for all text elements
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial', 'sans-serif']


# This reference truly has 1 cluster but is falsely clustered into 2
# when the observer bias is not corrected.
REF = "rna-140-1-15"
FOCUS = [83]
BLOCK = [80, 81, 82, 84, 85, 87]
FOCUS_COLOR = "#d55e00"
BLOCK_COLOR = "#e69f00"
OTHER_COLOR = "#56b4e9"
X_MAX = 0.105


def graph_mus():
    fig, axes = plt.subplots(nrows=2, figsize=(3.5, 3), gridspec_kw={'height_ratios': [3, 1]})

    # Mutation rates.
    ax = axes[0]
    mus_file = Path("sim", "params", REF, "full", "mus.csv")
    mus = pd.read_csv(mus_file, index_col=[0, 1], header=[0, 1]).squeeze()
    if not isinstance(mus, pd.Series):
        raise ValueError("mus is not a series")
    mus = mus.loc[np.isin(mus.index.get_level_values("Base"), ["A", "C"])]
    for (position, base), mu in mus.items():
        if position in FOCUS:
            color = FOCUS_COLOR
        elif position in BLOCK:
            color = BLOCK_COLOR
        else:
            color = OTHER_COLOR
        ax.bar(position, mu, color=color)
    ax.set_ylabel("Mutation Rate", fontsize=12)
    ax.set_ylim(bottom=0, top=X_MAX)
    ax.set_xticks(np.linspace(0, 140, 8))
    ax.set_xticklabels([])
    ax.set_xticklabels([])
    ax.set_yticks(np.linspace(0, 0.1, 6))
    ax.tick_params(axis='x', which='major', length=0)
    ax.tick_params(axis='y', which='major', labelsize=12)
    ax.grid(True, axis='both', clip_on=False, color='#E0E0E0')
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.set_axisbelow(True)

    # Structure.
    ax = axes[1]
    structure_file = Path("sim", "params", REF, "full", "simulated.ct")
    structures = list(from_ct(structure_file))
    if len(structures) != 1:
        raise ValueError("Expected 1 structure, got %d" % len(structures))
    structure = structures[0]
    for (p1, p2) in structure.pairs:
        # Calculate the center and dimensions of the ellipse
        center_x = (p1 + p2) / 2
        width = abs(p2 - p1)
        height = width / 140  # Adjust this factor to control the height of the arc
        # Create an arc patch
        arc = patches.Arc(
            xy=(center_x, 0),  # center point
            width=width,       # width of the arc
            height=height,     # height of the arc
            theta1=180,        # start angle
            theta2=0,          # end angle
            color='#808080',
            linewidth=1,
            zorder=0
        )
        ax.add_patch(arc)
    ax.set_xlim(left=0, right=140)
    ax.set_ylim(bottom=-0.225, top=0)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.set_xticks(np.linspace(0, 140, 8))
    ax.set_xlabel("Position", fontsize=12)
    ax.set_yticks([])
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_axisbelow(True)
    plt.subplots_adjust(hspace=0, left=0.2, right=0.95, top=0.95, bottom=0.2)
    plt.savefig("falsely_clustered_mus.svg")


def graph_clusters():
    clusters_file = Path("out-seismic",
                         "sample-biased-200000",
                         "graph_gap-0",
                         REF,
                         "amp-140",
                         "profile_clustered-2-x_m-ratio-q0.csv")
    clusters = pd.read_csv(clusters_file, skiprows=1, index_col=[0, 1], header=[0, 1])
    fig, ax = plt.subplots(figsize=(3, 3))
    # Plot each point with a different color
    ax.plot([0, X_MAX], [0, X_MAX], color="#e0e0e0", linewidth=1, zorder=0, label="y = x", clip_on=False)
    for (position, base), (x, y) in clusters.iterrows():
        if position in FOCUS:
            color = FOCUS_COLOR
        elif position in BLOCK:
            color = BLOCK_COLOR
        else:
            color = OTHER_COLOR
        ax.scatter(x, y, color=color, s=5, clip_on=False)        
        # Calculate the perpendicular line from point to y=x
        p = (x / X_MAX + X_MAX * y) / (1 + X_MAX**2) * X_MAX
        ax.plot([x, p], [y, p], color=color, linewidth=0.5, zorder=0)        
    ax.set_xlim(left=0, right=X_MAX)
    ax.set_ylim(bottom=0, top=1)
    ax.set_xticks(np.linspace(0, 0.1, 6))
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(fontsize=12)
    ax.set_xlabel("False Cluster 1", fontsize=12)
    ax.set_ylabel("False Cluster 2", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)
    ax.set_aspect(X_MAX)
    
    # Move spines slightly away from axes
    ax.spines['bottom'].set_position(('data', -0.02))
    ax.spines['left'].set_position(('data', -0.02 * X_MAX))
    
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.2)

    plt.savefig("falsely_clustered_scatter.svg")


if __name__ == "__main__":
    graph_mus()
    graph_clusters()
