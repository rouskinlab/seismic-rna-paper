from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import font_manager, pyplot as plt


# Configure matplotlib to write text as text in SVG files, not as paths
plt.rcParams['svg.fonttype'] = 'none'

# Add Helvetica Neue and Helvetica Neue Light to matplotlib's font list
font_dirs = []  # Add custom font directories if needed
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

# Set font family to Helvetica Neue for all text elements
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial', 'sans-serif']


def graph_mutdists():
    # Plot the observed and expected mutation distributions.
    fig, ax = plt.subplots()
    observed = dict()
    expected = dict()
    for rep in range(60):
        ref = f"rna-140-1-{rep}"
        file = Path("out-seismic",
                    "sample-biased-200000",
                    "graph_gap-0",
                    ref,
                    "amp-140",
                    "mutdist_filtered_m.csv")
        df = pd.read_csv(file, index_col=0)
        # Change the index from distance to gap.
        df.index = df.index - 1
        # Use reads with mutation gaps of 0 to 10.
        df = df.loc[0:10]
        observed[rep] = df["Mutated"]
        expected[rep] = df["Mutated-NULL"]
    observed = pd.DataFrame(observed)
    expected = pd.DataFrame(expected)
    style = {"marker": "o", "markersize": 3, "linewidth": 1}
    for rep in range(60):
        ax.plot(expected.index,
                expected[rep] / 1e3,
                color="#a6a6a6",
                alpha=0.1,
                linewidth=1)
        ax.plot(observed.index,
                observed[rep] / 1e3,
                color="#ff0000",
                alpha=0.05,
                linewidth=1)
    ax.plot(expected.index,
            expected.mean(axis=1) / 1e3,
            color="#a6a6a6",
            label="Theory",
            **style)
    ax.plot(observed.index,
            observed.mean(axis=1) / 1e3,
            color="#ff0000",
            label="Data",
            **style)
    ax.set_xlabel("Mutation Gap", fontsize=12)
    ax.set_ylabel("Reads (1000s)", fontsize=12)
    ax.set_xticks(np.arange(0, 11, 1))
    ax.set_ylim(bottom=0, top=10)
    ax.set_yticks(np.linspace(0, 10, 6))
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, clip_on=False, color='#E0E0E0')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()
    fig.set_size_inches(3, 2)
    plt.tight_layout()
    plt.savefig("mutdists.svg")


if __name__ == "__main__":
    graph_mutdists()
