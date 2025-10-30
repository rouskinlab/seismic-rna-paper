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


def graph_mutdist(file: str):
    df = pd.read_csv(file, index_col=0)
    # Change the index from distance to gap.
    df.index = df.index - 1
    # Use reads with mutation gaps of 0 to 10.
    df = df.loc[0:10]
    observed = df["Mutated"]
    expected = df["Mutated-NULL"]
    # Plot the observed and expected mutation distributions.
    fig, ax = plt.subplots()
    style = {"marker": "o", "markersize": 3, "linewidth": 1}
    ax.plot(expected.index, expected / 1e6, color="#a6a6a6", label="Theory", **style)
    ax.plot(observed.index, observed / 1e6, color="#ff0000", label="Data", **style)
    ax.set_xlabel("Mutation Gap", fontsize=12)
    ax.set_ylabel("Reads (millions)", fontsize=12)
    ax.set_xticks(np.arange(0, 11, 1))
    ax.set_ylim(bottom=0, top=2)
    ax.set_yticks(np.linspace(0, 2, 5))
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, clip_on=False, color='#E0E0E0')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()
    fig.set_size_inches(3, 2)
    plt.tight_layout()
    plt.savefig("mutdist.svg")


if __name__ == "__main__":
    file = "out-morandi-2021/pooled/graph_mutdist/sars-cov-2/full/mutdist_filtered_m.csv"
    graph_mutdist(file)
