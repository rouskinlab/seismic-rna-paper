"""Plot distributions of informative bases, mutations, and mutation rates per read/position.

Reads the SEISMIC-RNA graph outputs for the morandi-2021 and lan-2022 samples and
plots, as a 3-row x 2-column grid of histograms:
  1. Informative bases per read   (histread_filtered_n-count.csv)
  2. Mutations per read           (histread_filtered_m-count.csv)
  3. Mutation rate per position   (profile_filtered_m-ratio-q0.csv, "Mutated" column)

Each row shares its x- and y-axes across the two samples; rows are independent of
one another. Run from the directory that contains the morandi-2021/ and lan-2022/
sample directories (i.e. so that ./<sample>/out-seismic/pooled/graph/sars-cov-2/full/
resolves), or pass --base-dir.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

# Configure matplotlib to write text as text (not paths) in SVG output, and to
# match the house style used elsewhere in this repo (see parameterize/graph_mutdist.py).
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "sans-serif"]

GRAPH_SUBPATH = Path("out-seismic/pooled/graph/sars-cov-2/full")

# Sample directory name -> display label (column order is fixed).
SAMPLES = {
    "morandi-2021": "Morandi et al. (2021)",
    "lan-2022": "Lan et al. (2022)",
}

# One color per sample, held constant across all three rows. Orange and
# reddish purple from Wong, B. "Points of view: Color blindness." Nat Methods
# 8, 441 (2011).
COLORS = {
    "morandi-2021": "#E69F00",
    "lan-2022": "#CC79A7",
}

GRID_COLOR = "#e0e0e0"
TEXT_COLOR = "#333333"

MM_PER_IN = 25.4
FIG_WIDTH_IN = 180 / MM_PER_IN
FIG_HEIGHT_IN = 210 / MM_PER_IN


def load_read_hist(sample_dir: Path, rel_code: str, rel_name: str) -> pd.Series:
    """Load a pre-computed per-read histogram (Count -> number of reads)."""
    df = pd.read_csv(sample_dir / f"histread_filtered_{rel_code}-count.csv")
    return df.set_index("Count")[rel_name]


def load_mut_rates(sample_dir: Path) -> pd.Series:
    """Load per-position mutation rates (fraction of informative reads mutated)."""
    df = pd.read_csv(sample_dir / "profile_filtered_m-ratio-q0.csv")
    return df["Mutated"].dropna()


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=1)
    ax.tick_params(length=0)


def annotate_mean(ax, mean: float, fmt: str):
    ax.text(
        0.95, 0.92, f"Mean = {mean:{fmt}}",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9, color=TEXT_COLOR,
    )


def plot_count_hist(ax, counts: pd.Series, color: str, mean_fmt: str):
    # Counts are already tabulated per integer value (bin width 1), so density
    # is just each count divided by the total (no re-binning, no edge effects).
    total = counts.sum()
    density = counts.to_numpy(dtype=float) / total
    ax.bar(counts.index, density, width=1.0, color=color)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    mean = float((counts.index.to_numpy() * counts.to_numpy()).sum() / total)
    annotate_mean(ax, mean, mean_fmt)


def plot_rate_hist(ax, rates: pd.Series, color: str, rate_max: float, rate_step: float):
    bins = np.arange(0.0, rate_max + rate_step, rate_step)
    values = rates.to_numpy()
    weights = np.full_like(values, 1.0 / values.size)
    ax.hist(values, bins=bins, weights=weights, color=color)
    ax.set_xlim(0.0, rate_max)
    annotate_mean(ax, float(rates.mean()), ".3f")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("."),
        help="Directory containing the sample subdirectories (default: current directory)",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("distributions.svg"),
        help="Output image path (default: distributions.svg)",
    )
    parser.add_argument(
        "--rate-max",
        type=float,
        default=0.2,
        help="Upper limit of the mutation-rate axis (default: 0.2)",
    )
    parser.add_argument(
        "--rate-step",
        type=float,
        default=0.005,
        help="Bin width for the mutation-rate histogram (default: 0.005)",
    )
    args = parser.parse_args()

    sample_dirs = {}
    for sample in SAMPLES:
        sample_dir = args.base_dir / sample / GRAPH_SUBPATH
        if not sample_dir.is_dir():
            raise FileNotFoundError(f"Missing directory for sample {sample!r}: {sample_dir}")
        sample_dirs[sample] = sample_dir

    fig, axes = plt.subplots(
        nrows=3, ncols=len(SAMPLES),
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        sharex="row", sharey="row",
    )

    for col, sample in enumerate(SAMPLES):
        sample_dir = sample_dirs[sample]
        color = COLORS[sample]

        plot_count_hist(axes[0, col], load_read_hist(sample_dir, "n", "Informative"), color, ".1f")
        plot_count_hist(axes[1, col], load_read_hist(sample_dir, "m", "Mutated"), color, ".2f")
        plot_rate_hist(axes[2, col], load_mut_rates(sample_dir), color, args.rate_max, args.rate_step)

        axes[0, col].set_title(SAMPLES[sample], fontsize=12, fontweight="bold", pad=10)

    row_labels = [
        ("Informative bases per Read", "Fraction of Reads"),
        ("Mutations per read", "Fraction of Reads"),
        ("Mutation rate", "Fraction of Positions"),
    ]
    for row, (xlabel, ylabel) in enumerate(row_labels):
        for col in range(len(SAMPLES)):
            ax = axes[row, col]
            style_axis(ax)
            ax.set_xlabel(xlabel, fontsize=10)
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=10)

    fig.tight_layout()
    fig.savefig(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
