import json
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import font_manager, pyplot as plt

pd.set_option('display.max_rows', 500)

MAIN_DIR = Path(__file__).parent
RESULTS_FILE = MAIN_DIR.joinpath("jackpot-results.csv")
SIM_DIR = MAIN_DIR.joinpath("sim")
PARAMS_DIR = SIM_DIR.joinpath("params")
SAMPLES_DIR = SIM_DIR.joinpath("samples")
OUT_SEISMIC_DIR = MAIN_DIR.joinpath("out-seismic")
PRIMER_LENGTH = 20  # nt
MAX_AMPLICON_LENGTH = 292  # nt

# Configure matplotlib to write text as text in SVG files, not as paths
plt.rcParams['svg.fonttype'] = 'none'

# Add Helvetica Neue and Helvetica Neue Light to matplotlib's font list
font_dirs = []  # Add custom font directories if needed
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

# Set font family to Helvetica Neue for all text elements
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial', 'sans-serif']


def get_region_name(length: int):
    """ Get the region name for a given length. """
    if length < 0:
        raise ValueError(f"Invalid length: {length}")
    if length <= MAX_AMPLICON_LENGTH:
        return f"amp-{length}"
    return "full"


@cache
def get_refs():
    """ Map all references to their lengths, numbers of clusters, and
    replicates. """
    refs = dict()
    for length in [280]:
        for k in [1]:
            for rep in range(60):
                refs[f"rna-{length}-{k}-{rep}"] = (length, k, rep)
    return refs


@cache
def get_samples():
    """ Map all samples to their jackpotting iterations and total number of reads. """
    samples = dict()
    for jackpot in [0, 2, 4, 5, 6, 8, 10]:
        for reads in [200000]:
            samples[f"sample-jackpot-{jackpot}-{reads}"] = (jackpot, reads)
    return samples


def get_sample_results(sample: str, ref: str):
    print(f"Processing {sample} {ref}")
    length, true_k, rep = get_refs()[ref]
    jackpot, reads = get_samples()[sample]
    reg = get_region_name(length)
    cluster_dir = OUT_SEISMIC_DIR.joinpath(sample, "cluster", ref, reg)
    cluster_report_file = cluster_dir.joinpath("cluster-report.json")
    with open(cluster_report_file) as f:
        cluster_report_data = json.load(f)
    best_k = int(cluster_report_data["Best number of clusters"])
    summary_file = cluster_dir.joinpath("statistics", "summary.csv")
    summary = pd.read_csv(summary_file, header=[0], index_col=[0, 1])
    jackpot_quotient = summary.loc[(1, 0), "Jackpotting quotient"]
    bic_k1 = summary.loc[(1, 0), "Bayesian information criterion"]
    bic_k2 = summary.loc[(2, 0), "Bayesian information criterion"]
    return {
        "Length": length,
        "True K": true_k,
        "Rep": rep,
        "Jackpot Iterations": jackpot,
        "Reads": reads,
        "Best K": best_k,
        "Correct K": best_k == true_k,
        "Jackpot Quotient": jackpot_quotient,
        "BIC K1": bic_k1,
        "BIC K2": bic_k2,
        "BIC Diff": bic_k2 - bic_k1,
    }


def collect_results():
    return pd.DataFrame.from_records(
        get_sample_results(sample, ref)
        for sample in get_samples().keys()
        for ref in get_refs().keys()
    )


def graph_jackpot_quotient(results: pd.DataFrame):
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(4, 8))

    # Plot the distribution of jackpot quotients for each jackpot iteration using stripplot,
    # with each stripplot placed at its actual x value (numeric, not categorical)
    # This is achieved by passing the unique sorted x values as order, and using dodge=False.
    unique_jackpot_iters = np.sort(results["Jackpot Iterations"].unique())
    for x in unique_jackpot_iters:
        yvals = results.loc[results["Jackpot Iterations"] == x, "Jackpot Quotient"]
        # Add jitter to x for each point
        jittered_x = x + np.random.uniform(-0.15, 0.15, size=len(yvals))
        ax[0].scatter(
            jittered_x,
            yvals,
            color="#d55e00",
            alpha=0.5,
            s=10,
            edgecolor="none",
        )
    # Overlay the mean as a horizontal bar for each jackpot iteration at its actual x value
    means = results.groupby("Jackpot Iterations")["Jackpot Quotient"].mean()
    for jackpot, mean_val in means.items():
        ax[0].hlines(
            y=mean_val,
            xmin=jackpot - 0.2,
            xmax=jackpot + 0.2,
            color="black",
            linewidth=3,
            zorder=10,
        )
    ax[0].set_xlim(-0.5, 10.5)
    ax[0].set_xticks(np.linspace(0, 10, 6))
    ax[0].set_ylim(0.9, 2.6)
    ax[0].set_yticks(np.linspace(1.0, 2.5, 4))
    ax[0].set_ylabel("Jackpot Quotient")
    ax[0].grid(True, linestyle='-', alpha=0.5)

    # Plot the fraction of "Correct K" that is True vs. "Jackpot Iterations" on ax[1]
    correct_k_frac = results.groupby("Jackpot Iterations")["Correct K"].mean()
    ax[1].plot(
        correct_k_frac.index,
        correct_k_frac.values,
        marker="o",
        linestyle="-",
        color="#009e73",
        label="Fraction Correct K"
    )
    ax[1].set_xlim(-0.5, 10.5)
    ax[1].set_xticks(np.linspace(0, 10, 6))
    ax[1].set_xlabel("Rounds of Resampling")
    ax[1].set_ylabel("Accuracy of Clustering")
    ax[1].set_ylim(0, 1.05)
    ax[1].grid(True, linestyle='-', alpha=0.5)
    
    plt.savefig("jackpot-quotient.svg", dpi=300, bbox_inches='tight')
    plt.close()


def graph_bic_diff(results: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(4, 4))
    unique_jackpot_iters = np.sort(results["Jackpot Iterations"].unique())
    for x in unique_jackpot_iters:
        yvals = results.loc[results["Jackpot Iterations"] == x, "BIC Diff"]
        # Add jitter to x for each point
        jittered_x = x + np.random.uniform(-0.15, 0.15, size=len(yvals))
        ax.scatter(
            jittered_x,
            yvals,
            color="#0072b2",
            alpha=0.5,
            s=10,
            edgecolor="none",
        )
    # Overlay the mean as a horizontal bar for each jackpot iteration at its actual x value
    means = results.groupby("Jackpot Iterations")["BIC Diff"].mean()
    for jackpot, mean_val in means.items():
        ax.hlines(
            y=mean_val,
            xmin=jackpot - 0.2,
            xmax=jackpot + 0.2,
            color="black",
            linewidth=3,
            zorder=10,
        )
    ax.set_xlim(-0.5, 10.5)
    ax.set_xticks(np.linspace(0, 10, 6))
    #ax.set_ylim(0.9, 2.6)
    #ax.set_yticks(np.linspace(1.0, 2.5, 4))
    ax.set_ylabel("Difference in BIC")
    ax.grid(True, linestyle='-', alpha=0.5)
    plt.savefig("bic-diff.svg", dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    if not RESULTS_FILE.exists():
        results = collect_results()
        results.to_csv(RESULTS_FILE, index=False)
    else:
        results = pd.read_csv(RESULTS_FILE)
    graph_jackpot_quotient(results)
    graph_bic_diff(results)
