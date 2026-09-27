import json
import os
import re
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pandas as pd
from matplotlib import font_manager, pyplot as plt


MAIN_DIR = Path(__file__).parent
RESULTS_FILE = MAIN_DIR.joinpath("results_biased-vs-nobias.csv")
SIM_DIR = MAIN_DIR.joinpath("sim")
PARAMS_DIR = SIM_DIR.joinpath("params")
SAMPLES_DIR = SIM_DIR.joinpath("samples")
OUT_SEISMIC_DIR = MAIN_DIR.joinpath("out-seismic")
BIAS_TYPES = {
    True: "biased",
    #False: "nobias"
}
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


def ref_to_attrs(ref: str):
    """ Map a reference to its length, number of clusters, and
    replicate. """
    length, k, rep = map(int,
                         re.match(r"^rna-(\d+)-(\d+)-(\d+)$", ref).groups())
    return length, k, rep


def get_region_name(length: int):
    """ Get the region name for a given length. """
    if length < 0:
        raise ValueError(f"Invalid length: {length}")
    if length <= MAX_AMPLICON_LENGTH:
        return f"amp-{length}"
    return "full"


@cache
def get_real_params(ref: str):
    """ Get the real mutation rates for a given reference length, number
    of clusters (k), and replicate. """
    ref_reg_dir = PARAMS_DIR.joinpath(ref, "full")
    mus = pd.read_csv(ref_reg_dir.joinpath("mus.csv"),
                      index_col=[0, 1],
                      header=[0, 1])
    # Remove the bases from the index.
    mus.index = mus.index.get_level_values("Position")
    # Convert column labels from str to int.
    mus.columns = pd.MultiIndex.from_tuples([tuple(map(int, col))
                                             for col in mus.columns],
                                            names=mus.columns.names)
    pis = pd.read_csv(ref_reg_dir.joinpath("simulated.clusts.csv"),
                      index_col=[0, 1]).squeeze("columns")
    assert isinstance(mus, pd.DataFrame)
    assert isinstance(pis, pd.Series)
    assert mus.columns.equals(pis.index)
    return mus, pis


def calc_rmse(x: np.ndarray, y: np.ndarray):
    """ Calculate the RMSE of two arrays. """
    return np.sqrt(np.mean(np.square(x - y)))


def assign_clusterings(mus1: np.ndarray, mus2: np.ndarray):
    """ Optimally assign clusters from two groups to each other. """
    n1, k1 = mus1.shape
    n2, k2 = mus2.shape
    assert n1 == n2
    assert k1 >= 1
    assert k2 >= 1
    if n1 >= 1:
        # Match the clusters using linear_sum_assignment.
        costs = np.array([[calc_rmse(mus1[:, cluster1],
                                     mus2[:, cluster2])
                           for cluster2 in range(k2)]
                          for cluster1 in range(k1)]).reshape((k1, k2))
        from scipy.optimize import linear_sum_assignment
        rows, cols = linear_sum_assignment(costs)
    else:
        # If n1 == 0, then the costs matrix will contain NaN, which will
        # cause linear_sum_assignment to raise an error.
        rows = np.arange(k1)
        cols = np.arange(k2)
    return rows, cols


def calc_params_accuracy(mus: pd.DataFrame,
                         pis: pd.Series,
                         real_mus: pd.DataFrame,
                         real_pis: pd.Series):
    """ Calculate accuracies of the mutation rates and proportions. """
    # Make sure the indexes are aligned.
    if isinstance(mus.index, pd.MultiIndex):
        mus.index = mus.index.get_level_values("Position")
    mus = mus.reindex(real_mus.index, fill_value=np.nan)
    assert mus.index.equals(real_mus.index)
    if pis.size < real_pis.size:
        pis = pis.reindex(real_pis.index, fill_value=0.0)
    elif pis.size > real_pis.size:
        real_pis = real_pis.reindex(pis.index, fill_value=0.0)
    elif not pis.index.equals(real_pis.index):
        raise ValueError("Indexes of pis and real_pis do not match")
    assert pis.index.equals(real_pis.index)
    # Keep only the values.
    mus = mus.values
    real_mus = real_mus.values
    pis = pis.values
    real_pis = real_pis.values
    # Remove any NaN mutation rates.
    use_mus = ~np.any(np.isnan(mus), axis=1)
    mus = mus[use_mus]
    real_mus = real_mus[use_mus]
    # Match the clusters from the results and parameters.
    rows, cols = assign_clusterings(mus, real_mus)
    mus = mus[:, rows]
    real_mus = real_mus[:, cols]
    pis = pis[rows]
    real_pis = real_pis[cols]
    assert mus.shape == real_mus.shape
    assert pis.shape == real_pis.shape
    # Calculate the RMSE of the mutation rates.
    mus_rmse = calc_rmse(mus, real_mus)
    # Calculate the RMSE of the proportions.
    pis_rmse = calc_rmse(pis, real_pis)
    return mus_rmse, pis_rmse


def get_result_seismic(length: int,
                       true_k: int,
                       rep: int,
                       n_reads_total: int,
                       is_biased: bool,
                       is_corrected: bool,
                       real_mus: pd.DataFrame,
                       real_pis: pd.Series):
    # Get the best number of clusters.
    sample = f"sample-{BIAS_TYPES[is_biased]}-{n_reads_total}"
    ref = f"rna-{length}-{true_k}-{rep}"
    reg = get_region_name(length)
    if is_corrected:
        cluster_dir = "cluster"
        graph_dir = "graph"
    else:
        cluster_dir = "cluster_gap-0"
        graph_dir = "graph_gap-0"
    cluster_report_file = OUT_SEISMIC_DIR.joinpath(
        sample, cluster_dir, ref, reg, "cluster-report.json"
    )
    with open(cluster_report_file) as f:
        cluster_report_data = json.load(f)
    best_k = int(cluster_report_data["Best number of clusters"])
    # Get accuracies of the mutation rates and proportions.
    mus_file = OUT_SEISMIC_DIR.joinpath(
        sample, graph_dir, ref, reg,
        f"profile_clustered-{best_k}-x_m-ratio-q0.csv"
    )
    mus = pd.read_csv(mus_file,
                      index_col=[0, 1],
                      header=[0, 1, 2])
    pis_file = OUT_SEISMIC_DIR.joinpath(
        sample, graph_dir, ref, reg,
        "abundance_clustered.csv"
    )
    pis = pd.read_csv(pis_file,
                      index_col=[0],
                      header=[0])
    pis.columns = pis.columns.astype(int)
    pis = pis.T[best_k]
    mus_rmse, pis_rmse = calc_params_accuracy(mus, pis, real_mus, real_pis)
    return best_k, mus_rmse, pis_rmse


def get_sample_results(length: int,
                       true_k: int,
                       rep: int,
                       n_reads_total: int,
                       is_biased: bool,
                       is_corrected: bool,
                       real_mus: pd.DataFrame,
                       real_pis: pd.Series):
    print(f"Processing {length} {true_k} {n_reads_total} "
          f"{is_biased} {is_corrected} {rep}")
    # Get the best number of clusters.
    best_k, mus_rmse, pis_rmse = get_result_seismic(
        length, true_k, rep, n_reads_total, is_biased, is_corrected,
        real_mus, real_pis
    )
    return {"Length": length,
            "True K": true_k,
            "Total Reads": n_reads_total,
            "Biased": is_biased,
            "Corrected": is_corrected,
            "Replcate": rep,
            "Result K": best_k,
            "Mus RMSE": mus_rmse,
            "Pis RMSE": pis_rmse}


def get_ref_results(length: int,
                    true_k: int,
                    rep: int):
    ref = f"rna-{length}-{true_k}-{rep}"
    real_mus, real_pis = get_real_params(ref)
    real_pis = real_pis.loc[true_k]
    results = list()
    for n_reads_total in [200000]:
        for is_biased in BIAS_TYPES:
            for is_corrected in [True, False]:
                results.append(get_sample_results(
                    length, true_k, rep,
                    n_reads_total, is_biased, is_corrected,
                    real_mus, real_pis
                ))
    return results


def collect_results():
    results = list()
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        print(f"Opened {executor}")
        args = [(length, 1, rep)
                for length in [140]
                for rep in range(60)]
        for result in executor.map(get_ref_results, *zip(*args)):
            results.extend(result)
        print(f"Closed {executor}")
    return pd.DataFrame.from_records(results)


def graph_fraction_each_k(results: pd.DataFrame, fname: str):
    """ Graph the fraction of each K for each method. """
    max_k = 4

    fig, axs = plt.subplots(nrows=1,
                            ncols=2,
                            squeeze=False,
                            sharex=True,
                            sharey=True,
                            gridspec_kw={'hspace': 0.2,
                                         'wspace': 0.2},
                            figsize=(3, 2))
    row = 0
    for col, is_corrected in enumerate([False, True]):
        ax = axs[row, col]
        ax_data = results[results["Corrected"] == is_corrected]

        for k in range(1, max_k + 1):
            if k < max_k:
                fraction_k = np.mean(ax_data["Result K"] == k)
            else:
                fraction_k = np.mean(ax_data["Result K"] >= max_k)
            color = "#009e73" if k == 1 else "#d55e00"
            ax.bar(k, fraction_k, clip_on=False, color=color)
        
        ax.set_xlim(0.5, max_k + 0.5)
        ax.set_ylim(0, 1)
        ax.set_xticks(np.arange(1, max_k + 1))
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.set_axisbelow(True)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(
            lambda x, _: f"{x}+" if x >= max_k else str(x)
        ))
        ax.tick_params(axis='both', which='major', labelsize=12)
        
        x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        ax.set_aspect(1.1 * x_range / y_range)

        ax.grid(True, axis='y', clip_on=False, color='#E0E0E0')
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ax.get_legend():
            ax.get_legend().remove()
        
        if row == 0:
            ax.set_title("On" if is_corrected else "Off (Naive)",
                        fontsize=12,
                        pad=10)
    
    fig.text(0.45, 0.97, "SEISMIC-RNA Bias Correction Algorithm", ha='center', fontsize=12)
    fig.text(0.45, 0.03, "Number of Clusters Detected", ha='center', fontsize=12)
    fig.text(0.02, 0.50, "Proportion of\nSimulations", 
             va='center', rotation='vertical', fontsize=12)

    plt.subplots_adjust(left=0.175, right=0.925, top=0.9, bottom=0.15)

    plt.savefig(fname)
    plt.close()


if __name__ == "__main__":
    format = "svg"
    if not RESULTS_FILE.exists():
        results = collect_results()
        results.to_csv(RESULTS_FILE, index=False)
    else:
        results = pd.read_csv(RESULTS_FILE)
    graph_fraction_each_k(results[results["Biased"] == True], f"biased-vs-nobias_each-k.{format}")
    
