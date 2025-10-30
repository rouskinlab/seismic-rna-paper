import json
import os
import re
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import font_manager, pyplot as plt


OS = "redhat"
DRACO_MERGED = "merged"
MAIN_DIR = Path(__file__).parent
RESULTS_FILE = MAIN_DIR.joinpath(f"results-{OS}-{DRACO_MERGED}.csv")
SIM_DIR = MAIN_DIR.joinpath("sim")
PARAMS_DIR = SIM_DIR.joinpath("params")
SAMPLES_DIR = SIM_DIR.joinpath("samples")
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
def get_refs():
    """ Map all references to their lengths, numbers of clusters, and
    replicates. """
    refs = dict()
    for d in PARAMS_DIR.iterdir():
        if d.is_dir():
            ref = d.name
            assert ref not in refs
            refs[ref] = ref_to_attrs(ref)
    return refs


def sample_to_attrs(sample: str):
    """ Map a sample to its bias and number of reads. """
    label, n_reads = re.match(r"^sample-([a-z0-9-]+)-(\d+)$", sample).groups()
    return label, int(n_reads)


@cache
def get_samples():
    """ Map all samples to their biases and numbers of reads. """
    samples = dict()
    for d in SAMPLES_DIR.iterdir():
        if d.is_dir():
            sample = d.name
            assert sample not in samples
            label, n_reads = sample_to_attrs(sample)
            if label == "biased":
                samples[sample] = n_reads
    return samples


def get_result_seismic(sample: str,
                       ref: str,
                       reg: str,
                       cpus: int):
    seismic_dir = MAIN_DIR.joinpath(f"out-seismic-{cpus}-cpu-{OS}")
    # Get the best number of clusters.
    cluster_report_file = seismic_dir.joinpath(
        sample, "cluster", ref, reg, "cluster-report.json"
    )
    with open(cluster_report_file) as f:
        cluster_report_data = json.load(f)
    best_k = int(cluster_report_data["Best number of clusters"])
    # Get the time taken for each step.
    time_preprocess = 0.0
    for step in ["align", "relate", "mask"]:
        time_file = seismic_dir.joinpath(f"time_seismic-{step}_{sample}_{ref}.txt")
        with open(time_file) as f:
            time_preprocess += float(f.read().strip())
    time_file = seismic_dir.joinpath(f"time_seismic-cluster_{sample}_{ref}.txt")
    with open(time_file) as f:
        time_cluster = float(f.read().strip())
    return best_k, time_preprocess, time_cluster


def get_result_dreem(sample: str,
                     ref: str,
                     end5: int,
                     end3: int):
    dreem_dir = MAIN_DIR.joinpath(f"out-dreem-{OS}")
    sample_ref = f"{sample}_{ref}"
    bv_dir = dreem_dir.joinpath(sample_ref, "BitVector_Plots")
    bv_log_file = bv_dir.joinpath(f"{sample_ref}_log.txt")
    em_dir = dreem_dir.joinpath(
        sample_ref, "EM_Clustering", f"{sample_ref}_{end5}_{end3}"
    )
    em_log_file = em_dir.joinpath("log.txt")
    # Get the best number of clusters.
    best_k = 0
    with open(em_log_file) as f:
        for line in f:
            if line.startswith("Predicted number of clusters"):
                best_k = int(line.split(":")[1].strip())
                break
    if best_k == 0:
        raise ValueError(f"Best number of clusters not found in {em_log_file}")
    # Get the time taken for each step.
    bv_log_file_mtime = os.path.getmtime(bv_log_file)
    em_log_file_mtime = os.path.getmtime(em_log_file)
    if em_log_file_mtime < bv_log_file_mtime:
        raise ValueError(f"{em_log_file_mtime} < {bv_log_file_mtime}")
    time_cluster = em_log_file_mtime - bv_log_file_mtime
    time_file = dreem_dir.joinpath(f"time_dreem_{sample}_{ref}.txt")
    with open(time_file) as f:
        time_total = float(f.read().strip())
    if time_total < time_cluster:
        raise ValueError(f"{time_total} < {time_cluster}")
    time_preprocess = time_total - time_cluster
    return best_k, time_preprocess, time_cluster


def get_result_draco(sample: str,
                     ref: str,
                     cpus: int):
    draco_dir = MAIN_DIR.joinpath(f"out-draco-{DRACO_MERGED}-{cpus}-cpu-{OS}")
    # Determine the windows, number of clusters, and proportions of
    # clusters for each window.
    stoichs_file = draco_dir.joinpath(
        sample, ref, "rf_json2rc", "stoichiometries.txt"
    )
    stoichs_df = pd.read_csv(stoichs_file, header=[0], sep="\t", dtype={ref: str})
    if stoichs_df.columns.to_list() != ["#Transcript",
                                        "Start",
                                        "End",
                                        "extStart",
                                        "extEnd",
                                        ref]:
        raise ValueError(f"Invalid columns in {stoichs_file}")
    if np.any(stoichs_df["#Transcript"] != ref):
        raise ValueError(f"Reference mismatch in {stoichs_file}")
    if np.any(stoichs_df["Start"] > stoichs_df["End"]):
        raise ValueError(f"Start > End in {stoichs_file}")
    if np.any(stoichs_df["extStart"] > stoichs_df["extEnd"]):
        raise ValueError(f"extStart > extEnd in {stoichs_file}")
    # Get the number of clusters and their stoichiometries for each window.
    stoichs = {(end5, end3): list(map(float, window_stoichs.split(";")))
               for end5, end3, window_stoichs in zip(stoichs_df["Start"],
                                                     stoichs_df["End"],
                                                     stoichs_df[ref])}
    ks = {window: len(window_stoichs)
          for window, window_stoichs in stoichs.items()}
    if not ks:
        raise ValueError(stoichs_file)
    assert min(ks.values()) >= 1
    # To get a single number of clusters for DRACO, use the maximum
    # number of clusters among all windows.
    best_k = max(ks.values())
    # Get the time taken for each step.
    time_preprocess = 0.0
    for step in ["map", "count"]:
        time_file = draco_dir.joinpath(f"time_draco-{step}_{sample}_{ref}.txt")
        with open(time_file) as f:
            time_preprocess += float(f.read().strip())
    time_cluster = 0.0
    for step in ["cluster", "json2rc", "norm"]:
        time_file = draco_dir.joinpath(f"time_draco-{step}_{sample}_{ref}.txt")
        with open(time_file) as f:
            time_cluster += float(f.read().strip())
    return best_k, time_preprocess, time_cluster


def get_result_dance(sample: str,
                     ref: str,
                     cpus: int):
    dance_dir = MAIN_DIR.joinpath(f"out-dance-{cpus}-cpu-{OS}")
    reactivities_file = dance_dir.joinpath(
            "dancemapper", f"{sample}_{ref}-reactivities.txt"
    )
    with open(reactivities_file) as f:
        match = re.match("^([0-9]+) components; BIC=[0-9.-]+$",
                            f.readline().rstrip())
        if not match:
            raise ValueError(reactivities_file)
        best_k = int(match.groups()[0])
    time_file = dance_dir.joinpath(f"time_dance-sm_{sample}_{ref}.txt")
    with open(time_file) as f:
        time_preprocess = float(f.read().strip())
    time_file = dance_dir.joinpath(f"time_dance-dm_{sample}_{ref}.txt")
    with open(time_file) as f:
        time_cluster = float(f.read().strip())
    return best_k, time_preprocess, time_cluster


def get_sample_results(ref: str,
                       length: int,
                       true_k: int,
                       rep: int,
                       sample: str,
                       n_reads_total: int,
                       cpus: int,
                       only_true_k: bool = True):
    results = list()

    if length != 280 or true_k not in [2, 4] or n_reads_total != 200000:
        return results

    print(f"Processing {sample} {ref} {length} {true_k} {rep}")

    # SEISMIC
    print("  SEISMIC")
    best_k, time_preprocess, time_cluster = get_result_seismic(
        sample, ref, get_region_name(length), cpus
    )
    if best_k == true_k or not only_true_k:
        results.append({"Length": length,
                        "True K": true_k,
                        "Total Reads": n_reads_total,
                        "Replicate": rep,
                        "Method": "SEISMIC",
                        "CPUs": cpus,
                        "Best K": best_k,
                        "Preprocessing Time (s)": time_preprocess,
                        "Clustering Time (s)": time_cluster})

    # DANCE
    print("  DANCE")
    best_k, time_preprocess, time_cluster = get_result_dance(
        sample, ref, cpus
    )
    if best_k == true_k or not only_true_k:
        results.append({"Length": length,
                        "True K": true_k,
                        "Total Reads": n_reads_total,
                        "Replicate": rep,
                        "Method": "DANCE",
                        "CPUs": cpus,
                        "Best K": best_k,
                        "Preprocessing Time (s)": time_preprocess,
                        "Clustering Time (s)": time_cluster})

    # DRACO
    print("  DRACO")
    best_k, time_preprocess, time_cluster = get_result_draco(
            sample, ref, cpus
    )
    if best_k == true_k or not only_true_k:
        results.append({"Length": length,
                        "True K": true_k,
                        "Total Reads": n_reads_total,
                        "Replicate": rep,
                        "Method": "DRACO",
                        "CPUs": cpus,
                        "Best K": best_k,
                        "Preprocessing Time (s)": time_preprocess,
                        "Clustering Time (s)": time_cluster})

    # DREEM
    if length <= MAX_AMPLICON_LENGTH and cpus == 1:
        print("  DREEM")
        end5 = PRIMER_LENGTH + 1
        end3 = length - PRIMER_LENGTH
        best_k, time_preprocess, time_cluster = get_result_dreem(
            sample, ref, end5, end3
        )
        if best_k == true_k or not only_true_k:
            results.append({"Length": length,
                            "True K": true_k,
                            "Total Reads": n_reads_total,
                            "Replicate": rep,
                            "Method": "DREEM",
                            "CPUs": cpus,
                            "Best K": best_k,
                            "Preprocessing Time (s)": time_preprocess,
                            "Clustering Time (s)": time_cluster})
    
    return results


def get_ref_results(ref: str,
                    length: int,
                    true_k: int,
                    rep: int,
                    cpus: int):
    results = list()
    for sample, n_reads_total in get_samples().items():
        results.extend(get_sample_results(ref, length, true_k, rep,
                                          sample, n_reads_total, cpus))
    return results


def collect_results():
    results = list()
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        print(f"Opened {executor}")
        args = [(ref, length, true_k, rep, cpus)
                for ref, (length, true_k, rep) in get_refs().items()
                for cpus in [1]]
        for result in executor.map(get_ref_results, *zip(*args)):
            results.extend(result)
        print(f"Closed {executor}")
    return pd.DataFrame.from_records(results)


METHOD_STYLES = {"SEISMIC": {"color": "#0072b2",
                             "marker": "o",
                             "markersize": 3.0 * 1.125**0,
                             "linewidth": 1.0 * 1.125**0,
                             "zorder": 13},
                 "DANCE": {"color": "#009e73",
                           "marker": "o",
                           "markersize": 3.0 * 1.125**0,
                           "linewidth": 1.0 * 1.125**0,
                           "zorder": 12},
                 "DRACO": {"color": "#e69f00",
                           "marker": "o",
                           "markersize": 3.0 * 1.125**0,
                           "linewidth": 1.0 * 1.125**0,
                           "zorder": 11},
                 "DREEM": {"color": "#cc79a7",
                           "marker": "o",
                           "markersize": 3.0 * 1.125**0,
                           "linewidth": 1.0 * 1.125**0,
                           "zorder": 10}}


def graph_time(results: pd.DataFrame, time_col: str, fname: str):
    """ Graph the time for each method. """
    true_ks = np.sort(results["True K"].unique())
    #assert np.array_equal(true_ks, np.arange(1, true_ks.size + 1))

    fig, axs = plt.subplots(nrows=true_ks.size,
                            squeeze=True,                            
                            sharex=True,
                            sharey=True,
                            gridspec_kw={'hspace': 0.5},
                            figsize=(4, 4))
    
    for i, true_k in enumerate(true_ks):
        ax = axs[i]
        ax_data = results[(results["CPUs"] == 1) & (results["True K"] == true_k)]
        # Switch x and y axes: plot methods on y-axis, time on x-axis
        for i, (method, style) in enumerate(METHOD_STYLES.items()):
            method_data = ax_data[ax_data["Method"] == method][time_col]
            if method_data.size == 0:
                continue
            color = style['color']
            # Position each violin at the y-coordinate corresponding to the method index
            sns.stripplot(y=np.full(method_data.size, i), x=method_data, ax=ax, color=color, alpha=0.5, size=3.0, orient='h')
            # Add a vertical line at the mean
            mean = method_data.mean()
            print(f"{true_k} Clusters, {method}: {mean} minutes")
            ax.vlines(mean, i-0.2, i+0.2, color="black", linewidth=1)
        
        # Set y-ticks and labels manually
        ax.set_yticks(range(len(METHOD_STYLES)))
        ax.set_yticklabels(METHOD_STYLES.keys())
        plt.setp(ax.get_yticklabels(), rotation=0, va='center')
        ax.set_xlim(left=5, right=2000)
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(plt.LogLocator(base=10, numticks=4))
        ax.xaxis.set_major_formatter(plt.ScalarFormatter())
        ax.set_xticks([10, 100, 1000])

        ax.set_xlabel(time_col)
        ax.set_ylabel("")

        ax.grid(True, clip_on=False, color='#E0E0E0', axis='x')
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ax.get_legend():
            ax.get_legend().remove()
        
        ax.set_title(f"True Clusters: {true_k}",
                     fontsize=12,
                     pad=10)

    # hspace in plt.subplots_adjust is in fraction of the figure height, not absolute units.
    # For more vertical space, increase hspace to a value like 0.6 (default is 0.2).
    plt.subplots_adjust(left=0.25, right=0.95, top=0.9, bottom=0.15)

    plt.savefig(fname)
    plt.close()


def graph_time_total(results: pd.DataFrame, fname: str):
    graph_time(results, "Total Runtime (minutes)", fname)


if __name__ == "__main__":
    format = "svg"
    if not RESULTS_FILE.exists():
        results = collect_results()
        results["Preprocessing Time (minutes)"] = results["Preprocessing Time (s)"] / 60.0
        results["Clustering Time (minutes)"] = results["Clustering Time (s)"] / 60.0
        results["Total Runtime (minutes)"] = results["Preprocessing Time (minutes)"] + results["Clustering Time (minutes)"]
        results.to_csv(RESULTS_FILE, index=False)
    else:
        results = pd.read_csv(RESULTS_FILE)
    graph_time_total(results, f"total-time-{OS}-{DRACO_MERGED}.{format}")
    
