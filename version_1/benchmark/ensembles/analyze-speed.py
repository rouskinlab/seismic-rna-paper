import json
import os
import re
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from pathlib import Path
import time
from xml.etree import ElementTree

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import font_manager, pyplot as plt


MAIN_DIR = Path(__file__).parent
RESULTS_FILE = MAIN_DIR.joinpath(f"results-speed.csv")
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
    match = re.match(r"^long-rna-(\d+)$", ref)
    if not match:
        raise ValueError(ref)
    rep, = map(int, match.groups())
    return rep


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
            try:
                refs[ref] = ref_to_attrs(ref)
            except ValueError:
                pass
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


def get_result_seismic(sample: str, ref: str):
    seismic_dir = MAIN_DIR.joinpath(f"out-seismic")
    # Get the time taken for each step.
    time = 0.0
    for step in ["align", "relate", "ensembles"]:
        time_file = seismic_dir.joinpath(f"time_seismic-{step}_{sample}_{ref}.txt")
        with open(time_file) as f:
            time += float(f.read().strip())
    return time


def get_result_draco(sample: str, ref: str):
    draco_dir = MAIN_DIR.joinpath(f"out-draco")
    # Get the time taken for each step.
    time = 0.0
    for step in ["map", "count", "cluster-100", "json2rc-100", "norm-100"]:
        time_file = draco_dir.joinpath(f"time_draco-{step}_{sample}_{ref}.txt")
        with open(time_file) as f:
            time += float(f.read().strip())
    return time


def get_sample_results(ref: str,
                       rep: int,
                       sample: str,
                       n_reads_total: int):
    results = list()

    print(f"Processing {sample} {ref} {rep}")

    # SEISMIC
    print("  SEISMIC")
    time_seismic = get_result_seismic(sample, ref)
    results.append({"Total Reads": n_reads_total,
                    "Replicate": rep,
                    "Method": "SEISMIC",
                    "Time (s)": time_seismic})

    # DRACO
    print("  DRACO")
    time_draco = get_result_draco(sample, ref)
    results.append({"Total Reads": n_reads_total,
                    "Replicate": rep,
                    "Method": "DRACO",
                    "Time (s)": time_draco})
    
    return results


def get_ref_results(ref: str,
                    rep: int):
    results = list()
    for sample, n_reads_total in get_samples().items():
        results.extend(get_sample_results(ref, rep,
                                          sample, n_reads_total))
    return results


def collect_results():
    results = list()
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        print(f"Opened {executor}")
        args = [(ref, rep)
                for ref, rep in get_refs().items()]
        for result in executor.map(get_ref_results, *zip(*args)):
            results.extend(result)
        print(f"Closed {executor}")
    return pd.DataFrame.from_records(results)


METHOD_STYLES = {"SEISMIC": {"color": "#0072b2",
                             "marker": "o",
                             "markersize": 3.0 * 1.125**0,
                             "linewidth": 1.0 * 1.125**0,
                             "zorder": 13},
                 "DRACO": {"color": "#e69f00",
                           "marker": "o",
                           "markersize": 3.0 * 1.125**0,
                           "linewidth": 1.0 * 1.125**0,
                           "zorder": 11}}


def graph_time(results: pd.DataFrame, time_col: str, fname: str):
    """ Graph the time for each method. """
    fig, ax = plt.subplots(squeeze=True,                            
                           sharex=True,
                           sharey=True,
                           gridspec_kw={'hspace': 0.5},
                           figsize=(4, 1.25))

    # Switch x and y axes: plot methods on y-axis, time on x-axis
    for i, (method, style) in enumerate(METHOD_STYLES.items()):
        method_data = results[results["Method"] == method][time_col]
        if method_data.size == 0:
            continue
        color = style['color']
        # Position each violin at the y-coordinate corresponding to the method index
        sns.stripplot(y=np.full(method_data.size, i), x=method_data, ax=ax, color=color, alpha=0.5, size=3.0, orient='h')
        # Add a vertical line at the mean
        mean = method_data.mean()
        print(f"{method}: {mean} minutes")
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

    # hspace in plt.subplots_adjust is in fraction of the figure height, not absolute units.
    # For more vertical space, increase hspace to a value like 0.6 (default is 0.2).
    plt.subplots_adjust(left=0.25, right=0.95, top=0.9, bottom=0.4)

    plt.savefig(fname)
    plt.close()


def graph_time_total(results: pd.DataFrame, fname: str):
    graph_time(results, "Total Runtime (minutes)", fname)


if __name__ == "__main__":
    format = "svg"
    if not RESULTS_FILE.exists():
        results = collect_results()
        results["Total Runtime (minutes)"] = results["Time (s)"] / 60
        results.to_csv(RESULTS_FILE, index=False)
    else:
        results = pd.read_csv(RESULTS_FILE)
    graph_time_total(results, f"total-time.{format}")
    
