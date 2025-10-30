import json
import os
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import cache, reduce
from pathlib import Path
from xml.etree import ElementTree
from typing import Callable

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import font_manager, pyplot as plt
from scipy import stats

pd.set_option('display.max_rows', 500)

MAIN_DIR = Path(__file__).parent
RESULTS_FILE = MAIN_DIR.joinpath("results.csv")
SIM_DIR = MAIN_DIR.joinpath("sim")
PARAMS_DIR = SIM_DIR.joinpath("params")
SAMPLES_DIR = SIM_DIR.joinpath("samples")
OUT_SEISMIC_DIR = MAIN_DIR.joinpath("out-seismic")
OUT_DREEM_DIR = MAIN_DIR.joinpath("out-dreem")
OUT_DRACO_DIR = MAIN_DIR.joinpath("out-draco")
OUT_DANCE_DIR = MAIN_DIR.joinpath("out-dance")
BIAS_TYPES = {True: "biased", False: "nobias"}
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
plt.rcParams['font.size'] = 12


def mm_to_in(mm: float):
    """ Convert millimeters to inches. """
    return mm / 25.4


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
    match = re.match(r"^sample-([a-z]{6})-(\d+)$", sample)
    if not match:
        raise ValueError(sample)
    bias, n_reads = match.groups()
    for is_biased, bias_type in BIAS_TYPES.items():
        if bias == bias_type:
            break
    else:
        raise ValueError(f"Invalid bias: {bias}")
    if not is_biased:
        raise ValueError(bias)
    return is_biased, int(n_reads)


@cache
def get_samples():
    """ Map all samples to their biases and numbers of reads. """
    samples = dict()
    for d in SAMPLES_DIR.iterdir():
        if d.is_dir():
            sample = d.name
            assert sample not in samples
            try:
                samples[sample] = sample_to_attrs(sample)
            except ValueError:
                pass
    return samples


def get_data_series(sample: str, ref: str, full_region: bool,
                    csv_filename: str, colname: str, pos_base_index: bool):
    length, k, rep = ref_to_attrs(ref)
    region = "full" if full_region else get_region_name(length)
    csv = OUT_SEISMIC_DIR.joinpath(sample, "graph", ref, region, csv_filename)
    print(f"Reading {csv}")
    if pos_base_index:
        ni = 2
    else:
        ni = 1
    df = pd.read_csv(csv, index_col=list(range(ni)), header=[0])
    if pos_base_index:
        df.index = df.index.get_level_values("Position")
    return df[colname]


def graph_avg_of_refs(graph_filename: str,
                      fmt: str, *,
                      xlabel: str,
                      ylabel: str,
                      xticks: np.ndarray | None = None,
                      bins: np.ndarray | None = None,
                      **kwargs):
    accum = defaultdict(list)
    for ref, (length, true_k, rep) in get_refs().items():
        series = get_data_series(ref=ref, **kwargs)
        if bins is not None:
            if not series.max() <= bins.max():
                raise ValueError(f"Series max is {series.max()}, "
                                 f"greater than max bin {bins.max()}")
            values, edges = np.histogram(series, bins)
            centers = np.round((edges[:-1] + edges[1:]) / 2, 6)
            series = pd.Series(values, index=centers)
        accum[length].append(series)
    fig, axes = plt.subplots(nrows=len(accum), sharex=True, sharey=True)
    for ax, (length, data) in zip(axes, accum.items(), strict=True):
        df = pd.concat(data, axis=1).fillna(0)
        mean = pd.Series(np.nanmean(df, axis=1), df.index)
        for series in data:
            ax.plot(series.index, series, c="#0072b2", alpha=0.02)
        ax.plot(mean.index, mean, c="#000000", alpha=1.0)
        ax.grid(True, clip_on=False, color='#E0E0E0')
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if xticks is not None:
            ax.set_xticks(xticks)
        if ax is axes[-1]:
            ax.set_xlabel(xlabel)
        library = "Amplicon" if length <= MAX_AMPLICON_LENGTH else "Fragmented"
        ax.set_ylabel(f"Length: {length} nt\n({library})",
                      fontsize=12,
                      labelpad=-360,
                      rotation=0)
        df["Mean"] = mean
        df.to_csv(f"{graph_filename}-{length}.csv")
    fig.set_size_inches(mm_to_in(180), mm_to_in(180))
    fig.text(0.03, 0.50, ylabel, va='center', rotation='vertical', fontsize=12)
    plt.subplots_adjust(left=0.2, right=0.7, top=0.95, bottom=0.1)
    plt.savefig(f"{graph_filename}.{fmt}")
    plt.close()


if __name__ == "__main__":
    fmt = "svg"
    for reads in [200000]:
        sample = f"sample-biased-{reads}"
        graph_avg_of_refs(f"ncov-{reads}",
                          fmt,
                          xlabel="Position",
                          ylabel="Coverage per Position",
                          sample=sample,
                          full_region=True,
                          csv_filename="profile_all_n-count.csv",
                          colname="Informative",
                          pos_base_index=True)
        graph_avg_of_refs(f"mus-hist-{reads}",
                          fmt,
                          xlabel="Mutation Rate",
                          ylabel="Positions per RNA",
                          sample=sample,
                          full_region=False,
                          csv_filename="profile_filtered_m-ratio-q0.csv",
                          colname="Mutated",
                          pos_base_index=True,
                          bins=np.arange(0.0, 0.21, 0.005))
        graph_avg_of_refs(f"nmut-hist-{reads}",
                          fmt,
                          xlabel="Mutations per Read",
                          ylabel="Reads per RNA",
                          sample=sample,
                          full_region=False,
                          csv_filename="histread_filtered_m-count.csv",
                          colname="Mutated",
                          pos_base_index=False,
                          xticks=np.arange(0, 15, 2))
        

