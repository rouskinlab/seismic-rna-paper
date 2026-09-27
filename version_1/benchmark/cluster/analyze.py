import json
import os
import re
from concurrent.futures import ProcessPoolExecutor
from functools import cache
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
            #if ref != "rna-420-1-42":
            #     continue
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


def calc_pearson(x: np.ndarray, y: np.ndarray):
    """ Calculate the Pearson correlation coefficient of two arrays. """
    assert x.shape == y.shape, f"X shape: {x.shape}, Y shape: {y.shape}"
    assert x.size > 0, f"Size: {x.size}"
    # Use pandas instead of np.corrcoef to handle NaN values.
    df = pd.DataFrame({"x": x.flatten(), "y": y.flatten()})
    corr_matrix = df.corr("pearson")
    assert corr_matrix.shape == (2, 2), f"Corrcoef shape: {corr_matrix.shape}"
    return float(corr_matrix.loc["x", "y"])


def assign_clusterings(mus1: np.ndarray, mus2: np.ndarray):
    """ Optimally assign clusters from two groups to each other. """
    n1, k1 = mus1.shape
    n2, k2 = mus2.shape
    assert n1 == n2
    if n1 >= 1 and k1 >= 1 and k2 >= 1:
        # Match the clusters using linear_sum_assignment.
        costs = np.array([[1 - calc_pearson(mus1[:, cluster1],
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
    # Match the clusters from the results and parameters.
    rows, cols = assign_clusterings(mus, real_mus)
    mus = mus[:, rows]
    real_mus = real_mus[:, cols]
    pis = pis[rows]
    real_pis = real_pis[cols]
    assert mus.shape == real_mus.shape
    assert pis.shape == real_pis.shape
    # Calculate the RMSE of the mutation rates.
    mus_corr = calc_pearson(mus, real_mus)
    # Calculate the RMSE of the proportions.
    pis_rmse = calc_rmse(pis, real_pis)
    return mus_corr, pis_rmse


def get_mus_pis_seismic(sample: str, ref: str, reg: str, corrected: bool):
    branch = "" if corrected else "_gap-0"
    # Get the best number of clusters.
    cluster_dir = OUT_SEISMIC_DIR.joinpath(sample, f"cluster{branch}")
    cluster_report_file = cluster_dir.joinpath(ref, reg, "cluster-report.json")
    with open(cluster_report_file) as f:
        cluster_report_data = json.load(f)
    best_k = int(cluster_report_data["Best number of clusters"])
    # Get accuracies of the mutation rates and proportions.
    graph_dir = OUT_SEISMIC_DIR.joinpath(sample, f"graph{branch}")
    mus_file = graph_dir.joinpath(ref, reg, f"profile_clustered-{best_k}-x_m-ratio-q0.csv")
    mus = pd.read_csv(mus_file,
                      index_col=[0, 1],
                      header=[0, 1, 2])
    pis_file = graph_dir.joinpath(ref, reg, "abundance_clustered.csv")
    pis = pd.read_csv(pis_file,
                      index_col=[0],
                      header=[0])
    pis.columns = pis.columns.astype(int)
    pis = pis.T[best_k]
    return best_k, mus, pis


def get_accuracy_seismic(sample: str,
                        ref: str,
                        reg: str,
                        corrected: bool,
                        real_mus: pd.DataFrame,
                        real_pis: pd.Series):
    best_k, mus, pis = get_mus_pis_seismic(sample, ref, reg, corrected)
    mus_corr, pis_rmse = calc_params_accuracy(mus, pis, real_mus, real_pis)
    return best_k, mus_corr, pis_rmse


def get_mus_pis_dreem(sample: str,
                      ref: str,
                      end5: int,
                      end3: int):
    sample_ref = f"{sample}_{ref}"
    em_dir = OUT_DREEM_DIR.joinpath(
        sample_ref, "EM_Clustering", f"{sample_ref}_{end5}_{end3}"
    )
    # Get the best number of clusters.
    best_k = 0
    em_log_file = em_dir.joinpath("log.txt")
    with open(em_log_file) as f:
        for line in f:
            if line.startswith("Predicted number of clusters:"):
                best_k = int(line.split(":")[1].strip())
                break
    if not best_k:
        raise ValueError(f"Best number of clusters not found in {em_log_file}")
    k_dir = em_dir.joinpath(f"K_{best_k}")
    # Get accuracies of the mutation rates.
    em_runs_file = k_dir.joinpath("log_likelihoods.txt")
    em_runs = pd.read_csv(em_runs_file,
                          header=[0],
                          index_col=[0],
                          sep="\t")
    best_runs = [run for run in em_runs.index if run.endswith("-best")]
    if len(best_runs) != 1:
        raise ValueError(f"Expected 1 best run in {em_runs_file}, got {len(best_runs)}")
    best_run = best_runs[0]
    best_run_dir = k_dir.joinpath(f"run_{best_run}")
    mus_file = best_run_dir.joinpath("Clusters_Mu.txt")
    mus = pd.read_csv(mus_file,
                      index_col=[0],
                      header=[0],
                      skiprows=2,
                      sep="\t")
    # Get accuracies of the proportions.
    pis_file = best_run_dir.joinpath("Proportions.txt")
    pis = pd.read_csv(pis_file,
                      index_col=[0],
                      header=[0],
                      sep=",")
    pis = pis[" Real pi "]
    return best_k, mus, pis


def get_accuracy_dreem(sample: str,
                        ref: str,
                        end5: int,
                        end3: int,
                        real_mus: pd.DataFrame,
                        real_pis: pd.Series):
    best_k, mus, pis = get_mus_pis_dreem(sample, ref, end5, end3)
    mus_corr, pis_rmse = calc_params_accuracy(mus, pis, real_mus, real_pis)
    return best_k, mus_corr, pis_rmse


def get_mus_pis_draco(sample: str, ref: str, real_mus: pd.DataFrame):
    # Determine the windows, number of clusters, and proportions of
    # clusters for each window.
    _, true_k, _ = ref_to_attrs(ref)
    stoichs_file = OUT_DRACO_DIR.joinpath(
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
    
    mus = dict()
    pis = dict()
    for window_i, ((end5, end3), window_stoichs) in enumerate(stoichs.items()):
        # Skip windows with a sub-optimal number of clusters.
        if len(window_stoichs) != best_k:
            continue
        # Get the mutation rates for each cluster of this window.
        positions = list(range(end5 + 1, end3 + 2))
        mus[window_i] = dict()
        pis[window_i] = dict()
        for cluster, cluster_stoich in enumerate(window_stoichs, start=1):
            xml_file = OUT_DRACO_DIR.joinpath(
                sample, ref, "rf_norm",
                f"{ref}_{end5}-{end3}_c{cluster - 1}.xml"
            )
            if not xml_file.exists():
                print(f"WARNING: {xml_file} does not exist")
                continue
            tree = ElementTree.parse(xml_file)
            root = tree.getroot()
            tx = root.find('transcript')
            if tx is None:
                raise ValueError("No <transcript> element found in XML.")
            # Parse the transcript ID to get start/end (zero‑indexed, half-open)
            # e.g. id="rna-140-2-1_21-118_c0" → parts[1] == "21-118"
            tx_id = tx.get('id', '')
            parts = tx_id.split('_')
            if len(parts) < 2 or '-' not in parts[1]:
                raise ValueError(f"Unexpected transcript ID format: {tx_id}")
            if parts[0] != ref:
                raise ValueError(f"Reference mismatch: {tx_id} != {ref}")
            # Extract and clean the sequence
            seq_elem = tx.find('sequence')
            if seq_elem is None or seq_elem.text is None:
                raise ValueError("No <sequence> text found.")
            seq = ''.join(seq_elem.text.split())  # remove whitespace/newlines
            if len(seq) != len(positions):
                raise ValueError(f"Sequence length {len(seq)} ≠ number of positions {len(positions)}")
            # Extract and parse the reactivity values
            react_elem = tx.find('reactivity')
            if react_elem is None or react_elem.text is None:
                raise ValueError("No <reactivity> text found.")
            # Split on commas, strip whitespace, drop any empty trailing entries
            raw_vals = [v.strip() for v in react_elem.text.replace('\n', '').split(',') if v.strip()]
            if len(raw_vals) != len(positions):
                raise ValueError(f"Found {len(raw_vals)} reactivity entries but expected {len(positions)}")
            # Convert to a pandas Series.
            mus[window_i][cluster] = pd.Series(
                list(map(float, raw_vals)),
                index=pd.MultiIndex.from_arrays([positions, list(seq)],
                                                names=['Position', 'Base'])
            )
            pis[window_i][cluster] = cluster_stoich
        mus[window_i] = pd.DataFrame.from_dict(mus[window_i], orient="columns")
        pis[window_i] = pd.Series(pis[window_i])
        mus[window_i].columns = pd.Index(range(1, len(mus[window_i].columns) + 1))
        pis[window_i].index = mus[window_i].columns
        # Write the mutation rates and proportions to CSV files.
        mus_csv_file = OUT_DRACO_DIR.joinpath(
            sample, ref, "rf_norm",
            f"{ref}_{end5}-{end3}_mus.csv"
        )
        mus[window_i].to_csv(mus_csv_file)
        pis_csv_file = OUT_DRACO_DIR.joinpath(
            sample, ref, "rf_norm",
            f"{ref}_{end5}-{end3}_pis.csv"
        )
        pis[window_i].to_csv(pis_csv_file)
    # Average the mutation rates and proportions across windows.
    mus_avg = pd.DataFrame(np.nan,
                           index=real_mus.index,
                           columns=pd.MultiIndex.from_product([range(1, true_k + 1),
                                                               range(len(mus))],
                                                               names=["Cluster", "Window"]))
    pis_avg = pd.Series(np.nan, index=mus_avg.columns)
    for window_i, window_mus in mus.items():
        if window_mus.size == 0:
            continue
        window_mus.index = window_mus.index.get_level_values("Position")
        window_mus = window_mus.reindex(real_mus.index, fill_value=np.nan)
        use_mus = ~np.any(np.isnan(window_mus), axis=1)
        cxs, cys = assign_clusterings(real_mus.values[use_mus], window_mus.values[use_mus])
        for cx, cy in zip(cxs, cys, strict=True):
            mus_avg.loc[window_mus.index.get_level_values("Position"),
                        (cx + 1, window_i)] = window_mus.iloc[:, cy]
            pis_avg.loc[cx + 1] = pis[window_i].iloc[cy]
    mus_avg = mus_avg.T.groupby(level="Cluster").mean().T
    pis_avg = pis_avg.groupby(level="Cluster").mean()
    mus_avg = mus_avg.dropna(axis=1, how="all")
    mus_avg.columns = list(range(1, mus_avg.columns.size + 1))
    pis_avg = pis_avg.dropna()
    pis_avg.index = list(range(1, pis_avg.index.size + 1))
    assert pis_avg.index.equals(mus_avg.columns)
    return best_k, mus_avg, pis_avg


def get_accuracy_draco(sample: str,
                    ref: str,
                    real_mus: pd.DataFrame,
                    real_pis: pd.Series):
    best_k, mus, pis = get_mus_pis_draco(sample, ref, real_mus)
    if mus.size == 0 or pis.size == 0:
        return best_k, np.nan, np.nan
    mus_corr, pis_rmse = calc_params_accuracy(mus, pis, real_mus, real_pis)
    return best_k, mus_corr, pis_rmse


def read_mus_pis_dance(sample: str, ref: str):
    reactivities_file = OUT_DANCE_DIR.joinpath(
            "dancemapper", f"{sample}_{ref}-reactivities.txt"
    )
    with open(reactivities_file) as f:
        match = re.match("^([0-9]+) components; BIC=[0-9.-]+$",
                            f.readline().rstrip())
        if not match:
            raise ValueError(reactivities_file)
        best_k = int(match.groups()[0])
        pis_line = f.readline().rstrip()
        if not pis_line.startswith("p "):
            raise ValueError(reactivities_file)
        pis = pd.Series(list(map(float, pis_line.split()[1:])),
                        index=range(1, best_k + 1))
        f.readline()  # Skip the mutation rates header line.
        mus_data = dict()
        for line in f.readlines():
            line_data = line.rstrip().rstrip(" i").split()
            if len(line_data) != 3 + 2 * best_k:
                raise ValueError(reactivities_file)
            pos = int(line_data[0])
            base = line_data[1].replace("U", "T")
            mus_data[pos, base] = dict()
            for k in range(1, best_k + 1):
                mu = float(line_data[1 + 2 * k])
                mus_data[pos, base][k] = mu
    mus = pd.DataFrame.from_dict(mus_data, orient="index")
    mus.index.names = ["Position", "Base"]
    return best_k, mus, pis


def get_accuracy_dance(sample: str,
                    ref: str,
                    real_mus: pd.DataFrame,
                    real_pis: pd.Series):
    best_k, mus, pis = read_mus_pis_dance(sample, ref)
    mus_corr, pis_rmse = calc_params_accuracy(mus, pis, real_mus, real_pis)
    return best_k, mus_corr, pis_rmse


def get_sample_results(ref: str,
                       length: int,
                       true_k: int,
                       rep: int,
                       sample: str,
                       is_biased: bool,
                       n_reads_total: int,
                       real_mus: pd.DataFrame,
                       real_pis: pd.Series):
    print(f"Processing {sample} {ref} {length} {true_k} {rep}")
    results = list()

    # SEISMIC
    print("  SEISMIC")
    best_k, mus_corr, pis_rmse = get_accuracy_seismic(
        sample, ref, get_region_name(length), True, real_mus, real_pis
    )
    results.append({"Length": length,
                    "True K": true_k,
                    "Biased": is_biased,
                    "Total Reads": n_reads_total,
                    "Replicate": rep,
                    "Method": "SEISMIC",
                    "Result K": best_k,
                    "Correct K": best_k == true_k,
                    "Mus PCC": mus_corr,
                    "Pis RMSE": pis_rmse})

    # SEISMIC (no bias correction)
    print("  SEISMIC (no bias correction)")
    best_k, mus_corr, pis_rmse = get_accuracy_seismic(
        sample, ref, get_region_name(length), False, real_mus, real_pis
    )
    results.append({"Length": length,
                    "True K": true_k,
                    "Biased": is_biased,
                    "Total Reads": n_reads_total,
                    "Replicate": rep,
                    "Method": "SEISMIC (no bias correction)",
                    "Result K": best_k,
                    "Correct K": best_k == true_k,
                    "Mus PCC": mus_corr,
                    "Pis RMSE": pis_rmse})

    # DANCE
    print("  DANCE")
    best_k, mus_corr, pis_rmse = get_accuracy_dance(
        sample, ref, real_mus, real_pis
    )
    results.append({"Length": length,
                    "True K": true_k,
                    "Biased": is_biased,
                    "Total Reads": n_reads_total,
                    "Replicate": rep,
                    "Method": "DANCE",
                    "Result K": best_k,
                    "Correct K": best_k == true_k,
                    "Mus PCC": mus_corr,
                    "Pis RMSE": pis_rmse})

    # DRACO
    print("  DRACO")
    best_k, mus_corr, pis_rmse = get_accuracy_draco(
        sample, ref, real_mus, real_pis
    )
    results.append({"Length": length,
                    "True K": true_k,
                    "Biased": is_biased,
                    "Total Reads": n_reads_total,
                    "Replicate": rep,
                    "Method": "DRACO",
                    "Result K": best_k,
                    "Correct K": best_k == true_k,
                    "Mus PCC": mus_corr,
                    "Pis RMSE": pis_rmse})

    # DREEM
    if length <= MAX_AMPLICON_LENGTH:
        print("  DREEM")
        end5 = PRIMER_LENGTH + 1
        end3 = length - PRIMER_LENGTH
        best_k, mus_corr, pis_rmse = get_accuracy_dreem(
            sample, ref, end5, end3, real_mus, real_pis
        )
        results.append({"Length": length,
                        "True K": true_k,
                        "Biased": is_biased,
                        "Total Reads": n_reads_total,
                        "Replicate": rep,
                        "Method": "DREEM",
                        "Result K": best_k,
                        "Correct K": best_k == true_k,
                        "Mus PCC": mus_corr,
                        "Pis RMSE": pis_rmse})
    return results


def get_ref_results(ref: str,
                       length: int,
                       true_k: int,
                       rep: int):
    real_mus, real_pis = get_real_params(ref)
    real_pis = real_pis.loc[true_k]
    results = list()
    for sample, (is_biased, n_reads_total) in get_samples().items():
        results.extend(get_sample_results(ref, length, true_k, rep,
                                          sample, is_biased, n_reads_total,
                                          real_mus, real_pis))
    return results


def collect_results():
    results = list()
    args = [(ref, length, true_k, rep)
            for ref, (length, true_k, rep) in get_refs().items()]
    # for arg in args:
    #     results.extend(get_ref_results(*arg))
    max_workers = 1
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        print(f"Opened {executor}")
        for result in executor.map(get_ref_results, *zip(*args)):
            results.extend(result)
        print(f"Closed {executor}")
    return pd.DataFrame.from_records(results)


METHOD_STYLES = {"SEISMIC": {"color": "#0072b2",
                             "marker": "o",
                             "markersize": 3.0 * 1.125**0,
                             "linewidth": 1.0 * 1.125**0,
                             "zorder": 14},
                 "SEISMIC (no bias correction)": {"color": "#56b4e9",
                                                  "marker": "o",
                                                  "markersize": 3.0 * 1.125**0,
                                                  "linewidth": 1.0 * 1.125**0,
                                                  "zorder": 13},
                 "DANCE": {#"color": "#009e73",
                           "color": "#7fceb8",
                           "marker": "o",
                           "markersize": 3.0 * 1.125**0,
                           "linewidth": 1.0 * 1.125**0,
                           "zorder": 12},
                 "DRACO": {#"color": "#e69f00",
                           "color": "#f2ce7f",
                           "marker": "o",
                           "markersize": 3.0 * 1.125**0,
                           "linewidth": 1.0 * 1.125**0,
                           "zorder": 11},
                 "DREEM": {#"color": "#cc79a7",
                           "color": "#e5bbd2",
                           "marker": "o",
                           "markersize": 3.0 * 1.125**0,
                           "linewidth": 1.0 * 1.125**0,
                           "zorder": 10}}


def calc_ci(values: pd.Series, alpha: float):
    """ Calculate the confidence interval for a given series of values. """
    n = values.size
    if n == 0:
        return np.nan, 0.0, 0.0
    mean = values.mean()
    # Calculate the error depending on the type of the column.
    if pd.api.types.is_bool_dtype(values):
        # Calculate confidence interval using binomial (Clopper-Pearson) interval
        k = values.sum()
        # Clopper-Pearson interval (using beta distribution)
        error_lo = mean - stats.beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
        error_up = stats.beta.ppf(1 - alpha / 2, k + 1, n - k) - mean if k < n else 0.0
    elif pd.api.types.is_numeric_dtype(values):
        # Calculate confidence interval using t-distribution interval
        stderr = values.std(ddof=1) / np.sqrt(n)
        critical_value = stats.t.ppf(1 - alpha / 2, df=n-1)
        error_lo = critical_value * stderr
        error_up = error_lo
    else:
        raise TypeError(f"Unknown column type: {values.dtype}")
    return mean, error_lo, error_up


def graph_vs_length(results: pd.DataFrame, 
                    col_name: str,
                    title: str,
                    fname: str,
                    alpha: float = 0.05):
    """ Graph the proportion of correct K for each method. """
    lengths = np.sort(results["Length"].unique())
    true_ks = np.sort(results["True K"].unique())
    total_reads = np.sort(results["Total Reads"].unique())
    methods = results["Method"].unique()
    assert np.array_equal(true_ks, np.arange(1, true_ks.size + 1))

    fig, axs = plt.subplots(nrows=true_ks.size,
                            ncols=lengths.size,
                            squeeze=False,                            
                            sharex=True,
                            sharey=True,
                            gridspec_kw={'hspace': 0.2,
                                         'wspace': 0.2},
                            figsize=(mm_to_in(210),
                                     mm_to_in(195)))
    
    for row, true_k in enumerate(true_ks):
        for col, length in enumerate(lengths):
            ax = axs[row, col]

            ax_data = results[(results["Length"] == length) & (results["True K"] == true_k)]
            for method in methods:
                method_data = ax_data[ax_data["Method"] == method]
                x = list()
                y = list()
                l = list()
                u = list()
                for nreads in sorted(method_data["Total Reads"].unique()):
                    nreads_data = method_data[method_data["Total Reads"] == nreads]
                    value, error_lo, error_up = calc_ci(nreads_data[col_name], alpha)
                    x.append(nreads)
                    y.append(value)
                    l.append(error_lo)
                    u.append(error_up)
                ax.errorbar(x, y,
                            yerr=np.vstack([l, u]),
                            capsize=3,
                            clip_on=False,
                            label=method,
                            **METHOD_STYLES[method])
                
            ax.set_xscale("log")
            ax.set_xlim(np.min(total_reads), np.max(total_reads))
            if col_name == "Mus PCC":
                ax.set_ylim(0.5, 1.0)
                ax.set_yticks(np.linspace(0.5, 1.0, 6))
            else:
                ax.set_ylim(0.0, 1.0)
                ax.set_yticks(np.linspace(0, 1, 6))
            # Format x-axis with decimal numbers instead of powers
            ax.xaxis.set_major_formatter(plt.FuncFormatter(
                lambda x, _: f"{x / 1000:,.0f}"
            ))
            log_x_range = np.log10(ax.get_xlim()[1]) - np.log10(ax.get_xlim()[0])
            y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
            ax.set_aspect(0.8 * log_x_range / y_range)

            ax.grid(True, clip_on=False, color='#E0E0E0')
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if ax.get_legend():
                ax.get_legend().remove()
            
            if row == 0:
                library = "Amplicon" if length <= MAX_AMPLICON_LENGTH else "Fragmented"
                ax.set_title(f"Length: {length} nt\n({library})",
                             fontsize=12,
                             pad=10)
            if col == lengths.size - 1:                
                ax.set_ylabel(f"True Clusters: {true_k}",
                              fontsize=12,
                              labelpad=-180,
                              rotation=0)
    
    fig.text(0.45, 0.03, "Number of Reads (1,000s)", ha='center', fontsize=12)
    fig.text(0.02, 0.50, title, va='center', rotation='vertical', fontsize=12)

    plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)

    plt.savefig(fname)
    plt.close()


def graph_proportion_each_k(results: pd.DataFrame, fname: str):
    """ Graph the proportion of each K for each method. """
    lengths = np.sort(results["Length"].unique())
    true_ks = np.sort(results["True K"].unique())
    methods = results["Method"].unique()
    assert np.array_equal(true_ks, np.arange(1, true_ks.size + 1))

    fig, axs = plt.subplots(nrows=true_ks.size,
                            ncols=lengths.size,
                            squeeze=False,                            
                            sharex=True,
                            sharey=True,
                            gridspec_kw={'hspace': 0.2,
                                         'wspace': 0.2},
                            figsize=(mm_to_in(210),
                                     mm_to_in(195)))
    
    for row, true_k in enumerate(true_ks):
        for col, length in enumerate(lengths):
            ax = axs[row, col]

            ax_data = results[(results["Length"] == length) & (results["True K"] == true_k)]
            for k in range(1, true_ks.max() + 2):
                if k != true_k:
                    # Add a rectangle highlighting the true k value
                    rect = plt.Rectangle((k - 0.5, 0), 1, 1,
                                        facecolor='#F0F0F0',
                                        edgecolor='none', 
                                        zorder=0,
                                        clip_on=False)
                    ax.add_patch(rect)

            for method in methods:
                fraction_each_k = dict()
                method_data = ax_data[ax_data["Method"] == method]
                for k in range(1, true_ks.max() + 1):
                    fraction_each_k[k] = np.mean(method_data["Result K"] == k)
                fraction_each_k[true_ks.max() + 1] = np.mean(method_data["Result K"] > true_ks.max())
                fraction_each_k = pd.Series(fraction_each_k)
                ax.plot(fraction_each_k.index,
                        fraction_each_k.values,
                        clip_on=False,
                        label=method,
                        **METHOD_STYLES[method])
            
            ax.set_xlim(true_ks.min() - 0.5, true_ks.max() + 1.5)
            ax.set_ylim(0, 1)
            ax.set_xticks(np.arange(true_ks.min(), true_ks.max() + 2))
            ax.set_yticks(np.linspace(0, 1, 6))
            ax.xaxis.set_major_formatter(plt.FuncFormatter(
                lambda x, _: f"{x}+" if x > true_ks.max() else str(x)
            ))
            
            x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
            y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
            ax.set_aspect(0.8 * x_range / y_range)

            ax.grid(True, clip_on=False, color='#E0E0E0')
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if ax.get_legend():
                ax.get_legend().remove()
            
            if row == 0:
                library = "Amplicon" if length <= MAX_AMPLICON_LENGTH else "Fragmented"
                ax.set_title(f"Length: {length} nt\n({library})",
                             fontsize=12,
                             pad=10)
            if col == lengths.size - 1:
                ax.set_ylabel(f"True Clusters: {true_k}",
                              fontsize=12,
                              labelpad=-180,
                              rotation=0)
    
    fig.text(0.45, 0.03, "Number of Clusters Detected", ha='center', fontsize=12)
    fig.text(0.02, 0.50, "Proportion of Simulations", 
             va='center', rotation='vertical', fontsize=12)

    plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)

    plt.savefig(fname)
    plt.close()


def graph_accuracy(results: pd.DataFrame, col_name: str, fname: str):
    """ Graph the accuracy for each method. """
    lengths = np.sort(results["Length"].unique())
    true_ks = np.sort(results["True K"].unique())

    fig, axs = plt.subplots(nrows=lengths.size,
                            ncols=true_ks.size,
                            squeeze=False,                            
                            sharex=True,
                            sharey=True,
                            gridspec_kw={'hspace': 0.2,
                                         'wspace': 0.2},
                            figsize=(8, 4.5))
    
    for row, length in enumerate(lengths):
        for col, true_k in enumerate(true_ks):
            print(f"Graphing {col_name} for length={length} and K={true_k}")
            ax = axs[row, col]
            ax_data = results[(results["Length"] == length) & (results["True K"] == true_k)]

            # Plot each method separately to avoid the deprecation warning
            for i, (method, style) in enumerate(METHOD_STYLES.items()):
                method_data = ax_data[ax_data["Method"] == method][col_name]
                if method_data.size == 0:
                    continue
                color = style['color']
                # Position each stripplot at the x-coordinate corresponding to the method index
                sns.stripplot(x=np.full(method_data.size, i), y=method_data, ax=ax, color=color, 
                              zorder=5, size=2.0, alpha=0.5)
                # Add a horizontal line at the mean
                mean = method_data.mean()
                ax.hlines(mean, i-0.2, i+0.2, color="black", linewidth=1)
            
            # Set x-ticks and labels manually
            ax.set_xticks(range(len(METHOD_STYLES)))
            ax.set_xticklabels(METHOD_STYLES.keys())
            plt.setp(ax.get_xticklabels(), rotation=60, ha='right')
            ax.set_yticks(np.linspace(0., 1., 5))
            
            #ax.set_yscale("log")
            if col_name == "Mus PCC":
                ax.set_ylim(0.5, 1.05)
            elif col_name == "Pis RMSE":
                ax.set_ylim(0, 1.05)
            else:
                raise ValueError(f"Unknown column name: {col_name}")
            ax.set_ylabel("")
            ax.set_xlabel("")

            ax.grid(True, clip_on=False, color='#E0E0E0', zorder=0)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if ax.get_legend():
                ax.get_legend().remove()
            
            if row == 0:
                ax.set_title(f"True Clusters: {true_k}",
                             fontsize=12,
                             pad=10)
            if col == true_ks.size - 1:
                library = "Amplicon" if length <= MAX_AMPLICON_LENGTH else "Fragmented"
                ax.set_ylabel(f"Length: {length} nt\n({library})",
                              fontsize=12,
                              labelpad=-144,
                              rotation=0)
    
    fig.text(0.02, 0.50, col_name,
             va='center', rotation='vertical', fontsize=12)

    plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.15)

    plt.savefig(fname)
    plt.close()


def graph_mus_corr(results: pd.DataFrame, fname: str):
    graph_accuracy(results, "Mus PCC", fname)


def graph_pis_rmse(results: pd.DataFrame, fname: str):
    graph_accuracy(results, "Pis RMSE", fname)


def _line_graph_fraction(results: pd.DataFrame, 
                         fname: str, 
                         y_label: str, 
                         column: str,
                         alpha: float = 0.05):
    """ Graph the fraction of correct K predictions for each method. """

    lengths = np.sort(results["Length"].unique())
    true_ks = np.sort(results["True K"].unique())

    fig, axs = plt.subplots(ncols=lengths.size, sharex=True, sharey=True, figsize=(9, 4))

    for ax, length in zip(axs, lengths, strict=True):
        library = "Amplicon" if length <= MAX_AMPLICON_LENGTH else "Fragmented"
        ax.set_title(f"{length} nt\n({library})", fontsize=12)

        # Calculate fraction of correct K for each method and true K
        correct_fractions = []
        for method in METHOD_STYLES:
            for true_k in true_ks:
                subset = results.loc[(results["Method"] == method)
                                     & (results["True K"] == true_k)
                                     & (results["Length"] == length)]
                if len(subset) > 0:
                    value, error_lo, error_up = calc_ci(subset[column], alpha)
                    correct_fractions.append({
                        "Method": method,
                        "True K": true_k,
                        "Value": value,
                        "ErrorLo": error_lo,
                        "ErrorUp": error_up,
                    })
        # Convert to DataFrame for plotting
        plot_data = pd.DataFrame(correct_fractions)
        # Create line plot
        for method in reversed(METHOD_STYLES):
            method_data = plot_data[plot_data["Method"] == method]
            ax.errorbar(method_data["True K"],
                        method_data["Value"],
                        yerr=np.vstack([method_data["ErrorLo"],
                                        method_data["ErrorUp"]]),
                        label=method,
                        **METHOD_STYLES[method],
                        clip_on=False,
                        capsize=3)
        ax.set_xticks(true_ks)
        ax.set_xticklabels(true_ks)
        ax.set_xlim(true_ks.min() - 0.5, true_ks.max() + 0.5)
        if column == "Correct K":
            ax.set_ylim(0, 1)
            ax.set_yticks(np.linspace(0, 1, 6))
        elif column == "Mus PCC":
            ax.set_ylim(0.8, 1)
            ax.set_yticks(np.linspace(0.8, 1, 3))
        elif column == "Pis RMSE":
            ax.set_ylim(0, 0.5)
            ax.set_yticks(np.linspace(0, 0.5, 6))
        else:
            raise ValueError(f"Unknown column name: {column}")
        ax.set_xlabel("True Number of Clusters", fontsize=12)
        ax.set_ylabel(y_label, fontsize=12)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, clip_on=False, color='#E0E0E0')
        x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        ax.set_aspect(x_range / y_range)
    
    plt.tight_layout()
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.close()


def line_graph_fraction_correct_k(results: pd.DataFrame, fname: str):
    _line_graph_fraction(results,
                         fname,
                         "Fraction of Simulations with Correct Clusters",
                         "Correct K")


def line_graph_average_pcc(results: pd.DataFrame, fname: str):
    _line_graph_fraction(results,
                         fname,
                         "Average Correlation of Mutation Rates",
                         "Mus PCC")


def line_graph_average_pis_rmse(results: pd.DataFrame, fname: str):
    _line_graph_fraction(results,
                         fname,
                         "Average RMSE of Proportions",
                         "Pis RMSE")


if __name__ == "__main__":

    format = "svg"
    if not RESULTS_FILE.exists():
        results = collect_results()
        results.to_csv(RESULTS_FILE, index=False)
    else:
        results = pd.read_csv(RESULTS_FILE)

    is_biased = results["Biased"]
    only_main_methods = results["Method"].isin(["SEISMIC", "DANCE", "DRACO", "DREEM"])
    seismic_modes = results["Method"].isin(["SEISMIC", "SEISMIC (no bias correction)"])
    for selector, methods in [("main", only_main_methods), ("seismic", seismic_modes)]:
        is_biased_methods = is_biased & methods
        graph_vs_length(results.loc[is_biased_methods],
                        "Correct K",
                        "Proportion of Simulations with Correct Clusters",
                        f"vs-length-correct-k-{selector}.{format}")
        graph_vs_length(results.loc[is_biased_methods],
                        "Mus PCC",
                        "Average Pearson Correlation of Mutation Rates",
                        f"vs-length-mus-pcc-{selector}.{format}")
        graph_vs_length(results.loc[is_biased_methods],
                        "Pis RMSE",
                        "Average RMSE of Cluster Proportions",
                        f"vs-length-pis-rmse-{selector}.{format}")
        for reads in [5000, 200000]:
            label = f"{selector}-{reads}"
            rows = is_biased_methods & (results["Total Reads"] == reads)
            line_graph_fraction_correct_k(results.loc[rows],
                                          f"fraction-correct-k_{label}.{format}")
            line_graph_average_pcc(results.loc[rows],
                                   f"average-pcc_{label}.{format}")
            line_graph_average_pis_rmse(results.loc[rows],
                                        f"average-pis-rmse_{label}.{format}")
            #graph_mus_corr(results.loc[rows], f"mus-corr-{label}.{format}")
            #graph_pis_rmse(results.loc[rows], f"pis-rmse-{label}.{format}")
            graph_proportion_each_k(results.loc[rows], f"proportion-each-k-{label}.{format}")

