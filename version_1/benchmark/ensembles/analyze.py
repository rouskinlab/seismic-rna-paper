from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pandas as pd
from matplotlib import font_manager, pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from seismicrna.core.mu import compare_windows
from seismicrna.core.rna import from_ct, find_enclosing_pairs

# Configure matplotlib to write text as text in SVG files, not as paths
plt.rcParams['svg.fonttype'] = 'none'

# Add Helvetica Neue and Helvetica Neue Light to matplotlib's font list
font_dirs = []  # Add custom font directories if needed
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

# Set font family to Helvetica Neue for all text elements
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial', 'sans-serif']
plt.rcParams['font.size'] = 12


SIM_DIR = Path("sim")
PARAMS_DIR = SIM_DIR.joinpath("params")


def get_ref_name(refnum: int) -> str:
    return f"long-rna-{refnum}"


def get_ref_full_dir(refnum: int) -> Path:
    return PARAMS_DIR.joinpath(get_ref_name(refnum), "full")


def get_ct_path(refnum: int) -> Path:
    return get_ref_full_dir(refnum).joinpath("simulated.ct")


def get_structures_per_position_path(refnum: int) -> Path:
    return get_ref_full_dir(refnum).joinpath("structures_per_position.csv")


def get_structures_per_region_path(refnum: int) -> Path:
    return get_ref_full_dir(refnum).joinpath("structures_per_region.csv")


def count_structures_per_position(refnum: int):
    ct_path = get_ct_path(refnum)
    structures_per_position_path = get_structures_per_position_path(refnum)
    if structures_per_position_path.is_file():
        return pd.read_csv(structures_per_position_path, index_col=[0, 1])
    structures = list(from_ct(ct_path))
    if not structures:
        return pd.Series([0])
    index = structures[0].table.index
    enclosing_pairs_sets = defaultdict(set)
    for i, structure in enumerate(structures, start=1):
        print(f"Processing structure {i} of {len(structures)}")
        assert structure.table.index.equals(index)
        enclosing_pairs = find_enclosing_pairs(structure.table)
        for row in index:
            pair = (int(enclosing_pairs.loc[row, "5' End"]),
                    int(enclosing_pairs.loc[row, "3' End"]))
            enclosing_pairs_sets[row].add(pair)
    structures_per_position = pd.Series({row: len(pairs) for row, pairs in enclosing_pairs_sets.items()})
    structures_per_position.index.names = index.names
    structures_per_position.to_csv(structures_per_position_path)
    return structures_per_position


def count_structures_per_region(refnum: int):
    structures_per_region_path = get_structures_per_region_path(refnum)
    if structures_per_region_path.is_file():
        return pd.read_csv(structures_per_region_path, index_col=[0, 1])
    ct_files = list(PARAMS_DIR.glob(f"{get_ref_name(refnum)}-region-*/full/simulated.ct"))
    ct_files = sorted(ct_files, key=lambda ct_file: int(ct_file.parent.parent.name.split("-")[-1]))
    structures_per_region = list()
    offset = 0
    for ct_file in ct_files:
        structures = list(from_ct(ct_file))
        if len(structures) == 0:
            raise ValueError(f"No structures found in {ct_file}")
        index = pd.MultiIndex.from_arrays([
            structures[0].table.index.get_level_values("Position") + offset,
            structures[0].table.index.get_level_values("Base")
        ])
        offset += index.size
        structures_per_region.append(pd.Series(len(structures), index))
    structures_per_region = pd.concat(structures_per_region, axis=0)
    structures_per_region.to_csv(structures_per_region_path)
    return structures_per_region


def plot_structures_per_position(structures_per_position: pd.Series):
    plt.plot(structures_per_position.index.get_level_values("Position"),
             structures_per_position)
    plt.show()


def get_real_mus_per_position(refnum: int,
                              positions: pd.Index | None = None,
                              tol: float = 1e-6):
    mus_file = get_ref_full_dir(refnum).joinpath("mus.csv")
    mus = pd.read_csv(mus_file, index_col=[0, 1], header=[0, 1])
    if positions is not None:
        mus = mus.loc[positions]
    # Merge columns with very similar values (max abs diff <= tol)
    similar = pd.DataFrame(np.less_equal(np.max(np.abs(np.subtract(np.expand_dims(mus.to_numpy(), axis=1),
                                                               np.expand_dims(mus.to_numpy(), axis=2))),
                                                axis=0),
                                         tol),
                           index=mus.columns,
                           columns=mus.columns)
    mus_nonredundant = dict()
    for new, values in mus.items():
        for exist in mus_nonredundant.keys():
            if similar.loc[new, exist]:
                break
        else:
            mus_nonredundant[new] = values
    mus_nonredundant = pd.DataFrame.from_dict(mus_nonredundant, orient="columns").copy(deep=True)
    mus_nonredundant.index = mus_nonredundant.index.remove_unused_levels()
    mus_nonredundant.columns.names = mus.columns.names
    return mus_nonredundant


def get_calculated_mus_seismic(reads: int, refnum: int):
    ref_dir = Path("out-seismic", f"sample-biased-{reads}", "graph", get_ref_name(refnum))
    profiles = list(ref_dir.glob("*/profile_clustered-*-x_m-ratio-q0.csv"))
    all_mus = dict()
    for profile in profiles:
        region = profile.parent.name
        if region in all_mus:
            raise ValueError(f"Region {region} already in all_mus")
        mus = pd.read_csv(profile, index_col=[0, 1], header=[0, 1, 2])["Mutated"]
        all_mus[region] = mus
    return all_mus


def get_calculated_mus_draco(reads: int, refnum: int, winlen: int):
    ref_name = get_ref_name(refnum)
    rf_norm_dir = Path("out-draco", f"sample-biased-{reads}", ref_name, f"rf_norm-{winlen}")
    xml_files = list(rf_norm_dir.glob(f"{ref_name}_*-*_c*.xml"))
    all_mus = dict()
    for xml_file in xml_files:
        profile_name = xml_file.stem
        refstr, regstr, clstr = profile_name.split("_")
        assert refstr == ref_name
        start, end = map(int, regstr.split("-"))
        positions = np.arange(start + 1, end + 2)
        region = f"{start}-{end}"
        if region not in all_mus:
            all_mus[region] = dict()
        cluster = int(clstr[1:])
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
        if parts[0] != ref_name:
            raise ValueError(f"Reference mismatch: {tx_id} != {ref_name}")
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
        all_mus[region][cluster] = pd.Series(
            list(map(float, raw_vals)),
            index=pd.MultiIndex.from_arrays([positions, list(seq)],
                                            names=['Position', 'Base'])
        )
    return {region: pd.DataFrame.from_dict(clusters, orient="columns")
            for region, clusters in all_mus.items()}


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
    # Match the clusters using linear_sum_assignment.
    costs = np.array([[1 - calc_pearson(mus1[:, cluster1],
                                        mus2[:, cluster2])
                       for cluster2 in range(k2)]
                      for cluster1 in range(k1)]).reshape((k1, k2))
    if costs.size > 0 and np.all(~np.isnan(costs)):
        from scipy.optimize import linear_sum_assignment
        rows, cols = linear_sum_assignment(costs)
    else:
        # If n1 == 0, then the costs matrix will contain NaN, which will
        # cause linear_sum_assignment to raise an error.
        rows = np.arange(k1)
        cols = np.arange(k1)
    return rows, cols


def calc_params_accuracy(mus: pd.DataFrame,
                         real_mus: pd.DataFrame,
                         method: str = "pcc",
                         size: int = 45):
    """ Calculate accuracies of the mutation rates and proportions. """
    # Make sure the indexes are aligned.
    assert mus.index.equals(real_mus.index)
    # Match the clusters from the results and parameters.
    rows, cols = assign_clusterings(mus.values, real_mus.values)
    n = rows.size
    assert 1 <= n == cols.size
    # Calculate the rolling correlation of the mutation rates.
    for i in range(n):
        index1 = mus.iloc[:, rows[i]].index
        index2 = real_mus.iloc[:, cols[i]].index
        assert index1.equals(index2)
    mus_corr = pd.DataFrame.from_dict(
        {i: compare_windows(mus.iloc[:, rows[i]],
                            real_mus.iloc[:, cols[i]],
                            method=method,
                            size=size)
         for i in range(n)},
        orient="columns"
    )
    return mus_corr.mean(axis=1)


def get_accuracy_seismic(reads: int, refnum: int):
    all_mus = get_calculated_mus_seismic(reads, refnum)
    all_mus_corr = list()
    all_ks = list()
    for region, mus in all_mus.items():
        print("REGION", region)
        real_mus = get_real_mus_per_position(refnum, mus.index)
        mus_corr = calc_params_accuracy(mus, real_mus)
        mus_corr.index = mus_corr.index.get_level_values("Position")
        all_mus_corr.append(mus_corr)
        all_ks.append(pd.Series(mus.shape[1], index=mus_corr.index))
    all_mus_corr = pd.concat(all_mus_corr, axis=0)
    all_ks = pd.concat(all_ks, axis=0)
    return all_mus_corr, all_ks


def get_accuracies_seismic(reads: int, refnums: list[int]):
    all_mus_corr = dict()
    all_ks = dict()
    for refnum in refnums:
        print("REFERENCE", refnum)
        all_mus_corr[refnum], all_ks[refnum] = get_accuracy_seismic(reads, refnum)
    return (pd.DataFrame.from_dict(all_mus_corr, orient="columns"), 
            pd.DataFrame.from_dict(all_ks, orient="columns"))


def get_accuracy_draco(reads: int, refnum: int, winlen: int):
    all_mus = get_calculated_mus_draco(reads, refnum, winlen)
    all_mus_corr = defaultdict(dict)
    all_ks = defaultdict(dict)
    for region, region_mus in all_mus.items():
        print("REGION", region)
        end5, end3 = map(int, region.split("-"))
        center = (end5 + end3) // 2
        real_mus = get_real_mus_per_position(refnum, region_mus.index)
        mus_corr = calc_params_accuracy(region_mus, real_mus)
        k = region_mus.shape[1]
        for (pos, base), corr in mus_corr.items():
            distance = abs(pos - center)
            all_mus_corr[pos][distance] = corr
            all_ks[pos][distance] = k
    # For each position, pick the region with the closest center.
    for pos in list(all_mus_corr.keys()):
        min_distance = min(all_mus_corr[pos].keys())
        all_mus_corr[pos] = all_mus_corr[pos][min_distance]
        all_ks[pos] = all_ks[pos][min_distance]
    all_mus_corr = pd.Series(all_mus_corr)
    all_mus_corr.index.name = "Position"
    all_ks = pd.Series(all_ks)
    all_ks.index.name = "Position"
    return all_mus_corr, all_ks


def get_accuracies_draco(reads: int, refnums: list[int], winlen: int):
    all_mus_corr = dict()
    all_ks = dict()
    for refnum in refnums:
        print("REFERENCE", refnum)
        all_mus_corr[refnum], all_ks[refnum] = get_accuracy_draco(reads, refnum, winlen)
    return (pd.DataFrame.from_dict(all_mus_corr, orient="columns"),
            pd.DataFrame.from_dict(all_ks, orient="columns"))


def plot_r(df: pd.DataFrame,
                 color: str,
                 filename: str):
    fig, ax = plt.subplots()
    fig.set_size_inches(3.6, 3.0)
    alpha = 0.1
    for refnum in df.columns:
        col = df.loc[:, refnum]
        # Negative correlations (below the x-axis) are not shown.
        col = col.mask(col < 0, np.nan)
        ax.plot(col.index,
                col,
                c=color,
                alpha=alpha,
                linewidth=1,
                clip_on=False)
    ax.plot(df.index,
            np.nanmean(df, axis=1),
            c="#000000",
            alpha=1,
            linewidth=1,
            clip_on=False,
            label=f"Average")
    ax.legend(loc="lower left")
    ax.set_xlabel("Position")
    ax.set_xlim(0, df.index.size)
    ax.set_xticks(np.arange(0, df.index.size + 1, 200))
    ax.set_ylim(0, 1)
    ax.set_ylabel("Pearson Correlation of Mutation Rates")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.subplots_adjust(bottom=0.15, top=0.95, left=0.13, right=0.93)
    fig.savefig(filename)
    plt.close()



def plot_k(df: pd.DataFrame,
           color: str,
           filename: str,
           max_k: int = 4):
    fig, ax = plt.subplots()
    fig.set_size_inches(3.6, 3.6)
    counts = np.zeros((max_k, df.index.size))
    for refnum in df.columns:
        col = df.loc[:, refnum]
        for i, k in enumerate(range(max_k, 0, -1)):
            counts[i] += (col >= k).astype(int)
    img = ax.imshow(counts / df.columns.size,
                    cmap=LinearSegmentedColormap.from_list("custom",
                                                           ["#ffffff", color]),
                    interpolation="none",
                    vmin=0,
                    vmax=1,
                    extent=[0, df.index.size, 0, max_k])
    plt.colorbar(img, label="Fraction of Simulations", orientation="horizontal", ax=ax, pad=0.2)
    ax.set_aspect(2/3 * df.index.size / max_k)
    ax.set_xlabel("Position")
    ax.set_xlim(0, df.index.size)
    ax.set_xticks(np.arange(0, df.index.size + 1, 200))
    ax.set_ylim(0, max_k)
    ax.set_yticks(np.arange(0, max_k + 1))
    ax.set_ylabel("Number of Clusters Detected")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # Plot the ground truth.
    for x1, x2, y in [(0, 200, 1),
                      (200, 600, 2),
                      (600, 650, 1),
                      (650, 800, 3),
                      (800, 900, 2),
                      (900, 1000, 2),
                      (1000, 1200, 1)]:
        ax.plot([x1 + 4, x2 - 4], [y, y],
                c="#000000",
                linewidth=1,
                zorder=10,
                solid_capstyle="butt")
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.13, right=0.93)
    fig.savefig(filename)
    plt.close()


REGIONS = [(200, 1),
           (400, 2),
           ( 50, 1),
           (150, 3),
           (100, 2),
           (100, 2),
           (200, 1)]


COLORS = {
    "seismic": "#0072b2",
    "draco": "#e69f00"
}


if __name__ == "__main__":
    num_refs = 60
    refnums = list(range(num_refs))
    # Accuracies of SEISMIC-RNA.
    print("Getting accuracies of SEISMIC-RNA ...")
    accuracies_seismic_file = Path("accuracies-seismic.csv")
    ks_seismic_file = Path("ks-seismic.csv")
    try:
        accuracies_seismic = pd.read_csv(accuracies_seismic_file, index_col=0)
        ks_seismic = pd.read_csv(ks_seismic_file, index_col=0)
    except FileNotFoundError:
        accuracies_seismic, ks_seismic = get_accuracies_seismic(1000000, refnums)
        accuracies_seismic.to_csv(accuracies_seismic_file)
        ks_seismic.to_csv("ks-seismic.csv")
    # Accuracies of DRACO.
    print("Getting accuracies of DRACO ...")
    accuracies_draco_file = Path("accuracies-draco.csv")
    ks_draco_file = Path("ks-draco.csv")
    try:
        accuracies_draco = pd.read_csv(accuracies_draco_file, index_col=0)
        ks_draco = pd.read_csv(ks_draco_file, index_col=0)
    except FileNotFoundError:
        accuracies_draco, ks_draco = get_accuracies_draco(1000000, refnums, 100)
        accuracies_draco.to_csv(accuracies_draco_file)
        ks_draco.to_csv(ks_draco_file)
    # Plot the results.
    plot_r(accuracies_seismic, COLORS["seismic"], "accuracies-seismic.svg")
    plot_k(ks_seismic, COLORS["seismic"], "ks-seismic.svg")
    plot_r(accuracies_draco, COLORS["draco"], "accuracies-draco.svg")
    plot_k(ks_draco, COLORS["draco"], "ks-draco.svg")
