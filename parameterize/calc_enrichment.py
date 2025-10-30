#!python

import glob
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from seismicrna.core.batch import RegionMutsBatch
from seismicrna.core.dataset import RegionDataset
from seismicrna.core.rel import (NOCOV, MATCH, DELET, INSRT,
                                 SUB_A, SUB_C, SUB_G, SUB_T)
from seismicrna.core.task import as_list_of_tuples, dispatch
from seismicrna.mask.dataset import load_mask_dataset


MUTATION_DECODE = {SUB_A: "A",
                   SUB_C: "C",
                   SUB_G: "G",
                   SUB_T: "T",
                   DELET: "-",
                   INSRT: "+"}
MUTATION_TYPES = set(MUTATION_DECODE)
INFO_TYPES = {MATCH} | MUTATION_TYPES
COLUMNS = ["Focal Base", "Focal Mutation", "Second Base", "Distance"]


def calc_matrix(dataset: RegionDataset):
    return pd.concat([batch.matrix for batch in dataset.iter_batches()])


def calc_p_mut(matrix: pd.DataFrame, min_info: int = 1000):
    n_muts = pd.Series(np.count_nonzero(matrix.isin(MUTATION_TYPES), axis=0),
                       matrix.columns)
    n_matches = pd.Series(np.count_nonzero(matrix == MATCH, axis=0),
                          matrix.columns)
    n_info = n_muts + n_matches
    n_info.loc[n_info < min_info] = np.nan
    return n_muts / n_info


def calc_odds(p):
    return p / (1. - p)


def calc_odds_ratios_in_matrix(matrix: pd.DataFrame, max_dist: int):
    keys = list()
    ratios = list()
    for pos_j, base_j in matrix.columns:
        matrix_slice = matrix.loc[:, pos_j - max_dist: pos_j + max_dist]
        for mut_j_code in MUTATION_TYPES:
            mut_j_base = MUTATION_DECODE[mut_j_code]
            if mut_j_base == base_j:
                continue
            is_mut_j = matrix[pos_j, base_j] == mut_j_code
            p_mut_given_mut_j = calc_p_mut(matrix_slice.loc[is_mut_j])
            is_not_mut_j = matrix[pos_j, base_j].isin(INFO_TYPES - {mut_j_code})
            p_mut_given_not_mut_j = calc_p_mut(matrix_slice.loc[is_not_mut_j])
            odds_ratios = (calc_odds(p_mut_given_mut_j)
                           / calc_odds(p_mut_given_not_mut_j))
            for (pos_i, base_i), ratio in odds_ratios.items():
                if np.isfinite(ratio):
                    distance = pos_i - pos_j
                    keys.append((base_j, mut_j_base, base_i, distance))
                    ratios.append(ratio)
    return pd.Series(ratios,
                     index=pd.MultiIndex.from_tuples(keys, names=COLUMNS),
                     dtype=float).reset_index(name="Enrichment")


def calc_odds_ratios_in_dataset(dataset: RegionDataset, max_dist: int):
    print(dataset)
    return calc_odds_ratios_in_matrix(calc_matrix(dataset), max_dist)


def calc_mutation_enrichment(inputs: str, max_dist: int):
    """ Probability that position i is mutated given that position j
    (a distance d = j - i away) is mutated, divided by the probability
    that i is mutated regardless of whether j is mutated. """ 
    datasets = as_list_of_tuples(load_mask_dataset.iterate(glob.iglob(inputs)))
    return pd.concat(dispatch(calc_odds_ratios_in_dataset,
                              num_cpus=os.cpu_count(),
                              pass_num_cpus=False,
                              as_list=True,
                              ordered=False,
                              raise_on_error=True,
                              args=datasets,
                              kwargs=dict(max_dist=max_dist)))


def plot_enrichment(enrichment: pd.DataFrame, output_plot: str):
    
    # Set up the figure with 2x4 subplots
    fig, axes = plt.subplots(2, 4, figsize=(20, 10), sharex=True, sharey=True)
    
    # Define the plot configurations
    focal_bases = ['A', 'C']
    focal_mutations = {
        'A': ['G', 'C', 'T', '-'],
        'C': ['T', 'A', 'G', '-'],
    }
    
    xmin = enrichment["Distance"].min() - 0.5
    xmax = enrichment["Distance"].max() + 0.5

    log_enrichment = np.log10(enrichment["Enrichment"])
    enrichment["Log Enrichment"] = log_enrichment
    finite_log_enrichment = log_enrichment[np.isfinite(log_enrichment)]
    ymin = np.min(finite_log_enrichment)
    ymax = np.max(finite_log_enrichment)
    

    # Plot each subplot
    for row, focal_base in enumerate(focal_bases):
        for col, focal_mut in enumerate(focal_mutations[focal_base]):
            ax = axes[row, col]
            
            # Plot each combination of 5' base and 5' mutation
            mask = ((enrichment["Focal Base"] == focal_base) &
                    (enrichment["Focal Mutation"] == focal_mut))
            sns.violinplot(enrichment.loc[mask], 
                           x="Distance", 
                           y="Log Enrichment",
                           hue="Second Base",
                           hue_order=["A", "C"],
                           inner=None,
                           palette={"A": "#D44D5C",
                                    "C": "#046E8F"},
                           linewidth=0.1,
                           ax=ax)
            ax.set_title(f"{focal_base} → {focal_mut}")
            ax.grid(True, axis="y")
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.set_xlim(-1, xmax - xmin - 1)
            ax.set_ylim(ymin, ymax)
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()


def run(inputs: str, output_csv: str, output_plot: str, max_dist: int):
    if os.path.isfile(output_csv):
        enrichment = pd.read_csv(output_csv)
    else:
        enrichment = calc_mutation_enrichment(inputs, max_dist)
        enrichment.to_csv(output_csv, index=False)
    fig = plot_enrichment(enrichment, output_plot)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))

