"""
Count the average mutations per read using the histread files.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def process_histread_file(file: Path):
    counts = pd.read_csv(file, index_col="Count")["Mutated"]
    assert isinstance(counts, pd.Series)
    total_muts = np.vdot(counts.index, counts)
    total_reads = np.sum(counts)
    return float(total_muts / total_reads)


def count_avg_muts_per_read_pattern(pattern: str):
    print(pattern)
    files = list(Path().glob(pattern))
    if not files:
        print(f"No files found matching pattern: {pattern}")
    avgs_muts_per_read = list()
    for file in files:
        avgs_muts_per_read.append(process_histread_file(file))
    avgs_muts_per_read = np.array(avgs_muts_per_read)
    print(f"Number of files: {len(files)}")
    print(f"Mean: {np.mean(avgs_muts_per_read)}")
    print(f"Stdev: {np.std(avgs_muts_per_read)}")
    plt.hist(avgs_muts_per_read,
             bins=np.linspace(np.floor(np.min(avgs_muts_per_read)),
                              np.ceil(np.max(avgs_muts_per_read)),
                              len(files) + 1))
    plt.title(pattern)
    plt.xlabel("Average mutations per read")
    plt.ylabel("Frequency")
    plt.show()
    return avgs_muts_per_read
    
def count_avg_muts_per_read(patterns: list[str]):
    for pattern in patterns:
        count_avg_muts_per_read_pattern(pattern)


if __name__ == "__main__":
    patterns = [
        "out/sample-biased-200000/graph/rna-140-*/*/histread_filtered_m-count.csv",
        "out/sample-biased-200000/graph/rna-280-*/*/histread_filtered_m-count.csv",
        "out/sample-biased-200000/graph/rna-420-*/*/histread_filtered_m-count.csv",
    ]
    count_avg_muts_per_read(patterns)

