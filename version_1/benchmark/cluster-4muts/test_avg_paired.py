import os
from concurrent import futures

import numpy as np
from matplotlib import pyplot as plt

from seismicrna.core.seq import RNA


def simulate(reflen: int, number: int):
    rna = RNA.random(reflen)
    filename = f"rna-{reflen}-{number}"
    with open(f"{filename}.fa", "w") as f:
        f.write(f">Sequence\n{rna}")
    os.system(f"Fold --MFE {filename}.fa {filename}.ct")
    os.system(f"ct2dot {filename}.ct 0 {filename}.dot")
    with open(f"{filename}.dot", "r") as f:
        f.readline()
        f.readline()
        structure = f.readline().strip()
    os.remove(f"{filename}.fa")
    os.remove(f"{filename}.ct")
    os.remove(f"{filename}.dot")
    f_paired = sum(x != "." for x in structure) / len(structure)
    print(f"RESULT: {filename} {f_paired}")
    return f_paired


def run_simulations(reflen: int, num_sims: int):
    with futures.ProcessPoolExecutor() as executor:
        f_paireds = [
            future.result() for future in futures.as_completed(
                [executor.submit(simulate, reflen, number)
                 for number in range(num_sims)]
            )
        ]
    return f_paireds


def plot_results(reflen: int, num_sims: int):
    f_paireds = run_simulations(reflen, num_sims)
    print("Average fraction of paired bases:", np.mean(f_paireds))
    plt.hist(f_paireds, bins=20)
    plt.show()


if __name__ == "__main__":
    plot_results(140, 1000)  # 0.5607
    plot_results(280, 1000)  # 0.5889714285714285
    plot_results(420, 1000)  # 0.6014857142857142
