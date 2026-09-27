import sys

import pandas as pd


def calc_pearson(file: str):
    data = pd.read_csv(file, header=[0, 1], index_col=[0, 1])
    data.dropna(axis=0, how="any", inplace=True)
    pearson = data.corr()
    if pearson.shape != (2, 2):
        raise ValueError(pearson)
    return float(pearson.iloc[0, 1])


if __name__ == "__main__":
    file = sys.argv[1]
    min_pearson = float(sys.argv[2])
    if calc_pearson(file) >= min_pearson:
        print(1)
    else:
        print(0)

