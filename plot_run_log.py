#! /usr/bin/env python3
""" docstring """

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

#  import seaborn as sns

plt.rcParams["figure.figsize"] = [10, 16]
#  sns.set(rc={"figure.figsize": (10, 7.5)})


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", dest="input", required=True, help="run log file",
    )
    return parser.parse_args()


def main():
    """ do stuff """
    options = parse_options()
    df = pd.read_csv(options.input, sep="\t")
    print(df)

    with plt.style.context("fivethirtyeight"):
        ax = df.plot(
            x="step",
            y=["burn-between", "burn-in", "energy-err", "error-tot"],
            subplots=True,
        )
        plt.tight_layout()
        plt.savefig("figure.png")


if __name__ == "__main__":
    main()
