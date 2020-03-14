#! /usr/bin/env python3
""" docstring """

import argparse
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [12, 16]


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", dest="input", required=True, help="run log file",
    )
    parser.add_argument(
        "-o", "--output", dest="output", required=True, help="figure name",
    )
    return parser.parse_args()


def main():
    """ do stuff """
    options = parse_options()
    df = pd.read_csv(options.input, sep="\t")
    df["log-burn-in"] = np.log10(df["burn-in"])
    df["log-burn-between"] = np.log10(df["burn-between"])

    with plt.style.context("fivethirtyeight"):
        ax = df.plot(
            x="step",
            y=[
                "log-burn-between",
                "log-burn-in",
                "energy-err",
                "error-tot",
                "step-time",
            ],
            subplots=True,
        )
        plt.tight_layout()
        plt.savefig(options.output)


if __name__ == "__main__":
    main()
