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
    df = pd.read_csv(options.input, sep=",")
    df.columns = [
        "step",
        "step_importance",
        "burn-in",
        "burn-between",
        "auto-corr",
        "ergo-corr",
        "cross-corr",
        "energy-start",
        "energy-end",
        "ergo-err",
        "auto-err",
        "energy-err",
        "error-tot",
        "error-h",
        "error-j",
        "DELTA_T",
        "T_WAIT",
    ]

    df["log10-burn-in"] = np.log10(df["burn-in"])
    df["log10-burn-between"] = np.log10(df["burn-between"])

    with plt.style.context("fivethirtyeight"):
        ax = df.plot(
            x="step",
            y=[
                "log10-burn-between",
                "log10-burn-in",
                "energy-err",
                "error-tot",
            ],
            linewidth=1.0,
            subplots=True,
        )
        plt.tight_layout()
        plt.savefig(options.output)


if __name__ == "__main__":
    main()
