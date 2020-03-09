#! /usr/bin/env python3
""" docstring """

import argparse
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [10, 7.5]


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--msa_input",
        dest="msa",
        required=True,
        help="numerical msa alignment",
    )
    parser.add_argument(
        "-c",
        "--mcmc_input",
        dest="mcmc",
        required=True,
        help="numerical mcmc alignment",
    )
    parser.add_argument(
        "-t", "--title", dest="title", required=True, help="title"
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="figure output name",
    )
    return parser.parse_args()


def main():
    """ do stuff """
    options = parse_options()

    data_df = pd.DataFrame(
        data={
            "MCMC": np.loadtxt(options.mcmc, dtype="float"),
            "MSA": np.loadtxt(options.msa, dtype="float"),
        }
    )
    #  data_df.loc[~(data_df==0).all(axis=1)]

    with plt.style.context("fivethirtyeight"):
        ax = data_df.plot.hexbin(
            x="MSA", y="MCMC", bins="log", mincnt=1, cmap="viridis", zorder=1
        )
        x = np.linspace(*ax.get_xlim())
        ax.plot(x, x, "--k", alpha=0.25, zorder=0)
        plt.title(options.title)
        plt.tight_layout()
        plt.savefig(options.output + ".png")
        plt.close()


if __name__ == "__main__":
    main()
