#! /usr/bin/env python3
""" docstring """

import argparse
import os
from itertools import combinations
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [10, 7.5]


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--msa", dest="msa", required=True, help="msa file"
    )
    parser.add_argument(
        "-S",
        "--msa_label",
        dest="msa_label",
        required=True,
        help="msa label for plot",
    )
    parser.add_argument(
        "-c", "--mcmc", dest="mcmc", required=True, help="mcmc file"
    )
    parser.add_argument(
        "-C",
        "--mcmc_label",
        dest="mcmc_label",
        required=True,
        help="mcmc label for plot",
    )
    parser.add_argument(
        "-o", "--output", dest="output", required=True, help="output file name"
    )
    parser.add_argument(
        "-t", "--title", dest="title", required=True, help="plot title"
    )
    return parser.parse_args()


def load_distances(energy_file):
    """ load sequences """
    data = np.loadtxt(energy_file, dtype="double", skiprows=1)
    return data


def main():
    """ do stuff """
    options = parse_options()

    msa_distances = load_distances(options.msa)
    mcmc_distances = load_distances(options.mcmc)

    #  tmp = sum(msa_distances > .8)
    #  print(tmp)
    #  tmp = sum(mcmc_distances > .8)
    #  print(tmp)

    with plt.style.context("fivethirtyeight"):
        plt.hist(
            msa_distances, alpha=0.5, label=options.msa_label, density=True,
        )
        plt.hist(
            mcmc_distances, alpha=0.5, label=options.mcmc_label, density=True,
        )
        plt.legend(loc="upper right")
        plt.xlabel("Sequence Hamming Distance")
        plt.xlim(0, 1)
        plt.ylabel("Probability")
        plt.xlabel("Min Pairwise Distance")
        plt.title(options.title)
        plt.savefig(options.output)
        plt.close()


if __name__ == "__main__":
    main()
