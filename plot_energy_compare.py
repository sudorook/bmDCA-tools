#! /usr/bin/env python3
""" docstring """

import argparse
import os
from itertools import combinations
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import dcatools as tools

plt.rcParams["figure.figsize"] = [10, 12]


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--msa", dest="msa", required=True, help="msa file"
    )
    parser.add_argument(
        "-c", "--mcmc", dest="mcmc", required=True, help="mcmc file"
    )
    parser.add_argument(
        "-S",
        "--msa_label",
        dest="msa_label",
        required=True,
        help="msa plot label",
    )
    parser.add_argument(
        "-C",
        "--mcmc_label",
        dest="mcmc_label",
        required=True,
        help="mcmc plot label",
    )
    #  parser.add_argument(
    #      "-t", "--title", dest="title", required=True, help="title prefix"
    #  )
    #  parser.add_argument(
    #      "-m",
    #      "--offset_mode",
    #      dest="mode",
    #      default=0,
    #      help="offset type (default none)",
    #  )
    #  parser.add_argument(
    #      "-o", "--output", dest="output", required=True, help="output file name"
    #  )
    return parser.parse_args()


def main():
    """ do stuff """
    options = parse_options()

    prefix1 = os.path.splitext((options.msa).replace("/", "_"))[0]
    prefix2 = os.path.splitext((options.mcmc).replace("/", "_"))[0]

    energies1 = tools.load_energies(options.msa)
    energies2 = tools.load_energies(options.mcmc)

    energies1 = energies1 - energies1[0]
    energies2 = energies2 - energies2[0]

    energies1_sorted = np.sort(energies1)
    energies2_sorted = np.sort(energies2)

    with plt.style.context("fivethirtyeight"):
        fig, ax = plt.subplots(2, 1)

        ax[0].scatter(x=energies1, y=energies2)
        x = np.linspace(*(ax[0]).get_xlim())
        ax[0].plot(x, x, "--k", alpha=0.25, zorder=0)
        ax[0].set_xlabel(options.msa_label)
        ax[0].set_ylabel(options.mcmc_label)

        ax[1].scatter(x=energies1_sorted, y=energies2_sorted)
        x = np.linspace(*ax[1].get_xlim())
        ax[1].plot(x, x, "--k", alpha=0.25, zorder=0)
        ax[1].set_xlabel(options.msa_label + " (sorted)")
        ax[1].set_ylabel(options.mcmc_label + " (sorted)")

        fig.suptitle(options.msa_label + " vs " + options.mcmc_label)

        plt.tight_layout()
        print(prefix1)
        plt.savefig(prefix1 + "_" + prefix2 + "_energies.png")
        plt.close()


if __name__ == "__main__":
    main()
