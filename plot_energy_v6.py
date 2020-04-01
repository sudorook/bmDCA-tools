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
        "-c", "--mcmc", dest="mcmc", required=True, help="mcmc file"
    )
    parser.add_argument(
        "-c2", "--mcmc2", dest="mcmc2", required=True, help="mcmc file 2"
    )
    parser.add_argument(
        "-t", "--title", dest="title", required=True, help="title prefix"
    )
    parser.add_argument(
        "-o", "--output", dest="output", required=True, help="output file name"
    )
    return parser.parse_args()


def load_energies(energy_file):
    """ load sequences """
    data = np.loadtxt(energy_file, dtype="double", skiprows=1)
    return data


def main():
    """ do stuff """
    options = parse_options()

    msa_energies = load_energies(options.msa)
    mcmc_energies = load_energies(options.mcmc)
    mcmc_energies2 = load_energies(options.mcmc2)

    standard_energy = msa_energies[0]  # e coli energy
    #  standard_energy = np.mean(msa_energies)

    with plt.style.context("fivethirtyeight"):
        plt.hist(
            msa_energies - standard_energy,
            alpha=0.5,
            label="MSA",
            density=True,
        )
        plt.hist(
            mcmc_energies - standard_energy,
            alpha=0.5,
            label="MCMC",
            density=True,
        )
        plt.hist(
            mcmc_energies2 - standard_energy,
            alpha=0.5,
            label="MCMC (Martin)",
            density=True,
        )
        plt.legend(loc="upper right")
        plt.xlabel("Sequence Energy")
        #  plt.xlim(-100, 150)
        plt.ylabel("Probability")
        #  plt.title("Normed histogram of sequence energies (new bmDCA)")
        plt.title(options.title)
        plt.savefig(options.output)
        plt.close()


if __name__ == "__main__":
    main()
