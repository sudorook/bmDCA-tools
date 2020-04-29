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
        "-c", "--mcmc", dest="mcmc", action='append', required=True, help="mcmc file"
    )
    parser.add_argument(
        "-l", "--legend", dest="legend", action='append',  required=True, help="legends"
    )
    parser.add_argument(
        "-t", "--title", dest="title", required=True, help="title prefix"
    )
    parser.add_argument(
        "-o", "--output", dest="output", required=True, help="output file name"
    )
    parser.add_argument(
        "-m",
        "--offset_mode",
        dest="mode",
        default=0,
        help="offset type (default none)",
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
    mcmc_files = options.mcmc
    mcmc_energies = [load_energies(mcmc_file) for mcmc_file in mcmc_files]
    mcmc_labels = options.legend

    standard_energy = 0
    if options.mode == "1":
        standard_energy = msa_energies[0]
    elif options.mode == "2":
        standard_energy = np.mean(msa_energies)

    msa_energies = msa_energies - standard_energy
    mcmc_energies = [mcmc_energy - standard_energy for mcmc_energy in mcmc_energies]

    with plt.style.context("fivethirtyeight"):
        plt.hist(
            msa_energies,
            alpha=0.5,
            label="MSA",
            density=True,
        )
        for i, mcmc_energy in enumerate(mcmc_energies):
            plt.hist(
                mcmc_energy,
                alpha=0.5,
                label=mcmc_labels[i],
                density=True,
            )
        plt.legend(loc="upper right")
        plt.xlabel("Sequence Energy")
        plt.xlim(-270, -110)
        plt.ylabel("Probability")
        plt.title(options.title)
        plt.savefig(options.output)
        plt.close()
    
    #  with plt.style.context("fivethirtyeight"):
    #      plt.hist([msa_energies] + mcmc_energies,
    #          alpha=1,
    #          label=["MSA"] + mcmc_labels,
    #          density=True,
    #      )
    #      plt.legend(loc="upper right")
    #      plt.xlabel("Sequence Energy")
    #      plt.xlim(-270, -110)
    #      plt.ylabel("Probability")
    #      plt.title(options.title)
    #      plt.savefig(options.output)
    #      plt.close()


if __name__ == "__main__":
    main()
