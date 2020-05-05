#! /usr/bin/env python3
""" docstring """

import argparse
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import dcatools as tools

plt.rcParams["figure.figsize"] = [10, 7.5]


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--mcmc", dest="mcmc", required=True, help="mcmc file"
    )
    parser.add_argument(
        "-p", "--params", dest="params", required=True, help="parameters file"
    )
    parser.add_argument(
        "-o", "--output", dest="output", required=True, help="output file name"
    )
    return parser.parse_args()


def main():
    """ do stuff """
    options = parse_options()

    h, J = tools.load_params(options.params)

    mcmc_seqs = tools.load_sequences(options.mcmc)
    mcmc_energies = tools.compute_energies(mcmc_seqs, h, J)

    np.savetxt(options.output, mcmc_energies)


if __name__ == "__main__":
    main()
