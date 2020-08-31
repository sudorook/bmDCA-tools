#! /usr/bin/env python3
""" docstring """

import argparse
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import dcatools as tools

plt.rcParams["figure.figsize"] = [12, 9]
plt.rcParams.update({"mathtext.default": "regular"})


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n",
        "--msa",
        dest="msa",
        action="append",
        required=False,
        help="numerical msa file(s)",
    )
    parser.add_argument(
        "-e",
        "--energies",
        dest="energies",
        action="append",
        required=False,
        help="msa energy file(s)",
    )
    parser.add_argument(
        "-p",
        "--parameters",
        dest="params",
        required=False,
        help="msa energy file(s)",
    )
    parser.add_argument(
        "-l",
        "--labels",
        dest="labels",
        action="append",
        required=True,
        help="labels",
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


def main():
    """ do stuff """
    options = parse_options()

    if options.energies is not None:
        msa_energies = [
            tools.load_energies(msa_file) for msa_file in options.energies
        ]
    elif options.msa is not None and options.params is not None:
        h, J = tools.load_model(options.params)
        msa_energies = [
            tools.compute_energies(tools.load_sequences(msa_file), h, J)
            for msa_file in options.msa
        ]
    else:
        sys.exit("ERROR: missing input data")

    msa_labels = options.labels

    standard_energy = 0
    if options.mode == "1":
        standard_energy = msa_energies[0][0]
        msa_energies = [
            msa_energy - standard_energy for msa_energy in msa_energies
        ]
    elif options.mode == "2":
        standard_energy = np.mean(msa_energies[0])
        msa_energies = [
            msa_energy - standard_energy for msa_energy in msa_energies
        ]

    with plt.style.context("fivethirtyeight"):
        for i, msa_energy in enumerate(msa_energies):
            plt.hist(msa_energy, alpha=0.5, label=msa_labels[i], density=True)
        plt.legend(loc="upper right")
        plt.xlabel("Sequence Energy")
        #  plt.xlim(-270, -110)
        plt.ylabel("Probability")
        plt.title(options.title)
        plt.tight_layout()
        plt.savefig(options.output)
        plt.close()


if __name__ == "__main__":
    main()
