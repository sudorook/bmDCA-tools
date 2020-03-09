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
        "-p", "--params", dest="params", required=True, help="parameters file"
    )
    parser.add_argument(
        "-t", "--title", dest="title", required=True, help="title prefix"
    )
    parser.add_argument(
        "-o", "--output", dest="output", required=True, help="output file name"
    )
    return parser.parse_args()


def load_params(data_file):
    """ syntax checker is annoying me """

    tmpJ = []
    tmph = []
    with open(data_file, "r") as handle:
        for line in handle:
            if line[0] == "J":
                tmpJ.append(float(line.split(" ")[5].strip()))
            elif line[0] == "h":
                tmph.append(float(line.split(" ")[3].strip()))

    tmpJ = np.array(tmpJ)
    tmph = np.array(tmph)

    #  h = np.loadtxt(h_file)
    Naa = 21
    Npos = int(len(tmph) / Naa)
    h = tmph.reshape(Npos, Naa)

    # Jijs are saved as upper triangular matrix...
    J = np.zeros((Npos, Npos, Naa, Naa))
    Ntriu = int((Npos * Npos - Npos) / 2)

    #  tmp = np.loadtxt(J_file)
    tmpJ = tmpJ.reshape(Ntriu, Naa, Naa)
    idx = np.triu_indices(Npos, 1)
    J[idx] = tmpJ

    return h, J


def load_sequences(msa_file):
    """ load sequences """
    data = np.loadtxt(msa_file, dtype="int", skiprows=1)
    return data


def compute_energies(seqs, h, J):
    energies = np.zeros(seqs.shape[0])
    for i, seq in enumerate(seqs):
        N = len(seq)
        for n1 in range(0, N):
            energies[i] += -h[n1][seq[n1]]
            for n2 in range(n1 + 1, N):
                energies[i] += -J[n1][n2][seq[n1]][seq[n2]]
    return energies


def main():
    """ do stuff """
    options = parse_options()

    h, J = load_params(options.params)

    msa_seqs = load_sequences(options.msa)
    mcmc_seqs = load_sequences(options.mcmc)

    msa_energies = compute_energies(msa_seqs, h, J)
    mcmc_energies = compute_energies(mcmc_seqs, h, J)

    standard_energy = msa_energies[0]  # e coli energy

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
        plt.legend(loc="upper right")
        plt.xlabel("Sequence Energy")
        plt.xlim(-50, 200)
        plt.ylabel("Probability")
        plt.title(options.title)
        plt.savefig(options.output)
        plt.close()


if __name__ == "__main__":
    main()
