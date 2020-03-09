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
        "-o",
        "--mcmc_old",
        dest="mcmc_old",
        required=True,
        help="mcmc_old file",
    )
    parser.add_argument(
        "-n",
        "--mcmc_new",
        dest="mcmc_new",
        required=True,
        help="mcmc_new file",
    )
    parser.add_argument(
        "-p",
        "--params_new",
        dest="params_new",
        required=True,
        help="parameters file",
    )
    parser.add_argument(
        "-t",
        "--params_old",
        dest="params_old",
        required=True,
        help="params file",
    )
    #  parser.add_argument("-t", "--title", dest="title", required=True, help="title prefix")
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
        E = -sum(h[[list(range(len(seq))), seq]])
        #  E = sum(h[[list(range(len(seq))), seq]])
        tmp = list(combinations(list(range(len(seq))), 2))
        tmp2 = [(seq[x[0]], seq[x[1]]) for x in tmp]
        idx = np.c_[np.asarray(tmp), np.asarray(tmp2)]
        E = E - sum(J[idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]])
        #  E = E + sum(J[idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]])
        energies[i] = E
    return energies


def main():
    """ do stuff """
    options = parse_options()

    h_old, J_old = load_params(options.params_new)
    h_new, J_new = load_params(options.params_old)

    msa_seqs = load_sequences(options.msa)
    mcmc_old_seqs = load_sequences(options.mcmc_old)
    mcmc_new_seqs = load_sequences(options.mcmc_new)

    msa_old_energies = compute_energies(msa_seqs, h_old, J_old)
    msa_new_energies = compute_energies(msa_seqs, h_new, J_new)
    mcmc_old_energies = compute_energies(mcmc_old_seqs, h_old, J_old)
    mcmc_new_energies = compute_energies(mcmc_new_seqs, h_new, J_new)

    mcmc_old_energies_new = compute_energies(mcmc_old_seqs, h_new, J_new)
    mcmc_new_energies_old = compute_energies(mcmc_new_seqs, h_old, J_old)

    with plt.style.context("fivethirtyeight"):
        plt.hist(msa_new_energies, alpha=0.5, label="MSA", density=True)
        plt.hist(
            mcmc_new_energies, alpha=0.5, label="MCMC (New)", density=True
        )
        plt.legend(loc="upper right")
        plt.xlabel("Sequence Energy")
        plt.ylabel("Probability")
        plt.title("Normed histogram of sequence energies (new bmDCA)")
        plt.savefig("energies_new_hist_normed.svg")
        plt.close()

        plt.hist(msa_old_energies, alpha=0.5, label="MSA", density=True)
        plt.hist(
            mcmc_old_energies, alpha=0.5, label="MCMC (Old)", density=True
        )
        plt.legend(loc="upper right")
        plt.xlabel("Sequence Energy")
        plt.ylabel("Probability")
        plt.title("Normed histogram of sequence energies (old bmDCA)")
        plt.savefig("energies_old_hist_normed.svg")
        plt.close()

        plt.hist(msa_new_energies, alpha=0.5, label="MSA", density=True)
        plt.hist(
            mcmc_old_energies_new, alpha=0.5, label="MCMC (Old)", density=True
        )
        plt.hist(
            mcmc_new_energies, alpha=0.5, label="MCMC (New)", density=True
        )
        plt.legend(loc="upper right")
        plt.xlabel("Sequence Energy")
        plt.ylabel("Probability")
        plt.title("Combined normed hist of seq energies (new bmDCA params)")
        plt.savefig("energies_combined_new_hist_normed.svg")
        plt.close()

        plt.hist(msa_old_energies, alpha=0.5, label="MSA", density=True)
        plt.hist(
            mcmc_old_energies, alpha=0.5, label="MCMC (Old)", density=True
        )
        plt.hist(
            mcmc_new_energies_old, alpha=0.5, label="MCMC (New)", density=True
        )
        plt.legend(loc="upper right")
        plt.xlabel("Sequence Energy")
        plt.ylabel("Probability")
        plt.title("Combined normed hist of seq energies (old bmDCA params)")
        plt.savefig("energies_combined_old_hist_normed.svg")
        plt.close()


if __name__ == "__main__":
    main()
