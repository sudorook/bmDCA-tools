#! /usr/bin/env python3
""" docstring """

import argparse
import numpy as np
import pandas as pd
import random
from scipy.sparse import csr_matrix as sparsify
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pysca import scaTools as sca

plt.rcParams["figure.figsize"] = [10, 7.5]
sns.set(rc={"figure.figsize": (10, 7.5)})


random.seed(0)
np.random.seed(0)


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
        "-r",
        "--mcmc_reg",
        dest="mcmc_reg",
        required=True,
        help="numerical mcmc alignment (reg)",
    )
    parser.add_argument(
        "-n",
        "--mcmc_new",
        dest="mcmc_new",
        required=True,
        help="numerical mcmc alignment (new)",
    )
    parser.add_argument(
        "-o",
        "--mcmc_old",
        dest="mcmc_old",
        required=True,
        help="numerical mcmc alignment (old)",
    )
    parser.add_argument(
        "-w", "--weights", dest="weights", required=True, help="msa weights"
    )
    return parser.parse_args()


def load_sequences(seq_file):
    """ load sequences """
    data = np.loadtxt(seq_file, dtype="float", skiprows=1)
    return data


def load_sequence_weights(weight_file):
    """ load msa weights """
    return np.loadtxt(weight_file, dtype="float")


def compute_sequence_weights(msa, Q=21):
    """ comptue seq weights """
    X2d = sca.alg2bin(msa, Q)
    simMat = (X2d.dot(X2d.T)).todense() / msa.shape[1]
    seqw = np.array(1 / (simMat > 0.8).sum(axis=0))
    #  return seqw[0] # ???
    return seqw  # ???


def compute_sim_mat(msa, Q=21):
    """ Compute sequence similarity matrix """
    [M, N] = msa.shape
    Abin_tensor = np.zeros((Q, N, M))
    for ia in range(Q):
        Abin_tensor[ia, :, :] = (msa == ia + 1).T
    #  Abin = sparsify(Abin_tensor.reshape(Q * N, M, order="F").T)
    #  return (Abin.dot(Abin.T)).todense() / N
    Abin = Abin_tensor.reshape(Q * N, M, order="F").T
    return Abin.dot(Abin.T) / N


def alg2bin(alg, N_aa=21):
    """
    Translate an alignment of matrix of size M sequences by L positions where
    the amino acids are represented by numbers between 0 and N_aa (obtained
    using lett2num) to a sparse binary array of size M x (N_aa x L).

    **Example**::

      Abin = alg2bin(alg, N_aa=20)
    """

    [N_seq, N_pos] = alg.shape
    Abin_tensor = np.zeros((N_aa, N_pos, N_seq))
    for ia in range(N_aa):
        Abin_tensor[ia, :, :] = (alg == ia + 1).T
    Abin = sparsify(Abin_tensor.reshape(N_aa * N_pos, N_seq, order="F").T)
    return Abin


def sample_msa_weighted(msa, weights, N):
    """ subset the alignment """
    Ntot = min(int(sum(weights)), N)
    wlist = [w for w in weights]
    selection2 = np.random.choice(
        range(msa.shape[0]), Ntot, replace=False, p=weights / sum(weights)
    )
    return msa[selection2, :]


def main():
    """ do stuff """
    options = parse_options()

    msa = load_sequences(options.msa) + 1
    mcmc_new = load_sequences(options.mcmc_new) + 1
    mcmc_old = load_sequences(options.mcmc_old) + 1
    mcmc_reg = load_sequences(options.mcmc_reg) + 1

    msa_weights = load_sequence_weights(options.weights)

    # too many sequences, so subset the alignment and recompute the weights
    msa_subset = sample_msa_weighted(msa, msa_weights, 25000)
    msa_subset_2 = sample_msa_weighted(msa, msa_weights, 10000)
    msa_subset_weights = compute_sequence_weights(msa_subset)

    # compute eignvectors for the msa
    msa_subset_sim = compute_sim_mat(msa_subset)
    posw, Dia, Di = sca.posWeights(
        msa_subset, msa_subset_weights, 0, 21, np.ones(21) / 21
    )
    X2d = sca.alg2bin(msa_subset, 21)
    X2d_2 = sca.alg2bin(msa_subset_2, 21)
    X2dw = sparsify(np.diag(np.sqrt(msa_subset_weights[0]))).dot(X2d)
    u, s, v = sca.svdss(X2dw, k=250)
    print(s)
    print((s[0] ** 2 + s[1] ** 2) / sum(s ** 2))
    tmp = X2d_2.dot(v).dot(np.diag(1 / s))
    data_df = pd.DataFrame(data={"PC1": tmp[:, 0], "PC2": tmp[:, 1]})

    # plot the projection of the msa onto the first two eigenvectors
    with plt.style.context("fivethirtyeight"):
        ax = data_df.plot.hexbin(
            x="PC1", y="PC2", mincnt=1, cmap="viridis", zorder=1
        )
        [left, right] = ax.get_xlim()
        [top, bottom] = ax.get_ylim()
        plt.title("Projection of MSA samples onto top eigenvectors")
        plt.tight_layout()
        plt.savefig("msa_proj.svg")
        plt.close()

    # project new mcmc onto msa eigenvectors
    X2d_new = sca.alg2bin(mcmc_new, 21)
    tmp_new = X2d_new.dot(v).dot(np.diag(1 / s))

    data_df_new = pd.DataFrame(
        data={"PC1": tmp_new[:, 0], "PC2": tmp_new[:, 1]}
    )

    with plt.style.context("fivethirtyeight"):
        ax = data_df_new.plot.hexbin(
            x="PC1", y="PC2", mincnt=1, cmap="viridis", zorder=1
        )
        plt.title("Projection of new MCMC onto MSA space")
        ax.set_xlim(left, right)
        ax.set_ylim(top, bottom)
        plt.tight_layout()
        plt.savefig("mcmc_new_proj.svg")
        plt.close()

    # project reg mcmc onto msa eigenvectors
    X2d_reg = sca.alg2bin(mcmc_reg, 21)
    tmp_reg = X2d_reg.dot(v).dot(np.diag(1 / s))

    data_df_reg = pd.DataFrame(
        data={"PC1": tmp_reg[:, 0], "PC2": tmp_reg[:, 1]}
    )

    with plt.style.context("fivethirtyeight"):
        ax = data_df_reg.plot.hexbin(
            x="PC1", y="PC2", mincnt=1, cmap="viridis", zorder=1
        )
        plt.title("Projection of reg MCMC onto MSA space")
        ax.set_xlim(left, right)
        ax.set_ylim(top, bottom)
        plt.tight_layout()
        plt.savefig("mcmc_reg_proj.svg")
        plt.close()

    # project old mcmc onto msa eigenvectors
    X2d_old = sca.alg2bin(mcmc_old, 21)
    tmp_old = X2d_old.dot(v).dot(np.diag(1 / s))

    data_df_old = pd.DataFrame(
        data={"PC1": tmp_old[:, 0], "PC2": tmp_old[:, 1]}
    )

    with plt.style.context("fivethirtyeight"):
        ax = data_df_old.plot.hexbin(
            x="PC1", y="PC2", mincnt=1, cmap="viridis", zorder=1
        )
        plt.title("Projection of old MCMC onto MSA space")
        ax.set_xlim(left, right)
        ax.set_ylim(top, bottom)
        plt.tight_layout()
        plt.savefig("mcmc_old_proj.svg")
        plt.close()


if __name__ == "__main__":
    main()
