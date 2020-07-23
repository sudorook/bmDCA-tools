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

from pysca import scaTools as sca
import dcatools as tools

plt.rcParams["figure.figsize"] = [10, 7.5]

random.seed(0)
np.random.seed(0)


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n",
        "--msa_input",
        dest="msa_basis",
        required=True,
        help="numerical msa alignment for basis",
    )
    parser.add_argument(
        "-N",
        "--msa_proj",
        dest="msa_proj",
        required=True,
        help="numerical msa to project",
    )
    parser.add_argument(
        "-w", "--weights_input", dest="weights_input", help="msa basis weights"
    )
    parser.add_argument(
        "-W", "--weights_proj", dest="weights_proj", help="msa projection weights"
    )
    parser.add_argument(
        "-t", "--title", dest="title", help="figure title"
    )
    return parser.parse_args()


def load_sequences(seq_file):
    """ load sequences """
    data = np.loadtxt(seq_file, dtype="float", skiprows=1)
    return data


def load_sequence_weights_input(weight_file):
    """ load msa weights_input """
    return np.loadtxt(weight_file, dtype="float")


def compute_sequence_weights_input(msa, Q=21):
    """ comptue seq weights_input """
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


def sample_msa_weighted(msa, weights_input, N):
    """ subset the alignment """
    Ntot = min(int(sum(weights_input)), N)
    wlist = [w for w in weights_input]
    selection2 = np.random.choice(
        range(msa.shape[0]), Ntot, replace=False, p=weights_input / sum(weights_input)
    )
    return msa[selection2, :]


def main():
    """ do stuff """
    options = parse_options()

    msa_basis = load_sequences(options.msa_basis) + 1
    msa_proj = load_sequences(options.msa_proj) + 1

    outfile_basis = "msa_basis_" + tools.slugify(options.title) + ".png"
    outfile_proj = "msa_proj_" + tools.slugify(options.title) + ".png"

    [M, N] = msa_basis.shape
    if options.weights_input is not None:
        msa_weights_basis = load_sequence_weights_input(options.weights_input)
    else:
        msa_weights_basis = np.ones(M)
    
    [M2, N2] = msa_proj.shape
    if options.weights_proj is not None:
        msa_weights_proj = load_sequence_weights_input(options.weights_proj)
    else:
        msa_weights_proj = np.ones(M2)

    if M > 10000:
        # too many sequences, so subset the alignment and recompute the
        # weights_input
        print("subsetting the alignment...")
        msa_basis = sample_msa_weighted(msa_basis, msa_weights_basis, 10000)
        msa_weights_basis = compute_sequence_weights_input(msa_basis)

    # compute eignvectors for the msa
    msa_sim = compute_sim_mat(msa_basis)
    #  posw, Dia, Di = sca.posWeights(
    #      msa_basis, msa_weights_basis, 0, 21, np.ones(21) / 21
    #  )
    X2d = sca.alg2bin(msa_basis, 21)
    #  X2dw = sparsify(np.diag(np.sqrt(msa_weights_basis[0]))).dot(X2d)
    X2dw = sparsify(np.diag(np.sqrt(msa_weights_basis))).dot(X2d)
    u, s, v = sca.svdss(X2dw, k=250)

    #  print(s)
    #  print((s[0] ** 2 + s[1] ** 2) / sum(s ** 2))

    tmp = X2d.dot(v).dot(np.diag(1 / s))
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
        plt.savefig(outfile_basis)
        plt.close()

    # project new mcmc onto msa eigenvectors
    X2d_new = sca.alg2bin(msa_proj, 21)
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
        plt.savefig(outfile_proj)
        plt.close()


if __name__ == "__main__":
    main()
