#! /usr/bin/env python3
""" docstring """

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(rc={"figure.figsize": (10, 7.5)})


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--old",
        dest="params_old",
        required=True,
        help="input parameters file",
    )
    parser.add_argument(
        "-n",
        "--new",
        dest="params_new",
        required=True,
        help="input parameters file",
    )
    return parser.parse_args()


def reduce_tensor(J, norm=2):
    """ reduce rank 4 tensor using 2-norm """
    Npos = np.shape(J)[0]
    J_new = np.zeros((Npos, Npos))
    for i in range(Npos):
        for j in range(i + 1, Npos):
            J_new[i, j] = np.linalg.norm(J[i, j], norm)
    return J_new


def load_data(data_file):
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


def main():
    """ do stuff """
    options = parse_options()

    prefix_old = os.path.splitext((options.params_old).replace("/", "_"))[0]
    prefix_new = os.path.splitext((options.params_new).replace("/", "_"))[0]

    h_old, J_old = load_data(options.params_old)
    h_1d_old = h_old.flatten()
    J_1d_old = J_old.flatten()

    h_new, J_new = load_data(options.params_new)
    h_1d_new = h_new.flatten()
    J_1d_new = J_new.flatten()

    h_abs_error = abs(h_1d_old - h_1d_new)
    J_abs_error = abs(J_1d_old - J_1d_new)

    h_df = pd.DataFrame(data={"Old": h_old.flatten(), "New": h_new.flatten(),})
    J_df = pd.DataFrame(data={"Old": J_old.flatten(), "New": J_new.flatten(),})

    with sns.axes_style("whitegrid"):
        ax = sns.distplot(h_abs_error, kde=False)
        plt.title("h absolute error")
        plt.tight_layout()
        plt.savefig("h_abs_error.svg")
        plt.close()

    with sns.axes_style("whitegrid"):
        ax = sns.distplot(J_abs_error, kde=False)
        plt.title("J absolute error")
        plt.tight_layout()
        plt.savefig("J_abs_error.svg")
        plt.close()

    with plt.style.context("fivethirtyeight"):
        ax = J_df.plot.hexbin(
            x="Old", y="New", bins="log", mincnt=1, cmap="viridis", zorder=1
        )
        x = np.linspace(*ax.get_xlim())
        ax.plot(x, x, "--k", alpha=0.25, zorder=0)
        plt.title("J (old vs new)")
        plt.tight_layout()
        plt.savefig("J_hexbin.svg")
        plt.close()

    with plt.style.context("fivethirtyeight"):
        ax = h_df.plot.hexbin(
            x="Old", y="New", mincnt=1, cmap="viridis", zorder=1
        )
        x = np.linspace(*ax.get_xlim())
        ax.plot(x, x, "--k", alpha=0.25, zorder=0)
        plt.title("h (old vs new)")
        plt.tight_layout()
        plt.savefig("h_hexbin.svg")
        plt.close()


if __name__ == "__main__":
    main()
