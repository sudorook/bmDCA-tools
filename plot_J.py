#! /usr/bin/env python3
""" docstring """

import argparse
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(rc={"figure.figsize": (10, 7.5)})


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="params",
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


def reduce_h(h):
    """ reduce rank 4 tensor using 2-norm """
    Npos = np.shape(h)[0]
    Naa = np.shape(h)[1]
    h_new = np.zeros(Npos)
    for i in range(Npos):
        for j in range(Naa):
            h_new[i] += (h[i, j])**2
        h_new[i] = h_new[i]**.5
    return h_new


def load_data(data_file):
    """ syntax checker is annoying me """

    tmpJ = []
    tmph = []
    Naa = 0
    Npos = 0
    with open(data_file, "r") as handle:
        for line in handle:
            if line[0] == "J":
                tmpJ.append(float(line.split(" ")[5].strip()))
            elif line[0] == "h":
                tmp = line.split(" ")
                tmph.append(float(line.split(" ")[3].strip()))
                if tmp[1] == "0":
                    Naa += 1
                if tmp[2] == "0":
                    Npos += 1

    tmpJ = np.array(tmpJ)
    tmph = np.array(tmph)

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

    prefix = os.path.splitext((options.params).replace("/", "_"))[0]

    h, J = load_data(options.params)
    J_2d = reduce_tensor(J, "fro")
    h_1d = reduce_h(h)
    model = J_2d + np.transpose(J_2d) + np.diag(h_1d)

    # Make plots
    ax = sns.heatmap(model, cmap="Blues", linewidths=0.25, vmin=0, vmax=2)
    fig = ax.get_figure()
    ax.xaxis.set_ticks_position("top")
    ax.yaxis.set_ticks_position("left")
    ax.set(xlabel="Position", ylabel="Position", title="J Parameters")
    plt.tight_layout()
    fig.savefig(prefix + "_heatmap.svg")
    plt.close()


if __name__ == "__main__":
    main()
