#! /usr/bin/env python3
""" docstring """

import argparse
import os
import numpy as np
import pandas as pd
import statsmodels.formula.api as sm
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(rc={"figure.figsize": (8, 12)})


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="params",
        required=True,
        help="input parameters file",
    ),
    parser.add_argument(
        "-j",
        "--input2",
        dest="params2",
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
            h_new[i] += (h[i, j]) ** 2
        h_new[i] = h_new[i] ** 0.5
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
    prefix2 = os.path.splitext((options.params2).replace("/", "_"))[0]

    h1, J1 = load_data(options.params)
    h2, J2 = load_data(options.params2)

    J1_1d = J1.flatten()
    J2_1d = J2.flatten()
    h1_1d = h1.flatten()
    h2_1d = h2.flatten()

    df_h = pd.DataFrame(data={"h1": h1_1d, "h2": h2_1d,})
    res_h = sm.ols(formula="h2 ~ h1", data=df_h).fit()
    params_h = res_h.params

    df_J = pd.DataFrame(data={"J1": J1_1d, "J2": J2_1d,})
    res_J = sm.ols(formula="J2 ~ J1", data=df_J).fit()
    params_J = res_J.params

    with plt.style.context("fivethirtyeight"):
        fig, ax = plt.subplots(2, 1)

        df_h.plot.hexbin(
            ax=ax[0], x="h1", y="h2", mincnt=1, cmap="viridis", zorder=1
        )
        x = np.linspace(*(ax[0]).get_xlim())
        ax[0].plot(
            x,
            params_h.h1 * x + params_h.Intercept,
            label=r"y={:.3g}x+{:.3g}, $R^2$={:.3g}".format(
                params_h.h1, params_h.Intercept, res_h.rsquared
            ),
            alpha=0.5,
        )
        ax[0].plot(x, x, "--k", alpha=0.25, zorder=0, label=r"y=x")
        ax[0].set_xlabel("h value (old)")
        ax[0].set_ylabel("h value (new)")
        ax[0].legend(loc="lower right")

        df_J.plot.hexbin(
            ax=ax[1],
            x="J1",
            y="J2",
            bins="log",
            mincnt=1,
            cmap="viridis",
            zorder=1,
        )
        x = np.linspace(*ax[1].get_xlim())
        ax[1].plot(
            x,
            params_J.J1 * x + params_J.Intercept,
            label=r"y={:.3g}x+{:.3g}, $R^2$={:.3g}".format(
                params_J.J1, params_J.Intercept, res_J.rsquared
            ),
            alpha=0.5,
        )
        ax[1].plot(x, x, "--k", alpha=0.25, zorder=0, label=r"y=x")
        ax[1].set_xlabel("J value (old)")
        ax[1].set_ylabel("J value (new)")
        ax[1].legend(loc="lower right")

        fig.suptitle(prefix + " vs " + prefix2)

        plt.tight_layout()
        plt.savefig(prefix + "_" + prefix2 + ".png")
        plt.close()

    J1_1d = reduce_tensor(J1).flatten()
    J2_1d = reduce_tensor(J2).flatten()
    h1_1d = reduce_h(h1).flatten()
    h2_1d = reduce_h(h2).flatten()

    df_h = pd.DataFrame(data={"h1": h1_1d, "h2": h2_1d,})
    res_h = sm.ols(formula="h2 ~ h1", data=df_h).fit()
    params_h = res_h.params

    df_J = pd.DataFrame(data={"J1": J1_1d, "J2": J2_1d,})
    res_J = sm.ols(formula="J2 ~ J1", data=df_J).fit()
    params_J = res_J.params

    with plt.style.context("fivethirtyeight"):
        fig, ax = plt.subplots(2, 1)

        df_h.plot.hexbin(
            ax=ax[0], x="h1", y="h2", mincnt=1, cmap="viridis", zorder=1
        )
        x = np.linspace(*(ax[0]).get_xlim())
        ax[0].plot(
            x,
            params_h.h1 * x + params_h.Intercept,
            label=r"y={:.3g}x+{:.3g}, $R^2$={:.3g}".format(
                params_h.h1, params_h.Intercept, res_h.rsquared
            ),
            alpha=0.5,
        )
        ax[0].plot(x, x, "--k", alpha=0.25, zorder=0, label=r"y=x")
        ax[0].set_xlabel("h value (old)")
        ax[0].set_ylabel("h value (new)")
        ax[0].legend(loc="lower right")

        df_J.plot.hexbin(
            ax=ax[1],
            x="J1",
            y="J2",
            bins="log",
            mincnt=1,
            cmap="viridis",
            zorder=1,
        )
        x = np.linspace(*ax[1].get_xlim())
        ax[1].plot(
            x,
            params_J.J1 * x + params_J.Intercept,
            label=r"y={:.3g}x+{:.3g}, $R^2$={:.3g}".format(
                params_J.J1, params_J.Intercept, res_J.rsquared
            ),
            alpha=0.5,
        )
        ax[1].plot(x, x, "--k", alpha=0.25, zorder=0, label=r"y=x")
        ax[1].set_xlabel("J value (old)")
        ax[1].set_ylabel("J value (new)")
        ax[1].legend(loc="lower right")

        fig.suptitle(prefix + " vs " + prefix2)

        plt.tight_layout()
        plt.savefig(prefix + "_" + prefix2 + "_frob.png")
        plt.close()


if __name__ == "__main__":
    main()
