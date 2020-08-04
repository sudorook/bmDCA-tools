#! /usr/bin/env python3
""" docstring """

import argparse
import os
import numpy as np
import matplotlib
import pandas as pd
import statsmodels.formula.api as sm

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import dcatools as tools

plt.rcParams["figure.figsize"] = [10, 8]


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
        "-S",
        "--msa_label",
        dest="msa_label",
        required=True,
        help="msa plot label",
    )
    parser.add_argument(
        "-C",
        "--mcmc_label",
        dest="mcmc_label",
        required=True,
        help="mcmc plot label",
    )
    #  parser.add_argument(
    #      "-t", "--title", dest="title", required=True, help="title prefix"
    #  )
    #  parser.add_argument(
    #      "-m",
    #      "--offset_mode",
    #      dest="mode",
    #      default=0,
    #      help="offset type (default none)",
    #  )
    #  parser.add_argument(
    #      "-o", "--output", dest="output", required=True, help="output file name"
    #  )
    return parser.parse_args()


def main():
    """ do stuff """
    options = parse_options()

    prefix1 = os.path.splitext((options.msa).replace("/", "_"))[0]
    prefix2 = os.path.splitext((options.mcmc).replace("/", "_"))[0]

    energies1 = tools.load_energies(options.msa)
    energies2 = tools.load_energies(options.mcmc)

    #  energies1_sorted = np.sort(energies1)
    #  energies2_sorted = np.sort(energies2)

    df_e = pd.DataFrame(data={"e1": energies1, "e2": energies2,})
    res_e = sm.ols(formula="e2 ~ e1", data=df_e).fit()
    params_e = res_e.params

    #  df_es = pd.DataFrame(
    #      data={"e1": energies1_sorted, "e2": energies2_sorted,}
    #  )
    #  res_es = sm.ols(formula="e2 ~ e1", data=df_es).fit()
    #  params_es = res_es.params

    with plt.style.context("fivethirtyeight"):
        #  fig, ax = plt.subplots(2, 1)
        fig, ax = plt.subplots(1, 1)

        ax.scatter(x=energies1, y=energies2, color="midnightblue", s=5)
        x = np.linspace(*(ax).get_xlim())

        ax.plot(x, x, "--k", alpha=0.25, zorder=0)
        ax.plot(
            x,
            params_e.e1 * x + params_e.Intercept,
            label=r"y={:.3g}x+{:.3g}, $R^2$={:.3g}".format(
                params_e.e1, params_e.Intercept, res_e.rsquared
            ),
            alpha=0.5,
        )

        ax.set_xlabel(options.msa_label)
        ax.set_ylabel(options.mcmc_label)
        ax.legend(loc="lower right")

        #  ax[0].scatter(x=energies1, y=energies2, color="midnightblue", s=5)
        #  x = np.linspace(*(ax[0]).get_xlim())
        #
        #  ax[0].plot(x, x, "--k", alpha=0.25, zorder=0)
        #  ax[0].plot(
        #      x,
        #      params_e.e1 * x + params_e.Intercept,
        #      label=r"y={:.3g}x+{:.3g}, $R^2$={:.3g}".format(
        #          params_e.e1, params_e.Intercept, res_e.rsquared
        #      ),
        #      alpha=0.5,
        #  )
        #
        #  ax[0].set_xlabel(options.msa_label)
        #  ax[0].set_ylabel(options.mcmc_label)
        #  ax[0].legend(loc="lower right")

        #  ax[1].scatter(
        #      x=energies1_sorted, y=energies2_sorted, color="midnightblue", s=5
        #  )
        #  x = np.linspace(*(ax[1]).get_xlim())
        #
        #  ax[1].plot(x, x, "--k", alpha=0.25, zorder=0)
        #  ax[1].plot(
        #      x,
        #      params_es.e1 * x + params_es.Intercept,
        #      label=r"y={:.3g}x+{:.3g}, $R^2$={:.3g}".format(
        #          params_es.e1, params_es.Intercept, res_es.rsquared
        #      ),
        #      alpha=0.5,
        #  )
        #  ax[1].set_xlabel(options.msa_label + " (sorted)")
        #  ax[1].set_ylabel(options.mcmc_label + " (sorted)")
        #  ax[1].legend(loc="lower right")

        fig.suptitle(options.msa_label + " vs " + options.mcmc_label)

        plt.tight_layout()
        plt.savefig(prefix1 + "_" + prefix2 + "_energies.png")
        plt.close()


if __name__ == "__main__":
    main()
