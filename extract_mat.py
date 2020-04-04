#! /usr/bin/env python3
""" docstring """


import argparse
import os
from itertools import combinations
import numpy as np
from scipy.io import loadmat


def parse_options():
    """ cli parser """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", dest="input", required=True, help="mat file"
    )
    #  parser.add_argument(
    #      "-o", "--output", dest="output", required=True, help="output file name"
    #  )
    return parser.parse_args()


def load_model(data_file):
    """ load mat file with alignment, fields, and couplings """
    data = loadmat(data_file)
    alignment = data["align"]
    #  params_h = data['h'].T
    #  params_J = data['J'].T
    params_h = data["h"]
    params_J = data["J"]

    return alignment, params_h, params_J


def save_alignment(alignment, filename):
    """ save numerical alignment """
    np.savetxt(filename, alignment, fmt="%d")


def save_parameters(h, J, filename):
    """ save parameters (h, J) to file """
    Q, N = h.shape
    with open(filename, "w") as handle:
        for i in range(N):
            for j in range(i + 1, N):
                for a in range(Q):
                    for b in range(Q):
                        handle.write(
                            "J %d %d %d %d %lf\n" % (i, j, a, b, J[a, b, i, j])
                        )
        for i in range(N):
            for a in range(Q):
                handle.write("h %d %d %lf\n" % (i, a, h[a, i]))


def main():
    """ do stuff """
    options = parse_options()
    prefix = os.path.splitext((options.input).replace("/", "_"))[0]

    alignment, h, J = load_model(options.input)
    #  print(alignment)
    #  print(h)
    #  print(J)

    alignment_file = prefix + "_msa_numerical.txt"
    save_alignment(alignment, alignment_file)

    parameters_file = prefix + "_parameters.txt"
    save_parameters(h, J, parameters_file)


if __name__ == "__main__":
    main()
