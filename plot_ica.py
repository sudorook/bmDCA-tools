#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2020 - 2023 sudorook <daemon@nullcodon.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

""" docstring """

import argparse
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from pysca import scaTools as sca

np.random.seed(0)

sns.set(rc={"figure.figsize": (10, 7.5)})


def parse_options():
    """cli parser"""
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
    """reduce rank 4 tensor using 2-norm"""
    Npos = np.shape(J)[0]
    J_new = np.zeros((Npos, Npos))
    for i in range(Npos):
        for j in range(i + 1, Npos):
            J_new[i, j] = np.linalg.norm(J[i, j], norm)
    return J_new


def load_data(data_file):
    """syntax checker is annoying me"""

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
    """do stuff"""
    options = parse_options()

    prefix = os.path.splitext((options.params).replace("/", "_"))[0]

    h, J = load_data(options.params)
    #  J_2d = reduce_tensor(J)
    J_2d = reduce_tensor(J, "fro")
    J_2d = J_2d + np.transpose(J_2d)
    print(J_2d)

    u, s, v = sca.svdss(J_2d, k=6)
    #  u, s, v = np.linalg.svd(J_2d)
    print(u)
    print(s)
    print(v)

    [u_ica, w] = sca.rotICA(u, kmax=6)
    print(u_ica.shape)
    print(w.shape)


if __name__ == "__main__":
    main()
