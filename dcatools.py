#! /usr/bin/env python3
""" docstring """

import numpy as np


def reduce_J(J, norm=2):
    """ reduce J using Frobenius-norm """
    Npos = np.shape(J)[0]
    J_new = np.zeros((Npos, Npos))
    for i in range(Npos):
        for j in range(i + 1, Npos):
            J_new[i, j] = np.linalg.norm(J[i, j], norm)
    return J_new


def reduce_h(h, norm=2, gap_pos=-1):
    """ reduce h using Frobenius-norm """
    Npos, Naa = np.shape(h)
    pos_range = list(x for x in range(Npos) if x != gap_pos)
    h_new = np.zeros(len(pos_range))
    if norm == 2: # frob norm
        for i, pos in enumerate(pos_range):
            for aa in range(Naa):
                h_new[i] += (h[pos, aa]) ** 2
            h_new[i] = h_new[i] ** 0.5
    if norm == 0: # mean
        h_new = np.mean(h[pos_range, :], 1)
    return h_new


def compute_energies(seqs, h, J):
    energies = np.zeros(seqs.shape[0])
    N = len(seqs[0])
    for i, seq in enumerate(seqs):
        for n1 in range(0, N):
            energies[i] += -h[n1][seq[n1]]
            for n2 in range(n1 + 1, N):
                energies[i] += -J[n1][n2][seq[n1]][seq[n2]]
    return energies


def load_sequences(msa_file):
    """ load sequences """
    data = np.loadtxt(msa_file, dtype="int", skiprows=1)
    return data


def load_model(data_file):
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


def load_energies(energy_file):
    """ load sequences energies """
    data = np.loadtxt(energy_file, dtype="double", skiprows=1)
    return data
