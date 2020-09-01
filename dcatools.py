#! /usr/bin/env python3
""" docstring """

import re
import subprocess
from collections import defaultdict
import numpy as np
import pandas as pd
import unidecode
from scipy.io import loadmat

from Bio.PDB import PDBParser
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio import SeqIO
from Bio import AlignIO


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
    if norm == 2:  # frob norm
        for i, pos in enumerate(pos_range):
            for aa in range(Naa):
                h_new[i] += (h[pos, aa]) ** 2
            h_new[i] = h_new[i] ** 0.5
    if norm == 0:  # mean
        h_new = np.mean(h[pos_range, :], 1)
    return h_new


def compute_energies(seqs, h, J):
    """ compute statistical energies for sequences """
    energies = np.zeros(seqs.shape[0])
    N = len(seqs[0])
    for i, seq in enumerate(seqs):
        for n1 in range(0, N):
            energies[i] += -h[n1][seq[n1]]
            for n2 in range(n1 + 1, N):
                energies[i] += -J[n1][n2][seq[n1]][seq[n2]]
    return energies


def compute_pdb_distance_matrix(structure, chain_id="A", mode="CA"):
    """ compute PDB distance matrix """

    chain = structure[0][chain_id]
    labels = [
        ("".join(map(str, i.get_id()[1:3]))).strip()
        for i in chain.get_residues()
        if i.get_id()[0] == " "
    ]
    residues = [aa for aa in chain.get_residues() if aa.get_id()[0] == " "]

    Naa = len(residues)
    distmat = np.zeros((Naa, Naa))

    if mode == "CA":
        atoms = [
            atom
            for residue in residues
            for atom in residue.get_atoms()
            if atom.get_id() == "CA"
        ]
        for i in range(0, Naa):
            for j in range(i + 1, Naa):
                distmat[i, j] = atoms[i] - atoms[j]

    elif mode == "min":
        for i in range(0, Naa):
            for j in range(i + 1, Naa):
                distmat[i, j] = min(
                    [
                        atom0 - atom1
                        for atom0 in residues[i].get_atoms()
                        for atom1 in residues[j].get_atoms()
                    ]
                )

    return distmat, labels


def get_pdb_sequence(structure, chain_id="A"):
    """ compute PDB distance matrix """

    aa_table = defaultdict(lambda: "X")
    aa_table["ALA"] = "A"
    aa_table["ARG"] = "R"
    aa_table["ASN"] = "N"
    aa_table["ASP"] = "D"
    aa_table["CYS"] = "C"
    aa_table["GLN"] = "Q"
    aa_table["GLU"] = "E"
    aa_table["GLY"] = "G"
    aa_table["HIS"] = "H"
    aa_table["ILE"] = "I"
    aa_table["LEU"] = "L"
    aa_table["LYS"] = "K"
    aa_table["MET"] = "M"
    aa_table["PHE"] = "F"
    aa_table["PRO"] = "P"
    aa_table["SER"] = "S"
    aa_table["THR"] = "T"
    aa_table["TRP"] = "W"
    aa_table["TYR"] = "Y"
    aa_table["VAL"] = "V"

    chain = structure[0][chain_id]
    sequence = [
        aa_table[aa.resname]
        for aa in chain.get_residues()
        if aa.get_id()[0] == " "
    ]

    return Seq("".join(sequence), generic_protein)


def compute_contact_matrix(distmat, threshold=8.0):
    """ convert distance matrix to contact matrix """
    return distmat * (distmat != 0) * (distmat < threshold)


def load_sequences(msa_file):
    """ load sequences """
    data = np.loadtxt(msa_file, dtype="int", skiprows=1)
    return data


def load_sequences_fasta(sequence_file, sequence_format="fasta"):
    """ load FASTA-formatted sequences """
    sequences = list(SeqIO.parse(sequence_file, sequence_format))
    return sequences


def load_alignment_fasta(alignment_file, alignment_format="fasta"):
    """ Load the alignment and return it. """
    alignment = AlignIO.read(alignment_file, alignment_format)
    return alignment


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


def load_model_mat(data_file):
    """ load mat file with alignment, fields, and couplings """
    data = loadmat(data_file)
    alignment = data["align"]
    params_h = data["h"]
    params_J = data["J"]
    return alignment, params_h, params_J


def load_pdb_structure(pdb_id, pdb_file):
    """ load pdb structure """
    parser = PDBParser(PERMISSIVE=1)
    return parser.get_structure(pdb_id, pdb_file)


def load_energies(energy_file):
    """ load sequences energies """
    data = np.loadtxt(energy_file, dtype="double", skiprows=1)
    return data


def load_distances(distance_file):
    """ load sequence distances """
    data = np.loadtxt(distance_file, dtype="double", skiprows=1)
    return data


def save_alignment(alignment, M, N, Q, filename):
    """ save numerical alignment """
    with open(filename, "w") as handle:
        handle.write("%d %d %d\n" % (M, N, Q))
        for i in range(0, M):
            line = " ".join([str(x) for x in alignment[i, :]])
            handle.write(line)
            handle.write("\n")


def save_sequences_fasta(sequence, sequence_file, sequence_format="fasta"):
    """ Save sequences to disk. """
    with open(sequence_file, "w") as handle:
        SeqIO.write(
            SeqRecord(sequence, id="PDB sequence"), handle, sequence_format
        )


def save_parameters(h, J, M, N, Q, filename):
    """ save parameters (h, J) to file """
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


def find_reference_sequence_ggsearch(
    sequence_reference_fasta, msa_fasta, N_pos
):
    """ use ggsearch to find the reference sequence in an MSA """
    command = [
        "ggsearch36",
        "-M 1-" + str(N_pos),
        "-b",
        "1",
        "-m 8",
        sequence_reference_fasta,
        msa_fasta,
    ]
    res = subprocess.check_output(command)
    align_id = res.decode("ASCII").split("\t")[1]
    alignment = load_alignment_fasta(msa_fasta)
    sequence = [seq.seq for seq in alignment if seq.id == align_id][0]
    return sequence


def find_positions_to_keep(target_sequence, reference_sequence):
    """ align 2 sequences and return list of ungapped positons to keep """
    target_strip = str(target_sequence).replace("-", "")
    target_start_positions = [
        i for i in range(0, len(target_sequence)) if target_sequence[i] != "-"
    ]
    reference_strip = str(reference_sequence).replace("-", "")
    reference_start_positions = [
        i
        for i in range(0, len(reference_sequence))
        if reference_sequence[i] != "-"
    ]

    target_align, reference_align, _, _, _ = pairwise2.align.globalms(
        target_strip, reference_strip, 2, -1, -0.5, -0.1
    )[0]

    i_target = 0
    i_reference = 0
    target_positions = []
    reference_positions = []
    for i, aa in enumerate(target_align):
        if aa != "-" and reference_align[i] != "-":
            target_positions.append(i_target)
            reference_positions.append(i_reference)
        if aa != "-":
            i_target += 1
        if reference_align[i] != "-":
            i_reference += 1

    return_target_positions = [
        target_start_positions[i] for i in target_positions
    ]
    return_reference_positions = [
        reference_start_positions[i] for i in reference_positions
    ]

    return return_target_positions, return_reference_positions


def subset_model(h, J, position_list):
    """ subset the positions in a Potts model """
    N_pos = len(position_list)
    N_aa = h.shape[1]

    h_subset = h[position_list, :]

    J_subset = np.zeros((N_pos, N_pos, N_aa, N_aa))
    for i, pos1 in enumerate(position_list):
        for j, pos2 in enumerate(position_list):
            J_subset[i, j, :, :] = J[pos1, pos2, :, :]

    return h_subset, J_subset


def subset_distance_matrix(distmat, position_list):
    """ subset the distance matrix """
    N_pos = len(position_list)
    distmat_subset = np.zeros((N_pos, N_pos))
    for i, pos1 in enumerate(position_list):
        for j, pos2 in enumerate(position_list):
            distmat_subset[i, j] = distmat[pos1, pos2]
    return distmat_subset


def slugify(text):
    """ convert text string to slug """
    text = unidecode.unidecode(text).lower()
    return re.sub(r"[\W_]+", "_", text)


def arma2ascii(h_file, J_file):
    """ convert armadillo binary to text file """

    command = ["arma2ascii", "-p", h_file, "-P", J_file]
    res = subprocess.call(command)
    return res


def numeric2fasta(numeric_msa_file, output_fasta_file):
    """ convert numeric MSA to a FASTA file """
    command = [
        "numeric2fasta",
        "-n",
        numeric_msa_file,
        "-o",
        output_fasta_file,
    ]
    res = subprocess.call(command)
    return res


def load_run_log(log_file, offset=None):
    """ load the run log """

    if offset is not None:
        offset = int(offset)
        df = pd.read_csv(log_file, sep="\t")[offset:]
    else:
        df = pd.read_csv(log_file, sep="\t")

    df["log10-burn-in"] = np.log10(df["burn-in"])
    df["log10-burn-between"] = np.log10(df["burn-between"])

    return df
