/*
 * SPDX-FileCopyrightText: 2020 sudorook <daemon@nullcodon.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "msa_binary.hpp"

#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef AA_ALPHABET_SIZE
#define AA_ALPHABET_SIZE 21
#endif

MSABin::MSABin(std::string msa_file, bool is_numeric_msa, std::vector<int> keep)
  : keep_seq(keep)
{
  if (is_numeric_msa) {
    readInputNumericMSA(msa_file);
  } else {
    readInputMSA(msa_file);
    M = seq_records.size();
    N = getSequenceLength(seq_records.begin()->getSequence());
    Q = AA_ALPHABET_SIZE;
    makeNumericalMatrix();
  }
};

MSABin::MSABin(arma::Mat<unsigned int> alignment,
               int M,
               int N,
               int Q,
               std::vector<int> keep)
  : alignment(alignment)
  , M(M)
  , N(N)
  , Q(Q)
  , keep_seq(keep){};

void
MSABin::readInputNumericMSA(std::string numeric_msa_file)
{
  std::ifstream input_stream(numeric_msa_file);

  if (!input_stream) {
    std::cerr << "ERROR: couldn't open '" << numeric_msa_file
              << "' for reading." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  input_stream >> M >> N >> Q;
  alignment = arma::Mat<unsigned int>(M, N * Q, arma::fill::zeros);

  int counter = 0;
  int i = 0;
  std::string line;
  std::getline(input_stream, line);
  while (std::getline(input_stream, line)) {
    std::istringstream iss(line);
    int n;
    i = 0;

    while (iss >> n) {
      alignment(counter, (i * Q + n))++;
      i++;
    }
    counter++;
  }
}

void
MSABin::readInputMSA(std::string msa_file)
{
  std::ifstream input_stream(msa_file);

  if (!input_stream) {
    std::cerr << "ERROR: cannot read from '" << msa_file << "'." << std::endl;
    exit(2);
  }

  /*
   * Read a FASTA-formatted multiple sequence alignment. Each record from the
   * file is stored as a SeqRecord object and appended to the seq_records
   * vector.
   */
  std::string header, sequence, line;
  while (input_stream) {
    std::getline(input_stream, line);
    if (line[0] == '>') {
      if (sequence.length() > 0) {
        seq_records.push_back(SeqRecord(header, sequence));
        sequence.clear();
        header.clear();
      }
      header = line;
      header.erase(0, 1);
    } else {
      sequence += line;
      line.clear();
    }
  };
  seq_records.push_back(SeqRecord(header, sequence));
  input_stream.close();
};

void
MSABin::makeNumericalMatrix(void)
{
  alignment = arma::Mat<unsigned int>(M, N * Q, arma::fill::zeros);

  int row_idx = 0;
  for (auto seq = seq_records.begin(); seq != seq_records.end(); seq++) {
    std::string sequence = seq->getSequence();
    int col_idx = 0;
    for (auto aa = sequence.begin(); aa != sequence.end(); aa++) {
      switch (*aa) {
        case '-':
        case 'B':
        case 'J':
        case 'O':
        case 'U':
        case 'X':
        case 'Z':
          alignment(row_idx, (Q * col_idx + 0))++;
          col_idx++;
          break;
        case 'A':
          alignment(row_idx, (Q * col_idx + 1))++;
          col_idx++;
          break;
        case 'C':
          alignment(row_idx, (Q * col_idx + 2))++;
          col_idx++;
          break;
        case 'D':
          alignment(row_idx, (Q * col_idx + 3))++;
          col_idx++;
          break;
        case 'E':
          alignment(row_idx, (Q * col_idx + 4))++;
          col_idx++;
          break;
        case 'F':
          alignment(row_idx, (Q * col_idx + 5))++;
          col_idx++;
          break;
        case 'G':
          alignment(row_idx, (Q * col_idx + 6))++;
          col_idx++;
          break;
        case 'H':
          alignment(row_idx, (Q * col_idx + 7))++;
          col_idx++;
          break;
        case 'I':
          alignment(row_idx, (Q * col_idx + 8))++;
          col_idx++;
          break;
        case 'K':
          alignment(row_idx, (Q * col_idx + 9))++;
          col_idx++;
          break;
        case 'L':
          alignment(row_idx, (Q * col_idx + 10))++;
          col_idx++;
          break;
        case 'M':
          alignment(row_idx, (Q * col_idx + 11))++;
          col_idx++;
          break;
        case 'N':
          alignment(row_idx, (Q * col_idx + 12))++;
          col_idx++;
          break;
        case 'P':
          alignment(row_idx, (Q * col_idx + 13))++;
          col_idx++;
          break;
        case 'Q':
          alignment(row_idx, (Q * col_idx + 14))++;
          col_idx++;
          break;
        case 'R':
          alignment(row_idx, (Q * col_idx + 15))++;
          col_idx++;
          break;
        case 'S':
          alignment(row_idx, (Q * col_idx + 16))++;
          col_idx++;
          break;
        case 'T':
          alignment(row_idx, (Q * col_idx + 17))++;
          col_idx++;
          break;
        case 'V':
          alignment(row_idx, (Q * col_idx + 18))++;
          col_idx++;
          break;
        case 'W':
          alignment(row_idx, (Q * col_idx + 19))++;
          col_idx++;
          break;
        case 'Y':
          alignment(row_idx, (Q * col_idx + 20))++;
          col_idx++;
          break;
      }
    }
    row_idx++;
  }
};

void
MSABin::writeNumericAlignment(std::string output_file)
{
  std::ofstream output_stream(output_file);
  output_stream << M << " " << N << " " << Q << std::endl;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      int aa = 0;
      for (int a = 0; a < Q; a++) {
        aa += alignment(i, (j * Q + a)) * a;
      }
      if (j + 1 == N) {
        output_stream << aa << std::endl;
      } else {
        output_stream << aa << " ";
      }
    }
  }
};

int
MSABin::getSequenceLength(std::string sequence)
{
  int valid_aa_count = 0;
  for (std::string::iterator it = sequence.begin(); it != sequence.end();
       ++it) {
    switch (*it) {
      case '-':
      case 'B':
      case 'J':
      case 'O':
      case 'U':
      case 'X':
      case 'Z':
      case 'A':
      case 'C':
      case 'D':
      case 'E':
      case 'F':
      case 'G':
      case 'H':
      case 'I':
      case 'K':
      case 'L':
      case 'M':
      case 'N':
      case 'P':
      case 'Q':
      case 'R':
      case 'S':
      case 'T':
      case 'V':
      case 'W':
      case 'Y':
        valid_aa_count += 1;
        break;
    }
  }
  return valid_aa_count;
};

// void
// MSABin::filterGaps(double seq_threshold, double pos_threshold)
// {
//   arma::Col<unsigned int> seq_gap_counts =
//     arma::Col<unsigned int>(M, arma::fill::zeros);
//   arma::Col<unsigned int> pos_gap_counts =
//     arma::Col<unsigned int>(N, arma::fill::zeros);
// 
//   unsigned int seq_gap_cutoff = (unsigned int)(seq_threshold * N);
//   unsigned int pos_gap_cutoff = (unsigned int)(pos_threshold * M);
// 
// #pragma omp parallel
//   {
// #pragma omp for
//     for (int i = 0; i < M; i++) {
//       for (int j = 0; j < N; j++ ) {
//         seq_gap_counts(i) += alignment(i, Q * j);
//       }
//     }
//   }
// 
// #pragma omp parallel
//   {
// #pragma omp for
//     for (int j = 0; j < N; j++) {
//       pos_gap_counts(j) = arma::sum(alignment.col(Q * j));
//     }
//   }
// 
//   arma::uvec bad_sequences = arma::find(seq_gap_counts > seq_gap_cutoff);
//   for (size_t i = 0; i < keep_seq.size(); i++) {
//     bad_sequences.shed_rows(arma::find(bad_sequences == keep_seq[i]));
//   }
//   alignment.shed_rows(bad_sequences);
//   
//   arma::uvec bad_positions = arma::find(pos_gap_counts > pos_gap_cutoff);
//   arma::uvec all_bad_positions;
//   arma::uvec tmp = arma::linspace<arma::uvec>(0, Q - 1, Q);
//   for (size_t i = 0; i < bad_positions.n_elem; i++) {
//     all_bad_positions =
//       arma::join_cols(all_bad_positions, tmp + bad_positions(i) * Q);
//   }
//   alignment.shed_cols(all_bad_positions);
// 
//   M = alignment.n_rows;
//   N = (int)alignment.n_cols / Q;
// };

void
MSABin::filterSequenceGaps(double seq_threshold)
{
  arma::Col<unsigned int> seq_gap_counts =
    arma::Col<unsigned int>(M, arma::fill::zeros);

  unsigned int seq_gap_cutoff = (unsigned int)(seq_threshold * N);

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++ ) {
        seq_gap_counts(i) += alignment(i, Q * j);
      }
    }
  }

  arma::uvec bad_sequences = arma::find(seq_gap_counts > seq_gap_cutoff);
  for (size_t i = 0; i < keep_seq.size(); i++) {
    bad_sequences.shed_rows(arma::find(bad_sequences == keep_seq[i]));
  }
  alignment.shed_rows(bad_sequences);
  M = alignment.n_rows;
};

void
MSABin::filterPositionGaps(double pos_threshold)
{
  arma::Col<unsigned int> pos_gap_counts =
    arma::Col<unsigned int>(N, arma::fill::zeros);

  unsigned int pos_gap_cutoff = (unsigned int)(pos_threshold * M);

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      pos_gap_counts(i) = arma::sum(alignment.col(Q * i));
    }
  }
  
  arma::uvec bad_positions = arma::find(pos_gap_counts > pos_gap_cutoff);
  arma::uvec all_bad_positions;
  arma::uvec tmp = arma::linspace<arma::uvec>(0, Q - 1, Q);
  for (size_t i = 0; i < bad_positions.n_elem; i++) {
    all_bad_positions =
      arma::join_cols(all_bad_positions, tmp + bad_positions(i) * Q);
  }
  alignment.shed_cols(all_bad_positions);
  N = (int)alignment.n_cols / Q;
};

void
MSABin::filterSimilarSequences(double threshold)
{
  arma::Col<int> sequence_status = arma::Col<int>(M, arma::fill::zeros);

  unsigned int sim_cutoff = (unsigned int)(N * threshold);
  arma::Mat<unsigned int> alignment_T = alignment.t();
#pragma omp parallel
  {
#pragma omp for schedule(dynamic, 1)
    for (int i = 0; i < M; i++) {
      arma::Col<unsigned int> col = alignment_T.col(i);
      for (int j = 0; j < M; j++) {
        if (j != i) {
          unsigned int sim = arma::dot(col, alignment_T.col(j));
          if (sim > sim_cutoff) {
            sequence_status(i) = 1;
            break;
          }
        }
      }
    }
  }
  arma::uvec bad_sequences = arma::find(sequence_status == 1);
  for (size_t i = 0; i < keep_seq.size(); i++) {
    bad_sequences.shed_rows(arma::find(bad_sequences == keep_seq[i]));
  }
  alignment.shed_rows(bad_sequences);
  M = alignment.n_rows;
};

void
MSABin::computeSequenceWeights(double threshold)
{
  sequence_weights = arma::Col<double>(M, arma::fill::zeros);

  unsigned int sim_cutoff = (unsigned int)(N * threshold);
  arma::Mat<unsigned int> alignment_T = alignment.t();

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < M; i++) {
      arma::Col<unsigned int> col = alignment_T.col(i);
      for (int j = 0; j < M; j++) {
        unsigned int sim = arma::dot(col, alignment_T.col(j));
        if (sim > sim_cutoff) {
          sequence_weights(i)++;
        }
      }
    }
  }
  sequence_weights = 1. / sequence_weights;
};

void
MSABin::writeSequenceWeights(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < M; i++) {
    output_stream << sequence_weights(i) << std::endl;
  }
};
