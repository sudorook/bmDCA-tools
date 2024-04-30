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

#ifndef MSA_HPP
#define MSA_HPP

#include <armadillo>
#include <string>
#include <vector>

#include "utils.hpp"

class MSA
{
public:
  arma::Mat<int> alignment;           // numerical multiple sequence alignment
  arma::Col<double> sequence_weights; // weights for each sequence
  int M;                              // number of sequences
  int N;                              // number of positions
  int Q;                              // number of amino acids
  const bool reweight;                // whether reweighting performed
  const double threshold;             // reweighting threshold

  arma::Col<double> max_similarity;
  arma::Col<double> mean_similarity;

  MSA(std::string, std::string = "", bool = true, bool = false, double = 0.8);
  MSA(arma::Mat<int>, int, int, int, bool = true, double = 0.8);
  void printAlignment();
  void writeMatrix(std::string);
  void writeSequenceWeights(std::string);
  void writeSequenceSimilarity(std::string, std::string);

  void computeSequenceSimilarity(void);

private:
  std::vector<SeqRecord> seq_records;
  int getSequenceLength(std::string);
  void readInputMSA(std::string);
  void readInputNumericMSA(std::string);
  void readSequenceWeights(std::string);
  void makeNumericalMatrix(void);
  void computeSequenceWeights(double);
};

#endif
