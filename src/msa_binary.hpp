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

#ifndef MSA_BINARY_HPP
#define MSA_BINARY_HPP

#include <armadillo>
#include <string>
#include <vector>

#include "utils.hpp"

class MSABin
{
public:
  arma::Mat<unsigned int> alignment;
  arma::Col<double> sequence_weights;
  int M;                              // number of sequences
  int N;                              // number of positions
  int Q;                              // number of amino acids

  MSABin(std::string, bool, std::vector<int>);
  MSABin(arma::Mat<unsigned int>, int, int, int, std::vector<int>);

  void writeNumericAlignment(std::string);
  void writeSequenceWeights(std::string);

  // void filterGaps(double = 0.2, double = 0.2);
  void filterSequenceGaps(double = 0.2);
  void filterPositionGaps(double = 0.2);
  void filterSimilarSequences(double = 0.8);
  void computeSequenceWeights(double = 0.8);

private:
  std::vector<SeqRecord> seq_records;
  int getSequenceLength(std::string);
  void readInputMSA(std::string);
  void readInputNumericMSA(std::string);
  void makeNumericalMatrix(void);

  std::vector<int> keep_seq;
};

#endif
