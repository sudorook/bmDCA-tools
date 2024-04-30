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

#ifndef MSA_STATS_HPP
#define MSA_STATS_HPP

#include "msa.hpp"

#include <armadillo>

class MSAStats
{
public:
  MSAStats(MSA*);

  double getEffectiveM();
  double getN();
  double getM();
  double getQ();

  void computeErrorMSA(int = 100, long int = 0);

  void writeRelEntropy(std::string);
  void writeRelEntropyAscii(std::string);
  // void writeRelEntropyPos(std::string);
  void writeRelEntropyPosAscii(std::string);
  void writeRelEntropyGradient(std::string);
  void writeRelEntropyGradientAscii(std::string);
  void writeFrequency1p(std::string);
  void writeFrequency2p(std::string);
  void writeCorrelation2p(std::string);
  void writeFrequency1pAscii(std::string);
  void writeFrequency2pAscii(std::string);
  void writeCorrelation2pAscii(std::string);

  arma::Mat<double> frequency_1p;
  arma::field<arma::Mat<double>> frequency_2p;
  arma::field<arma::Mat<double>> correlation_2p;
  arma::Mat<double> rel_entropy_1p;
  arma::Col<double> rel_entropy_pos_1p;
  arma::Mat<double> rel_entropy_grad_1p;

  double freq_rms;
  arma::Col<double> msa_rms;

private:
  double pseudocount;
  int M;              // number of sequences
  int N;              // number of positions
  int Q;              // amino acid alphabet size
  double M_effective; // effect number of sequences

  MSA* msa;

  arma::Col<double> aa_background_frequencies;
};

#endif
