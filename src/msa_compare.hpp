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

#ifndef MSA_COMPARE_HPP
#define MSA_COMPARE_HPP

#include "msa.hpp"

#include <armadillo>

class MSACompare
{
public:
  MSACompare(MSA*, MSA*, int = BINS);

  void computeFrequency1p(void);
  void computeFrequency2p(void);
  void computeFrequency3p(void);
  void computeFrequency4p(void);
  void computeCorrelation2p(void);
  void computeCorrelation3p(void);
  // void computeCorrelation4p(void);

  void makeFrequency1pHistogram(void);
  void makeFrequency2pHistogram(void);
  void makeFrequency3pHistogram(void);
  void makeCorrelation2pHistogram(void);
  void makeCorrelation3pHistogram(void);
  void makeEfficient3pHistograms(void);
  void makeEfficient4pHistograms(void);

  void writeFrequency1p(void);
  void writeFrequency2p(void);

private:
  int M_msa;              // number of sequences from MSA
  int M_mc;               // number of sequences from MC
  int N;                  // number of positions
  int Q;                  // amino acid alphabet size
  double M_effective_msa; // effect number of sequences from MSA
  double M_effective_mc;  // effect number of sequences from MC

  MSA* msa;
  MSA* mc;

  int bins;

  bool frequency_1p_set = false;
  bool frequency_2p_set = false;
  bool frequency_3p_set = false;
  // bool frequency_4p_set = false;
  bool correlation_2p_set = false;
  bool correlation_3p_set = false;
  // bool correlation_4p_set = false;

  double freq_1p_max = 0;
  double freq_1p_min = 1;
  double freq_2p_max = 0;
  double freq_2p_min = 1;
  double freq_3p_max = 0;
  double freq_3p_min = 1;
  double freq_4p_max = 0;
  double freq_4p_min = 1;
  double corr_2p_max = -1;
  double corr_2p_min = 1;
  double corr_3p_max = -1;
  double corr_3p_min = 1;
  double corr_4p_max = -1;
  double corr_4p_min = 1;

  arma::Col<double> msa_frequency_1p;
  arma::Col<double> mc_frequency_1p;
  arma::Col<double> msa_frequency_2p;
  arma::Col<double> mc_frequency_2p;
  arma::Col<double> msa_frequency_3p;
  arma::Col<double> mc_frequency_3p;
  arma::Col<double> msa_correlation_2p;
  arma::Col<double> mc_correlation_2p;
  arma::Col<double> msa_correlation_3p;
  arma::Col<double> mc_correlation_3p;
};

#endif
