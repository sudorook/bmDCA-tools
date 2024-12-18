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

#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "utils.hpp"

int
main(int argc, char* argv[])
{
  std::string params_file;
  std::string params_J_file;
  bool compat_mode = true;
  bool drop_gaps = false;

  std::string output_h_file;
  std::string output_J_file;

  char c;
  while ((c = getopt(argc, argv, "p:P:g")) != -1) {
    switch (c) {
      case 'p':
        params_file = optarg;
        break;
      case 'P':
        params_J_file = optarg;
        compat_mode = false;
        break;
      case 'g':
        drop_gaps = true;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  std::cout << "reading parameters... " << std::flush;
  potts_model params;
  if (compat_mode) {
    params = loadPottsModelAscii(params_file);
    int idx = params_file.find_last_of(".");
    int idx2 = params_file.find_last_of("_");
    output_h_file = params_file.substr(0, idx2) + "_h" +
                    params_file.substr(idx2, idx - idx2) + "_zg.bin";
    output_J_file = params_file.substr(0, idx2) + "_J" +
                    params_file.substr(idx2, idx - idx2) + "_zg.bin";
  } else {
    params = loadPottsModel(params_file, params_J_file);
    int idx_h = params_file.find_last_of(".");
    int idx_J = params_file.find_last_of(".");
    output_h_file = params_file.substr(0, idx_h) + "_zg.bin";
    output_J_file = params_J_file.substr(0, idx_J) + "_zg.bin";
  }
  std::cout << "done" << std::endl;

  int Q = params.h.n_rows;
  int N = params.h.n_cols;

  std::cout << "rescaling gauge... " << std::flush;
  potts_model params_zg;

  // initialize
  params_zg.h = arma::Mat<double>(Q, N, arma::fill::zeros);
  params_zg.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params_zg.J(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

  // rescale couplings
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i < j) {
        double J_ij_mean = arma::mean(arma::mean(params.J(i, j)));
        arma::Mat<double> J_ija_mean = arma::mean(params.J(i, j), 1);
        arma::Mat<double> J_ijb_mean = arma::mean(params.J(i, j), 0);
        for (int a = 0; a < Q; a++) {
          for (int b = 0; b < Q; b++) {
            params_zg.J(i, j)(a, b) =
              params.J(i, j)(a, b) - J_ija_mean(a) - J_ijb_mean(b) + J_ij_mean;
          }
          params_zg.h(a, i) += J_ija_mean(a) - J_ij_mean;
        }
      } else if (i > j) {
        double J_ij_mean = arma::mean(arma::mean(params.J(j, i)));
        arma::Mat<double> J_ija_mean = arma::mean(params.J(j, i), 0);
        arma::Mat<double> J_ijb_mean = arma::mean(params.J(j, i), 1);
        for (int a = 0; a < Q; a++) {
          params_zg.h(a, i) += J_ija_mean(a) - J_ij_mean;
        }
      }
    }
  }

  // rescale fields
  arma::Row<double> h_i_mean = arma::mean(params.h, 0);
  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      params_zg.h(a, i) += params.h(a, i) - h_i_mean(i);
    }
  }
  std::cout << "done" << std::endl;

  // saving parameters
  if (drop_gaps == false) {
    std::cout << "writing output... " << std::flush;
    params_zg.h.save(output_h_file, arma::arma_binary);
    params_zg.J.save(output_J_file, arma::arma_binary);
    std::cout << "done" << std::endl;
  } else {
    std::cout << "dropping gaps... " << std::flush;
    potts_model params_zg_nogap;
    params_zg_nogap.h = arma::Mat<double>(Q - 1, N, arma::fill::zeros);
    params_zg_nogap.J = arma::field<arma::Mat<double>>(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        params_zg_nogap.J(i, j) =
          arma::Mat<double>(Q - 1, Q - 1, arma::fill::zeros);
      }
    }
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int a = 1; a < Q; a++) {
          for (int b = 1; b < Q; b++) {
            params_zg_nogap.J(i, j)(a - 1, b - 1) = params_zg.J(i, j)(a, b);
          }
        }
      }
    }
    for (int i = 0; i < N; i++) {
      for (int a = 1; a < Q; a++) {
        params_zg_nogap.h(a - 1, i) = params_zg.h(a, i);
      }
    }
    std::cout << "done" << std::endl;

    std::cout << "writing output... " << std::flush;
    params_zg_nogap.h.save(output_h_file, arma::arma_binary);
    params_zg_nogap.J.save(output_J_file, arma::arma_binary);
    std::cout << "done" << std::endl;
  }
};
