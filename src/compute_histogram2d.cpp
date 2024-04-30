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
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.hpp"
#include "utils.hpp"

#define NBINS 101

std::string
generateOutputString(std::string name1, std::string name2)
{

  std::vector<std::string> dirlist1;
  {
    std::stringstream test(name1);
    std::string segment;
    while (std::getline(test, segment, '/')) {
      dirlist1.push_back(segment);
    }
  }

  std::vector<std::string> seglist1;
  {
    std::stringstream test(dirlist1[dirlist1.size() - 1]);
    std::string segment;
    while (std::getline(test, segment, '_')) {
      seglist1.push_back(segment);
    }
  }

  std::vector<std::string> dirlist2;
  {
    std::stringstream test(name2);
    std::string segment;
    while (std::getline(test, segment, '/')) {
      dirlist2.push_back(segment);
    }
  }

  std::vector<std::string> seglist2;
  {
    std::stringstream test(dirlist2[dirlist2.size() - 1]);
    std::string segment;
    while (std::getline(test, segment, '_')) {
      seglist2.push_back(segment);
    }
  }

  std::string output_file;
  for (int i = 0; i < Min(seglist2.size(), seglist1.size()); i++) {
    if (seglist1[i] == seglist2[i]) {
      output_file += seglist1[i] + "_";
    }
  }
  if (output_file.back() == '_')
    output_file.pop_back();

  return output_file;
};

int
main(int argc, char* argv[])
{
  std::string params_file;
  std::string params_J_file;
  bool compat_mode = true;

  std::string params2_file;
  std::string params2_J_file;
  bool compat2_mode = true;

  std::string hist_h_file;
  std::string hist_J_file;
  std::string model_h_file;
  std::string model_J_file;

  char c;
  while ((c = getopt(argc, argv, "p:P:q:Q:")) != -1) {
    switch (c) {
      case 'p':
        params_file = optarg;
        break;
      case 'P':
        params_J_file = optarg;
        compat_mode = false;
        break;
      case 'q':
        params2_file = optarg;
        break;
      case 'Q':
        params2_J_file = optarg;
        compat2_mode = false;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  potts_model params1;
  if (compat_mode) {
    params1 = loadPottsModelAscii(params_file);

    // Generate output histogram file names
    int idx = params_file.find_last_of(".");
    int idx2 = params_file.find_last_of("_");
    hist_h_file = params_file.substr(0, idx2) + "_h" +
                  params_file.substr(idx2, idx - idx2) + "_hist.tsv";
    hist_J_file = params_file.substr(0, idx2) + "_J" +
                  params_file.substr(idx2, idx - idx2) + "_hist.tsv";
    model_h_file = params_file.substr(0, idx2) + "_h" +
                   params_file.substr(idx2, idx - idx2) + "_model.tsv";
    model_J_file = params_file.substr(0, idx2) + "_J" +
                   params_file.substr(idx2, idx - idx2) + "_model.tsv";
  } else {
    params1 = loadPottsModel(params_file, params_J_file);

    // Generate output histogram file names
    int idx_h = params_file.find_last_of(".");
    int idx_J = params_file.find_last_of(".");
    hist_h_file = params_file.substr(0, idx_h) + "_hist.tsv";
    hist_J_file = params_J_file.substr(0, idx_J) + "_hist.tsv";
    model_h_file = params_file.substr(0, idx_h) + "_model.tsv";
    model_J_file = params_J_file.substr(0, idx_J) + "_model.tsv";
  }

  potts_model params2;
  if (compat2_mode) {
    params2 = loadPottsModelAscii(params2_file);

    // Generate output histogram file names
    int idx = params_file.find_last_of(".");
    int idx2 = params_file.find_last_of("_");

    std::string hist2_h_file = params2_file.substr(0, idx2) + "_h" +
                               params2_file.substr(idx2, idx - idx2) +
                               "_hist.tsv";
    std::string hist2_J_file = params2_file.substr(0, idx2) + "_J" +
                               params2_file.substr(idx2, idx - idx2) +
                               "_hist.tsv";
    hist_h_file = generateOutputString(hist_h_file, hist2_h_file);
    hist_J_file = generateOutputString(hist_J_file, hist2_J_file);

    std::string model2_h_file = params2_file.substr(0, idx2) + "_h" +
                                params2_file.substr(idx2, idx - idx2) +
                                "_model.tsv";
    std::string model2_J_file = params2_file.substr(0, idx2) + "_J" +
                                params2_file.substr(idx2, idx - idx2) +
                                "_model.tsv";
    model_h_file = generateOutputString(model_h_file, model2_h_file);
    model_J_file = generateOutputString(model_J_file, model2_J_file);
  } else {
    params2 = loadPottsModel(params2_file, params2_J_file);

    // Generate output histogram file names
    int idx_h = params2_file.find_last_of(".");
    int idx_J = params2_file.find_last_of(".");

    std::string hist2_h_file = params2_file.substr(0, idx_h) + "_hist.tsv";
    std::string hist2_J_file = params2_J_file.substr(0, idx_J) + "_hist.tsv";
    hist_h_file = generateOutputString(hist_h_file, hist2_h_file);
    hist_J_file = generateOutputString(hist_J_file, hist2_J_file);

    std::string model2_h_file = params2_file.substr(0, idx_h) + "_model.tsv";
    std::string model2_J_file = params2_J_file.substr(0, idx_J) + "_model.tsv";
    model_h_file = generateOutputString(model_h_file, model2_h_file);
    model_J_file = generateOutputString(model_J_file, model2_J_file);
  }

  int Q = params1.h.n_rows;
  int N = params1.h.n_cols;

  // Check that fields and couplings have matching dimensions
  {
    if ((N != (int)params1.J.n_rows) | (N != (int)params1.J.n_cols)) {
      std::cerr << "ERROR: parameters dimension mismatch." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if ((Q != (int)params1.J(0, 1).n_cols) |
        (Q != (int)params1.J(0, 1).n_rows)) {
      std::cerr << "ERROR: parameters dimension mismatch." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    int Q2 = params2.h.n_rows;
    int N2 = params2.h.n_cols;

    if ((N != N2) | (Q != Q2)) {
      std::cerr << "ERROR: two parameter sets have mismatched dimensions."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    if ((N2 != (int)params2.J.n_rows) | (N2 != (int)params2.J.n_cols)) {
      std::cerr << "ERROR: parameters dimension mismatch." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if ((Q2 != (int)params2.J(0, 1).n_cols) |
        (Q2 != (int)params2.J(0, 1).n_rows)) {
      std::cerr << "ERROR: parameters dimension mismatch." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // Compute h histogram
  histogram2d hist_h;
  linear_model model_h;
  {
    double h_min = Min(params1.h.min(), params2.h.min());
    double h_max = Max(params1.h.max(), params2.h.max());

    hist_h.grid =
      arma::Mat<unsigned long long int>(NBINS, NBINS, arma::fill::zeros);
    double bin_width = 1.05 * (h_max - h_min) / (double)(NBINS - 1);
    hist_h.bin_width = bin_width;
    hist_h.max = h_max + .025 * (h_max - h_min);
    hist_h.min = h_min - .025 * (h_max - h_min);

    double x_mean = arma::mean(arma::mean(params1.h));
    double y_mean = arma::mean(arma::mean(params2.h));

    double xx = arma::accu(arma::pow((params1.h - x_mean), 2));
    double xy = arma::accu((params1.h - x_mean) % (params2.h - y_mean));
    double yy = arma::accu(arma::pow((params2.h - y_mean), 2));
    double b = xy / xx;
    double a = y_mean - b * x_mean;
    double r = b * sqrt(xx / yy);

    model_h.a = a;
    model_h.b = b;
    model_h.R2 = pow(r, 2);

    std::cout << "model h: y=" << b << "x+" << a << " (r=" << r << ")"
              << std::endl;

    for (int i = 0; i < N; i++) {
      for (int a = 0; a < Q; a++) {
        hist_h.grid(floor((params1.h(a, i) - hist_h.min) / bin_width + .5),
                    floor((params2.h(a, i) - hist_h.min) / bin_width + .5))++;
      }
    }
  }

  // Compute J histogram
  histogram2d hist_J;
  linear_model model_J;
  {
    double J_min = 0;
    double J_max = 0;
    arma::Mat<double> x_mean_mat = arma::Mat<double>(N, N, arma::fill::zeros);
    arma::Mat<double> y_mean_mat = arma::Mat<double>(N, N, arma::fill::zeros);
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        x_mean_mat(i, j) = arma::mean(arma::mean(params1.J(i, j)));
        y_mean_mat(i, j) = arma::mean(arma::mean(params1.J(i, j)));
        if ((i == 0) & (j == 1)) {
          J_min = Min(params1.J(i, j).min(), params2.J(i, j).min());
          J_max = Max(params1.J(i, j).max(), params2.J(i, j).max());
        } else {
          double tmp_min = Min(params1.J(i, j).min(), params2.J(i, j).min());
          double tmp_max = Max(params1.J(i, j).max(), params2.J(i, j).max());

          if (tmp_max > J_max) {
            J_max = tmp_max;
          }
          if (tmp_min < J_min) {
            J_min = tmp_min;
          }
        }
      }
    }
    double x_mean = arma::accu(x_mean_mat) / (double)(N * (N - 1.) / 2.);
    double y_mean = arma::accu(y_mean_mat) / (double)(N * (N - 1.) / 2.);
    ;
    arma::Mat<double> xx_mat = arma::Mat<double>(N, N, arma::fill::zeros);
    arma::Mat<double> xy_mat = arma::Mat<double>(N, N, arma::fill::zeros);
    arma::Mat<double> yy_mat = arma::Mat<double>(N, N, arma::fill::zeros);
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        xx_mat(i, j) = arma::accu(arma::pow((params1.J(i, j) - x_mean), 2));
        xy_mat(i, j) =
          arma::accu((params1.J(i, j) - x_mean) % (params2.J(i, j) - y_mean));
        yy_mat(i, j) = arma::accu(arma::pow((params2.J(i, j) - y_mean), 2));
      }
    }
    double xx = arma::accu(xx_mat);
    double xy = arma::accu(xy_mat);
    double yy = arma::accu(yy_mat);

    double b = xy / xx;
    double a = y_mean - b * x_mean;
    double r = b * sqrt(xx / yy);

    model_J.a = a;
    model_J.b = b;
    model_J.R2 = pow(r, 2);

    std::cout << "model J: y=" << b << "x+" << a << " (r=" << r << ")"
              << std::endl;

    hist_J.grid =
      arma::Mat<unsigned long long int>(NBINS, NBINS, arma::fill::zeros);
    double bin_width = 1.05 * (J_max - J_min) / (double)(NBINS - 1);
    hist_J.bin_width = bin_width;
    hist_J.max = J_max + .025 * (J_max - J_min);
    hist_J.min = J_min - .025 * (J_max - J_min);

    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int a = 0; a < Q; a++) {
          for (int b = 0; b < Q; b++) {
            hist_J.grid(
              floor((params1.J(i, j)(a, b) - hist_J.min) / bin_width + .5),
              floor((params2.J(i, j)(a, b) - hist_J.min) / bin_width + .5))++;
          }
        }
      }
    }
  }

  // Write output
  writeHistogram2D(hist_h_file, hist_h);
  writeLinearModel(model_h_file, model_h);
  writeHistogram2D(hist_J_file, hist_J);
  writeLinearModel(model_J_file, model_J);
};
