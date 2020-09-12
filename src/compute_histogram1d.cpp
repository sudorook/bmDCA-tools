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

#define NBINS 51

int
main(int argc, char* argv[])
{
  std::string params_file;
  std::string params_J_file;
  bool compat_mode = true;

  std::string output_h_file;
  std::string output_J_file;

  char c;
  while ((c = getopt(argc, argv, "p:P:")) != -1) {
    switch (c) {
      case 'p':
        params_file = optarg;
        break;
      case 'P':
        params_J_file = optarg;
        compat_mode = false;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  potts_model params;
  if (compat_mode) {
    params = loadPottsModelAscii(params_file);

    // Generate output histogram file names
    int idx = params_file.find_last_of(".");
    int idx2 = params_file.find_last_of("_");
    output_h_file = params_file.substr(0, idx2) + "_h" +
                    params_file.substr(idx2, idx - idx2) + "_hist.tsv";
    output_J_file = params_file.substr(0, idx2) + "_J" +
                    params_file.substr(idx2, idx - idx2) + "_hist.tsv";
  } else {
    params = loadPottsModel(params_file, params_J_file);

    // Generate output histogram file names
    int idx_h = params_file.find_last_of(".");
    int idx_J = params_file.find_last_of(".");
    output_h_file = params_file.substr(0, idx_h) + "_hist.tsv";
    output_J_file = params_J_file.substr(0, idx_J) + "_hist.tsv";
  }

  int Q = params.h.n_rows;
  int N = params.h.n_cols;

  // Check that fields and couplings have matching dimensions
  if ((N != (int)params.J.n_rows) & (N != (int)params.J.n_cols)) {
    std::cerr << "ERROR: parameters dimension mismatch." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((Q != (int)params.J(0, 1).n_cols) & (Q != (int)params.J(0, 1).n_rows)) {
    std::cerr << "ERROR: parameters dimension mismatch." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Compute h histogram
  histogram1d hist_h;
  {
    double h_min = params.h.min();
    double h_max = params.h.max();

    hist_h.range = arma::Col<unsigned long long int>(NBINS, arma::fill::zeros);
    double bin_width = 1.05 * (h_max - h_min) / (double)(NBINS - 1);
    hist_h.bin_width = bin_width;
    hist_h.max = h_max + .025 * (h_max - h_min);
    hist_h.min = h_min - .025 * (h_max - h_min);

    for (int i = 0; i < N; i++) {
      for (int a = 0; a < Q; a++) {
        hist_h.range(floor((params.h(a, i) - hist_h.min) / bin_width + .5))++;
      }
    }
  }

  // Compute J histogram
  histogram1d hist_J;
  {
    double J_min = 0;
    double J_max = 0;
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        if ((i == 0) & (j == 1)) {
          J_min = params.J(i, j).min();
          J_max = params.J(i, j).max();
        } else {
          double tmp_min = params.J(i, j).min();
          double tmp_max = params.J(i, j).max();

          if (tmp_max > J_max) {
            J_max = tmp_max;
          }
          if (tmp_min < J_min) {
            J_min = tmp_min;
          }
        }
      }
    }

    hist_J.range = arma::Col<unsigned long long int>(NBINS, arma::fill::zeros);
    double bin_width = 1.05 * (J_max - J_min) / (double)(NBINS - 1);
    hist_J.bin_width = bin_width;
    hist_J.max = J_max + .025 * (J_max - J_min);
    hist_J.min = J_min - .025 * (J_max - J_min);

    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int a = 0; a < Q; a++) {
          for (int b = 0; b < Q; b++) {
            hist_J.range(
              floor((params.J(i, j)(a, b) - hist_J.min) / bin_width + .5))++;
          }
        }
      }
    }
  }

  // Write output
  writeHistogram1D(output_h_file, hist_h);
  writeHistogram1D(output_J_file, hist_J);
};
