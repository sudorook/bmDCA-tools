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

  std::cout << "reading parameters... " << std::flush;
  potts_model params;
  if (compat_mode) {
    params = loadPottsModelAscii(params_file);
  } else {
    params = loadPottsModel(params_file, params_J_file);
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
    for (int j = i + 1; j < N; j++) {
      double J_ij_mean = arma::mean(arma::mean(params.J(i, j)));
      arma::Mat<double> J_ija_mean = arma::mean(params.J(i, j), 1);
      arma::Mat<double> J_ijb_mean = arma::mean(params.J(i, j), 0);
      for (int a = 0; a < Q; a++) {
        for (int b = 0; b < Q; b++) {
          params_zg.J(i, j)(a, b) =
            params.J(i, j)(a, b) - J_ija_mean(a) - J_ijb_mean(b) + J_ij_mean;
        }
      }
    }
  }

  // rescale fields
  params_zg.h = params.h;
  arma::Mat<double> h_i_mean = arma::mean(params.h, 0);
  for (int i = 0; i < N; i++) {
    arma::Col<double> J_ia = arma::Col<double>(N, arma::fill::zeros);
    for (int j = 0; j < N; j++) {
      if (j > i) {
        J_ia(j) = arma::mean(arma::mean(params.J(i, j)));
      } else if (i > j) {
        J_ia(j) = arma::mean(arma::mean(params.J(j, i)));
      }
    }
    arma::Col<double> J_i = arma::Col<double>(Q, arma::fill::zeros);
    for (int b = 0; b < Q; b++) {
      arma::Col<double> J_ib = arma::Col<double>(N, arma::fill::zeros);
      for (int j = 0; j < N; j++) {
        arma::Col<double> J_ijb = arma::Col<double>(Q, arma::fill::zeros);
        for (int a = 0; a < Q; a++) {
          if (j > i) {
            J_ijb(a) = params.J(i, j)(a, b);
          } else if (i > j) {
            J_ijb(a) = params.J(j, i)(a, b);
          }
        }
        J_ib(j) = arma::mean(J_ijb);
      }
      J_i(b) = arma::sum(J_ib);
    }
    for (int a = 0; a < Q; a++) {
      params_zg.h(a, i) += J_i(a) - arma::sum(J_ia) - h_i_mean(i);
    }
  }
  std::cout << "done" << std::endl;

  // checking means
  double max = 0;
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int a = 0; a < Q; a++) {
        double sum = 0;
        for (int b = 0; b < Q; b++) {
          sum += params_zg.J(i, j)(a, b);
        }
        std::cout << "J(" << i << "," << j << ")(" << a << ",:) mean = " << sum
                  << std::endl;
        if (sum > max) max = sum;
      }
      for (int b = 0; b < Q; b++) {
        double sum = 0;
        for (int a = 0; a < Q; a++) {
          sum += params_zg.J(i, j)(a, b);
        }
        std::cout << "J(" << i << "," << j << ")(:," << b << ") mean = " << sum
                  << std::endl;
        if (sum > max) max = sum;
      }
    }
  }
  std::cout << "Jmax = " << max << std::flush;

  arma::sum(params_zg.h, 1).print("h_zg_means:");
  arma::sum(params_zg.h, 0).print("h_zg_means:");

  // saving parameters
  std::cout << "writing output... " << std::flush;
  params_zg.h.save("zero_gauge_h.bin", arma::arma_binary);
  params_zg.J.save("zero_gauge_J.bin", arma::arma_binary);
  std::cout << "done" << std::endl;
}
