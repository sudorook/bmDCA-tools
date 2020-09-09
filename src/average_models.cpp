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
  std::vector<std::string> params_files;
  std::vector<std::string> params_J_files;
  bool compat_mode = true;

  char c;
  while ((c = getopt(argc, argv, "p:P:")) != -1) {
    switch (c) {
      case 'p':
        params_files.push_back(optarg);
        break;
      case 'P':
        params_J_files.push_back(optarg);
        compat_mode = false;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  if (compat_mode == false) {
    if (params_J_files.size() != params_files.size()) {
      std::cerr << "ERROR: different numbers of fields and couplings provided."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (params_files.size() == 0) {
    std::cerr << "ERROR: no parameters given." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  potts_model params_avg;
  if (compat_mode) {
    params_avg = loadPottsModelAscii(params_files[0]);
  } else {
    params_avg = loadPottsModel(params_files[0], params_J_files[0]);
  }
  int Q = params_avg.h.n_rows;
  int N = params_avg.h.n_cols;

  potts_model params_tmp;
  for (size_t i = 1; i < params_files.size(); i++) {
    if (compat_mode) {
      params_tmp = loadPottsModelAscii(params_files[i]);
    } else {
      params_tmp = loadPottsModel(params_files[i], params_J_files[i]);
    }
    params_avg.h = (double)(i) / (double)(i + 1) * params_avg.h +
                   params_tmp.h / (double)(i + 1);
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        params_avg.J(i, j) =
          (double)(i) / (double)(i + 1) * params_avg.J(i, j) +
          params_tmp.J(i, j) / (double)(i + 1);
      }
    }
  }

  if (compat_mode) {
    std::ofstream output_stream("parameters_avg.bin");
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int aa1 = 0; aa1 < Q; aa1++) {
          for (int aa2 = 0; aa2 < Q; aa2++) {
            output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                          << " " << params_avg.J(i, j)(aa1, aa2) << std::endl;
          }
        }
      }
    }
    for (int i = 0; i < N; i++) {
      for (int aa = 0; aa < Q; aa++) {
        output_stream << "h " << i << " " << aa << " " << params_avg.h(aa, i)
                      << std::endl;
      }
    }
  } else {
    params_avg.h.save("parameters_avg_h.bin", arma::arma_binary);
    params_avg.J.save("parameters_avg_J.bin", arma::arma_binary);
  }
}
