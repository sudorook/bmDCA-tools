#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.hpp"
#include "utils.hpp"

int main(int argc, char* argv[]) {
  std::string msa_file;
  std::string energy_file = "sequence_energies.txt";
  std::string params_file;
  std::string params_J_file;
  bool compat_mode = true;

  char c;
  while ((c = getopt(argc, argv, "i:o:p:P:")) != -1) {
    switch (c) {
      case 'i':
        msa_file = optarg;
        break;
      case 'o':
        energy_file = optarg;
        break;
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

  int idx = msa_file.find_last_of("."); 
  std::string prefix = msa_file.substr(0, idx);

  std::cout << "reading sequences... " << std::flush;
  MSA msa = MSA(msa_file, true, true, 0.8);
  std::cout << "done" << std::endl;

  std::cout << "reading parameters... " << std::flush;
  potts_model params;
  if (compat_mode) {
    params = loadPottsModelCompat(params_file);
  } else {
    params = loadPottsModel(params_file, params_J_file);
  }
  std::cout << "done" << std::endl;

  std::cout << "computing energies... " << std::flush;
  arma::Mat<double> energies = arma::Mat<double>(1, msa.M, arma::fill::zeros);
  double E;
  for (int rep = 0; rep < 1; rep++) {
    for (int seq = 0; seq < msa.M; seq++) {
      E = 0;
      for (int i = 0; i < msa.N; i++) {
        E -= params.h.at(msa.alignment.at(seq, i), i);
        for (int j = i + 1; j < msa.N; j++) {
          E -= params.J.at(i, j).at(msa.alignment.at(seq, i),
                                    msa.alignment.at(seq, j));
        }
      }
      energies.at(rep, seq) = E;
    }
  }
  std::cout << "done" << std::endl;

  std::ofstream output_stream(energy_file);

  for (int rep = 0; rep < 1; rep++) {
    for (int m = 0; m < msa.M; m++) {
      output_stream << energies.at(rep, m) << std::endl;
    }
  }
}