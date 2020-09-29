#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.hpp"
#include "utils.hpp"

int
main(int argc, char* argv[])
{
  std::string msa_file;
  std::string energy_file = "sequence_energies.txt";
  std::string params_file;
  std::string params_J_file;
  bool compat_mode = true;
  bool is_numeric = false;
  bool ignore_gaps = false;

  char c;
  while ((c = getopt(argc, argv, "i:o:n:p:P:g")) != -1) {
    switch (c) {
      case 'i':
        msa_file = optarg;
        break;
      case 'n':
        msa_file = optarg;
        is_numeric = true;
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
      case 'g':
        ignore_gaps = true;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
        std::exit(EXIT_FAILURE);
        break;
    }
  }

  int idx = msa_file.find_last_of(".");
  std::string prefix = msa_file.substr(0, idx);

  std::cout << "reading sequences... " << std::flush;
  MSA msa = MSA(msa_file, "", false, is_numeric, 0.8);
  std::cout << "done" << std::endl;

  std::cout << "reading parameters... " << std::flush;
  potts_model params;
  if (compat_mode) {
    params = loadPottsModelAscii(params_file);
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
        if ((ignore_gaps) & (msa.alignment(seq, i) == 0))
          continue;
        E -= params.h(msa.alignment(seq, i), i);
        for (int j = i + 1; j < msa.N; j++) {
          if ((ignore_gaps) & (msa.alignment(seq, j) == 0))
            continue;
          E -= params.J(i, j)(msa.alignment(seq, i), msa.alignment(seq, j));
        }
      }
      energies(rep, seq) = E;
    }
  }
  std::cout << "done" << std::endl;

  std::cout << "writing energies... " << std::flush;
  std::ofstream output_stream(energy_file);

  // output_stream.unsetf(std::ios::floatfield);
  // output_stream.precision(50);

  for (int rep = 0; rep < 1; rep++) {
    for (int m = 0; m < msa.M; m++) {
      output_stream << energies(rep, m) << std::endl;
    }
  }
  output_stream.close();
  std::cout << "done" << std::endl;
}
