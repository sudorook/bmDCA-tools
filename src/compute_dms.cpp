#include <armadillo>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>

#include "msa.hpp"
#include "utils.hpp"

char
convertAA(int n)
{
  char aa = '\0';
  switch (n) {
    case 0:
      aa = '-';
      break;
    case 1:
      aa = 'A';
      break;
    case 2:
      aa = 'C';
      break;
    case 3:
      aa = 'D';
      break;
    case 4:
      aa = 'E';
      break;
    case 5:
      aa = 'F';
      break;
    case 6:
      aa = 'G';
      break;
    case 7:
      aa = 'H';
      break;
    case 8:
      aa = 'I';
      break;
    case 9:
      aa = 'K';
      break;
    case 10:
      aa = 'L';
      break;
    case 11:
      aa = 'M';
      break;
    case 12:
      aa = 'N';
      break;
    case 13:
      aa = 'P';
      break;
    case 14:
      aa = 'Q';
      break;
    case 15:
      aa = 'R';
      break;
    case 16:
      aa = 'S';
      break;
    case 17:
      aa = 'T';
      break;
    case 18:
      aa = 'V';
      break;
    case 19:
      aa = 'W';
      break;
    case 20:
      aa = 'Y';
      break;
  }
  return aa;
};


void
writeDMS(arma::Mat<double> dms, std::string filename) {
  std::ofstream output_stream(filename);
 
  int Q = dms.n_rows;
  int N = dms.n_cols;
  
  for (int i = 0; i < N; i++) {
    output_stream << "\t" << i;
  }
  output_stream << std::endl;

  for (int aa = 0; aa < Q; aa++) {
    output_stream << convertAA(aa);
    for (int i = 0; i < N; i++) {
      output_stream << "\t" << dms(aa, i);
    }
    output_stream << std::endl;
  }
  output_stream.close();

  return;
};


int
main(int argc, char* argv[])
{
  std::string msa_file;
  std::string weight_file;
  // std::string energy_file = "muence_energies.txt";
  std::string params_file;
  std::string params_J_file;
  bool compat_mode = true;
  bool is_numeric = false;

  int reference_idx = -1;

  char c;
  while ((c = getopt(argc, argv, "i:o:n:p:P:r:")) != -1) {
    switch (c) {
      case 'i':
        msa_file = optarg;
        break;
      case 'n':
        msa_file = optarg;
        is_numeric = true;
        break;
      case 'w':
        weight_file = optarg;
        break;
      case 'p':
        params_file = optarg;
        break;
      case 'P':
        params_J_file = optarg;
        compat_mode = false;
        break;
      case 'r':
        reference_idx = std::stoi(optarg);
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  int idx = msa_file.find_last_of(".");
  std::string prefix = msa_file.substr(0, idx);

  std::cout << "reading sequences... " << std::flush;
  MSA msa = MSA(msa_file, weight_file, false, is_numeric, 1.0);
  std::cout << "done" << std::endl;

  std::cout << "reading parameters... " << std::flush;
  potts_model params;
  if (compat_mode) {
    params = loadPottsModelCompat(params_file);
  } else {
    params = loadPottsModel(params_file, params_J_file);
  }
  std::cout << "done" << std::endl;

  int M = msa.M;
  int N = msa.N;
  int Q = msa.Q;

  // Compute initial neighborhood.
  std::cout << "computing dms... " << std::flush;
  arma::Cube<double> de = arma::Cube<double>(Q, N, M, arma::fill::zeros);
  for (int m = 0; m < M; m++) {
    for (int i = 0; i < N; i++) {
      int q0 = msa.alignment(m, i);
      double e0 = -params.h(q0, i);
      for (int j = 0; j < N; ++j) {
        if (i > j) {
          e0 -= params.J(j, i)(msa.alignment(m, j), q0);
        } else if (i < j) {
          e0 -= params.J(i, j)(q0, msa.alignment(m, j));
        }
      }

      for (int q1 = 0; q1 < Q; q1++) {
        if (q0 != q1) {
          double e1 = -params.h(q1, i);
          for (int j = 0; j < N; ++j) {
            if (i > j) {
              e1 -= params.J(j, i)(msa.alignment(m, j), q1);
            } else if (i < j) {
              e1 -= params.J(i, j)(q1, msa.alignment(m, j));
            }
          }
          de(q1, i, m) = e1 - e0;
        }
      }
    }
  }
  std::cout << "done" << std::endl;

  std::cout << "writing output... " << std::flush;
  if (reference_idx != -1) {
    arma::Mat<double> ref_dms = -de.slice(reference_idx);
    writeDMS(ref_dms,
             prefix + "_dms_" + std::to_string(reference_idx) + ".tsv");
  }

  arma::Mat<double> mean_dms = arma::mean(-de, 2);
  writeDMS(mean_dms, prefix + "_dms_mean.tsv");
  
  arma::Mat<double> min_dms = arma::min(-de, 2);
  writeDMS(min_dms, prefix + "_dms_min.tsv");

  arma::Mat<double> max_dms = arma::max(-de, 2);
  writeDMS(max_dms, prefix + "_dms_max.tsv");
  std::cout << "done" << std::endl;

  return 0;
}
