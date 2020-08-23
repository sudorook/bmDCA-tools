#include <armadillo>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>

#include "msa.hpp"
#include "msa_stats.hpp"
#include "utils.hpp"

int
main(int argc, char* argv[])
{
  std::string msa_file;
  std::string out_file;
  std::string params_file;
  std::string params_J_file;
  bool compat_mode = true;
  bool is_numeric = false;

  int reference_seq = -1;

  char c;
  while ((c = getopt(argc, argv, "i:o:n:p:P:s:")) != -1) {
    switch (c) {
      case 'i':
        msa_file = optarg;
        break;
      case 'o':
        out_file = optarg;
        break;
      case 'n':
        msa_file = optarg;
        is_numeric = true;
        break;
      case 'p':
        params_file = optarg;
        break;
      case 'P':
        params_J_file = optarg;
        compat_mode = false;
        break;
      case 's':
        reference_seq = std::stoi(optarg);
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  int idx = msa_file.find_last_of(".");
  std::string prefix = msa_file.substr(0, idx);

  std::cout << "reading sequences... " << std::flush;
  MSA msa = MSA(msa_file, "", false, is_numeric);
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
#pragma omp parallel
    {
#pragma omp for
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
            de(q1, i, m) = (e1 - e0);
          }
        }
      }
    }
  }
  std::cout << "done" << std::endl;

  std::cout << "computing dcs... " << std::flush;
  
    arma::field<arma::Mat<double>> de_2p =
      arma::field<arma::Mat<double>>(N, N); // upper triangular

  if (reference_seq == -1) {
    arma::field<arma::Mat<double>> de_2p_m =
      arma::field<arma::Mat<double>>(N, N, M); // upper triangular
#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1)
      for (int m  = 0; m < M; m++ ) {
        std::cout << "m = " << m << std::endl;
        for (int n1 = 0; n1 < N; n1++) {
          arma::Mat<double> de_n1 = de.slice(m);
          for (int n2 = n1 + 1; n2 < N; n2++) {
            de_2p_m(n1, n2, m) =
              arma::Mat<double>(Q, Q, arma::fill::zeros); // upper triangular
            for (int aa1 = 0; aa1 < Q ; aa1++) {
              if (aa1 != msa.alignment(m, n1)) {
                // arma::Mat<double> de_n1_step = arma::Mat<double>(Q, N, arma::fill::zeros);
                arma::Mat<double> de_n1_step = de.slice(m);
                int aa_msa = msa.alignment(m, n1);
                
                // Compute neighborhood after moving to aa1
                double tmp = de(m, n1, aa1);
                for (int pos = 0; pos < N; pos++) {
                  for (int aa = 0; aa < Q; aa++) {
                    if (pos < n1) {
                      de_n1_step(pos, aa) +=
                        params.J(pos, n1)(msa.alignment(m, pos), aa1) -
                        params.J(pos, n1)(msa.alignment(m, pos), aa_msa) -
                        params.J(pos, n1)(aa, aa1) + params.J(pos, n1)(aa, aa_msa);
                    } else if (pos > n1) {
                      de_n1_step(pos, aa) +=
                        params.J(n1, pos)(aa1, msa.alignment(m, pos)) -
                        params.J(n1, pos)(aa_msa, msa.alignment(m, pos)) -
                        params.J(n1, pos)(aa1, aa) + params.J(n1, pos)(aa_msa, aa);
                    } else {
                      if (aa1 == aa) {
                        de_n1_step(pos, aa) = 0;
                      } else if (aa_msa == aa) {
                        de_n1_step(pos, aa) = -tmp;
                      } else {
                        de_n1_step(pos, aa) += params.h(aa1, pos) - params.h(aa_msa, pos);
                        for (int pos2 = 0; pos2 < N; pos2++) {
                          if (pos2 < n1) {
                            de_n1_step(pos, aa) +=
                              params.J(pos2, n1)(msa.alignment(m, pos2), aa1) -
                              params.J(pos2, n1)(msa.alignment(m, pos2), aa_msa);
                          } else if (pos2 > n1) {
                            de_n1_step(pos, aa) +=
                              params.J(n1, pos2)(aa1, msa.alignment(m, pos2)) -
                              params.J(n1, pos2)(aa_msa, msa.alignment(m, pos2));
                          }
                        }
                      }
                    }
                  }
                }

                for (int aa2 = aa1 + 1; aa2 < Q ; aa2++) {
                  if (aa2 != msa.alignment(m, n2)) {
                    de_2p_m(n1, n2, m)(aa1, aa2) +=
                      de_n1(n1, aa1) + de_n1(n2, aa2) - de_n1_step(n2, aa2);
                  }
                }
              }
            }
          }
        }
      }
    }

    // average across homologs
    for (int n1 = 0; n1 < N; n1++) {
      for (int n2 = n1 + 1; n2 < N; n2++) {
        de_2p(n1, n2) = arma::Mat<double>(Q, Q, arma::fill::zeros);
        for (int aa1 = 0; aa1 < Q; aa1++) {
          for (int aa2 = aa1 + 1; aa2 < Q; aa2++) {
            for (int m = 0; m < M; m++) {
              de_2p(n1, n2)(aa1, aa2) += de_2p_m(n1, n2, m)(aa1, aa2);
            }
          }
        }
        de_2p(n1, n2) = de_2p(n1, n2) / M;
      }
    }
  } else {
    int m = reference_seq;
#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1)
      for (int n1 = 0; n1 < N; n1++) {
        arma::Mat<double> de_n1 = de.slice(m);
        for (int n2 = n1 + 1; n2 < N; n2++) {
          de_2p(n1, n2) =
            arma::Mat<double>(Q, Q, arma::fill::zeros); // upper triangular
          for (int aa1 = 0; aa1 < Q ; aa1++) {
            if (aa1 != msa.alignment(m, n1)) {
              arma::Mat<double> de_n1_step = de.slice(m);
              int aa_msa = msa.alignment(m, n1);
              
              // Compute neighborhood after moving to aa1
              double tmp = de(m, n1, aa1);
              for (int pos = 0; pos < N; pos++) {
                for (int aa = 0; aa < Q; aa++) {
                  if (pos < n1) {
                    de_n1_step(pos, aa) +=
                      params.J(pos, n1)(msa.alignment(m, pos), aa1) -
                      params.J(pos, n1)(msa.alignment(m, pos), aa_msa) -
                      params.J(pos, n1)(aa, aa1) + params.J(pos, n1)(aa, aa_msa);
                  } else if (pos > n1) {
                    de_n1_step(pos, aa) +=
                      params.J(n1, pos)(aa1, msa.alignment(m, pos)) -
                      params.J(n1, pos)(aa_msa, msa.alignment(m, pos)) -
                      params.J(n1, pos)(aa1, aa) + params.J(n1, pos)(aa_msa, aa);
                  } else {
                    if (aa1 == aa) {
                      de_n1_step(pos, aa) = 0;
                    } else if (aa_msa == aa) {
                      de_n1_step(pos, aa) = -tmp;
                    } else {
                      de_n1_step(pos, aa) += params.h(aa1, pos) - params.h(aa_msa, pos);
                      for (int pos2 = 0; pos2 < N; pos2++) {
                        if (pos2 < n1) {
                          de_n1_step(pos, aa) +=
                            params.J(pos2, n1)(msa.alignment(m, pos2), aa1) -
                            params.J(pos2, n1)(msa.alignment(m, pos2), aa_msa);
                        } else if (pos2 > n1) {
                          de_n1_step(pos, aa) +=
                            params.J(n1, pos2)(aa1, msa.alignment(m, pos2)) -
                            params.J(n1, pos2)(aa_msa, msa.alignment(m, pos2));
                        }
                      }
                    }
                  }
                }
              }

              for (int aa2 = aa1 + 1; aa2 < Q ; aa2++) {
                if (aa2 != msa.alignment(m, n2)) {
                  de_2p(n1, n2)(aa1, aa2) +=
                    de_n1(n1, aa1) + de_n1(n2, aa2) - de_n1_step(n2, aa2);
                }
              }
            }
          }
        }
      }
    }
  }
  std::cout << "done" << std::endl;
  
  std::cout << "writing output... " << std::flush;
  std::ofstream output_stream(out_file);

  output_stream << "N1"
                << "\t"
                << "N2"
                << "\t"
                << "AA1"
                << "\t"
                << "AA2"
                << "\t"
                << "ddG" << std::endl;
  for (int n1 = 0; n1 < N; n1++) {
    for (int n2 = n1 + 1; n2 < N; n2++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = aa1 + 1; aa2 < Q; aa2++) {
          output_stream << n1 << "\t" << n2 << "\t" << aa1 << "\t" << aa2
                        << "\t" << de_2p(n1, n2)(aa1, aa2) << std::endl;
        }
      }
    }
  }
  output_stream.close();
  std::cout << "done" << std::endl;

  return 0;
}
