#include "msa_stats.hpp"

#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>

MSAStats::MSAStats(MSA msa): msa(msa)
{
  // Initialize
  N = msa.N;
  M = msa.M;
  Q = msa.Q;
  M_effective = sum(msa.sequence_weights);

  std::cout << M << " sequences" << std::endl;
  std::cout << N << " positions" << std::endl;
  std::cout << Q << " amino acids (including gaps)" << std::endl;
  std::cout << M_effective << " effective sequences" << std::endl;

  computeFrequency1p();
  computeFrequency2p();
  // computeCorrelation2p();
};

void MSAStats::computeFrequency1p(void) {
  arma::wall_clock timer;
  frequency_1p = arma::Col<double>((int)Q * N, arma::fill::zeros);

  // std::cout << "computing 1p frequencies... " << std::flush;
  // timer.tic();
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      int* align_ptr = msa.alignment.colptr(i);
      double* weight_ptr = msa.sequence_weights.memptr();
      for (int m = 0; m < M; m++) {
        frequency_1p.at(i * Q + *(align_ptr + m)) += *(weight_ptr + m);
      }
    }
  }
  frequency_1p = frequency_1p / M_effective;
  // std::cout << timer.toc() << " sec" << std::endl;
};

void MSAStats::computeFrequency2p(void) {
  arma::wall_clock timer;
  frequency_2p =
    arma::Col<double>((int)N * (N - 1) / 2 * Q * Q, arma::fill::zeros);

  // std::cout << "computing 2p frequencies... " << std::flush;
  // timer.tic();
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        double* weight_ptr = msa.sequence_weights.memptr();
        int* align_ptr1 = msa.alignment.colptr(i);
        int* align_ptr2 = msa.alignment.colptr(j);
        for (int m = 0; m < M; m++) {
          frequency_2p.at(*(align_ptr2 + m) + *(align_ptr1 + m) * Q +
                          Q * Q *
                            (j - i - 1 + i * (N - 1) - i * (i - 1) / 2)) +=
            *(weight_ptr + m);
        }
      }
    }
  }
  frequency_2p = frequency_2p / M_effective;
  // std::cout << timer.toc() << " sec" << std::endl;
};

void MSAStats::computeCorrelation2p(void) {
  arma::wall_clock timer;
  correlation_2p = arma::Col<double>((int)N*(N-1)/2*Q*Q, arma::fill::zeros);

  // std::cout << "computing 2p correlations... " << std::flush;
  // timer.tic();
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int aa1 = 0; aa1 < Q; aa1++) {
          for (int aa2 = 0; aa2 < Q; aa2++) {
            int idx = aa2 + aa1 * Q +
                      Q * Q * (j - i - 1 + i * (N - 1) - i * (i - 1) / 2);
            correlation_2p.at(idx) =
              frequency_2p.at(idx) -
              frequency_1p.at(i * Q + aa1) * frequency_1p.at(j * Q + aa2);
          }
        }
      }
    }
  }
  // std::cout << timer.toc() << " sec" << std::endl;
};

double
MSAStats::getQ(void)
{
  return Q;
};

double
MSAStats::getM(void)
{
  return M;
};

double
MSAStats::getN(void)
{
  return N;
};

double
MSAStats::getEffectiveM(void)
{
  return M_effective;
};

void
MSAStats::writeFrequency1p(std::string output_file)
{
  frequency_1p.save(output_file, arma::raw_binary);
};

void
MSAStats::writeFrequency1pAscii(std::string output_file)
{
  frequency_1p.save(output_file, arma::raw_ascii);
};

void
MSAStats::writeFrequency2p(std::string output_file)
{
  frequency_2p.save(output_file, arma::raw_binary);
};

void
MSAStats::writeFrequency2pAscii(std::string output_file)
{
  frequency_2p.save(output_file, arma::raw_ascii);
};

void
MSAStats::writeFrequency3p(std::string output_file)
{
  frequency_3p.save(output_file, arma::raw_binary);
};

void
MSAStats::writeFrequency3pAscii(std::string output_file)
{
  frequency_3p.save(output_file, arma::raw_ascii);
};

void
MSAStats::writeCorrelation2p(std::string output_file)
{
  correlation_2p.save(output_file, arma::raw_binary);
};

void
MSAStats::writeCorrelation2pAscii(std::string output_file)
{
  correlation_2p.save(output_file, arma::raw_ascii);
};

// void
// MSAStats::writeCorrelation3p(std::string output_file)
// {
//   correlation_3p = arma::field<arma::Cube<double>>(N, N, N);
//   for (int i = 0; i < N; i++) {
//     for (int j = i + 1; j < N; j++) {
//       for (int k = j + 1; k < N; k++) {
//         for (int aa1 = 0; aa1 < Q; aa1++) {
//           for (int aa2 = 0; aa2 < Q; aa2++) {
//             for (int aa3 = 0; aa3 < Q; aa3++) {
//               correlation_3p.at(i, j, k).at(aa1, aa2, aa3) =
//                 frequency_3p.at(i, j, k).at(aa1, aa2, aa3) -
//                 frequency_2p.at(i, j).at(aa1, aa2) * frequency_1p.at(aa3, k) -
//                 frequency_2p.at(i, k).at(aa1, aa3) * frequency_1p.at(aa2, j) -
//                 frequency_2p.at(j, k).at(aa2, aa3) * frequency_1p.at(aa1, i) +
//                 2.0 * frequency_1p.at(aa1, i) * frequency_1p.at(aa2, j) *
//                   frequency_1p.at(aa3, k);
//             }
//           }
//         }
//       }
//     }
//   }
//   correlation_3p.save(output_file, arma::raw_binary);
// };
//
// void
// MSAStats::writeCorrelation3pAscii(std::string output_file)
// {
//   std::ofstream output_stream(output_file);
//
//   for (int i = 0; i < N; i++) {
//     for (int j = i + 1; j < N; j++) {
//       for (int k = j + 1; k < N; k++) {
//         output_stream << i << " " << j << " " << k;
//         for (int aa1 = 0; aa1 < Q; aa1++) {
//           for (int aa2 = 0; aa2 < Q; aa2++) {
//             for (int aa3 = 0; aa3 < Q; aa3++) {
//               output_stream << frequency_3p.at(i, j, k).at(aa1, aa2, aa3) -
//                 frequency_2p.at(i, j).at(aa1, aa2)*frequency_1p.at(aa3, k) -
//                 frequency_2p.at(i, k).at(aa1, aa3)*frequency_1p.at(aa2, j) -
//                 frequency_2p.at(j, k).at(aa2, aa3)*frequency_1p.at(aa1, i) +
//                 2.0*frequency_1p.at(aa1, i)*frequency_1p.at(aa2,
//                     j)*frequency_1p.at(aa3, k);
//             }
//           }
//         }
//         output_stream << std::endl;;
//       }
//     }
//   }
// };
