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

#include "msa.hpp"

int
main(int argc, char* argv[])
{
  std::string infile1;
  std::string infile2;
  std::string outfile;
  bool is_numeric1 = false;
  bool is_numeric2 = false;
  double between_threshold = 0.8;
  double within_threshold = 0.8;

  char c;
  while ((c = getopt(argc, argv, "i:I:n:N:o:w:b:")) != -1) {
    switch (c) {
      case 'i':
        infile1 = optarg;
        break;
      case 'I':
        infile2 = optarg;
        break;
      case 'n':
        infile1 = optarg;
        is_numeric1 = true;
        break;
      case 'N':
        infile2 = optarg;
        is_numeric2 = true;
        break;
      case 'o':
        outfile = optarg;
        break;
      case 'w':
        within_threshold = std::stod(optarg);
        break;
      case 'b':
        between_threshold = std::stod(optarg);
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  std::cout << "reading msa 1 sequences... " << std::flush;
  MSA msa1 = MSA(infile1, "", false, is_numeric1, 0.8);
  std::cout << "done" << std::endl;

  std::cout << "reading msa 2 sequences... " << std::flush;
  MSA msa2 = MSA(infile2, "", false, is_numeric2, 0.8);
  std::cout << "done" << std::endl;

  int M1 = msa1.M;
  int M2 = msa2.M;

  int Q1 = msa1.Q;
  int Q2 = msa2.Q;

  int N1 = msa1.N;
  int N2 = msa2.N;

  if (N1 != N2) {
    std::cerr << "ERROR: alignments are of different length." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (Q1 != Q2) {
    std::cerr << "ERROR: alignments have different number of states."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::Mat<int> alignment1_T = msa1.alignment.t();
  arma::Mat<int> alignment2_T = msa2.alignment.t();
  arma::Mat<int> results = arma::Mat<int>(M1, N1, arma::fill::zeros);
  int result_counter = 0;

  arma::wall_clock timer;
  timer.tic();
  std::cout << "finding distant sequences... " << std::flush;
  {
    int* m1_ptr = nullptr;
    int* m2_ptr = nullptr;
    int count = 0;
    double id = 0;
    for (int m1 = 0; m1 < M1; m1++) {
      m1_ptr = alignment1_T.colptr(m1);

      volatile bool too_close = false;
      if ((between_threshold < 1) & (between_threshold >= 0)) {
#pragma omp parallel for shared(too_close)
        for (int m2 = 0; m2 < M2; m2++) {
          if (too_close)
            continue;
          count = 0;
          m2_ptr = alignment2_T.colptr(m2);
          for (int n = 0; n < N1; n++) {
            if (*(m1_ptr + n) == *(m2_ptr + n)) {
              count++;
            }
          }
          id = (double)count / N1;
          if (id > between_threshold) {
            too_close = true;
          }
        }
      }

      if ((within_threshold < 1) & (within_threshold >= 0)) {
        if (too_close == false) {
#pragma omp parallel for shared(too_close)
          for (int m2 = 0; m2 < result_counter; m2++) {
            if (too_close)
              continue;
            arma::Mat<int> results_T = results.t();
            count = 0;
            m2_ptr = results_T.colptr(m2);
            for (int n = 0; n < N1; n++) {
              if (*(m1_ptr + n) == *(m2_ptr + n)) {
                count++;
              }
            }
            id = (double)count / N1;
            if (id > within_threshold) {
              too_close = true;
            }
          }
        }
      }

      if (too_close == false) {
        for (int n = 0; n < N1; n++) {
          results(result_counter, n) = *(m1_ptr + n);
        }
        result_counter++;
      }
    }
  }
  std::cout << timer.toc() << " sec" << std::endl;

  std::cout << "writing output... " << std::flush;
  {
    std::ofstream output_stream(outfile);
    output_stream << result_counter << " " << N1 << " " << Q1 << std::endl;
    for (int m = 0; m < result_counter; m++) {
      output_stream << results(m, 0);
      for (int n = 1; n < N1; n++) {
        output_stream << " " << results(m, n);
      }
      output_stream << std::endl;
    }
    output_stream.close();
  }

  // {
  //   std::ofstream output_stream(outfile);
  //   char aa;
  //   for (int m = 0; m < result_counter; m++) {
  //     output_stream << ">sample" << m << std::endl;
  //     for (int n = 0; n < N1; n++) {
  //       aa = convertAA(results(m, n));
  //       if (aa != '\0') {
  //         output_stream << aa;
  //       }
  //     }
  //     output_stream << std::endl << std::endl;
  //   }
  //   output_stream.close();
  // }

  std::cout << "done" << std::endl;
}
