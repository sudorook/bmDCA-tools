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
  // std::string outfile;
  bool is_numeric1 = false;
  bool is_numeric2 = false;
  int reference_idx1 = -1;
  int reference_idx2 = -1;

  char c;
  while ((c = getopt(argc, argv, "i:I:n:N:r:R:")) != -1) {
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
      case 'r':
        reference_idx1 = std::stoi(optarg);
        break;
      case 'R':
        reference_idx2 = std::stoi(optarg);
        break;
      // case 'o':
      //   outfile = optarg;
      //   break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  std::cout << "reading msa 1 sequences... " << std::flush;
  MSA msa1 = MSA(infile1, "", false, is_numeric1, 0.8);
  msa1.computeSequenceSimilarity();
  msa1.writeSequenceSimilarity("msa1_max_similarity.txt",
                               "msa1_mean_similarity.txt");
  std::cout << "done" << std::endl;

  std::cout << "reading msa 2 sequences... " << std::flush;
  MSA msa2 = MSA(infile2, "", false, is_numeric2, 0.8);
  msa2.computeSequenceSimilarity();
  msa2.writeSequenceSimilarity("msa2_max_similarity.txt",
                               "msa2_mean_similarity.txt");
  std::cout << "done" << std::endl;

  int M1 = msa1.M;
  int M2 = msa2.M;

  int N1 = msa1.N;
  int N2 = msa2.N;

  if (N1 != N2) {
    std::cerr << "ERROR: alignments are of different length." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::Mat<int> alignment1_T = msa1.alignment.t();
  arma::Mat<int> alignment2_T = msa2.alignment.t();

  arma::Col<double> msa1_msa2_max_distances =
    arma::Col<double>(M1, arma::fill::zeros);
  arma::Col<double> msa2_msa1_max_distances =
    arma::Col<double>(M2, arma::fill::zeros);

  arma::Col<double> msa1_msa2_mean_distances =
    arma::Col<double>(M1, arma::fill::zeros);
  arma::Col<double> msa2_msa1_mean_distances =
    arma::Col<double>(M2, arma::fill::zeros);

  arma::Col<double> msa1_ref1_max_distances =
    arma::Col<double>(M1, arma::fill::zeros);
  arma::Col<double> msa2_ref1_max_distances =
    arma::Col<double>(M2, arma::fill::zeros);

  arma::Col<double> msa1_ref2_max_distances =
    arma::Col<double>(M1, arma::fill::zeros);
  arma::Col<double> msa2_ref2_max_distances =
    arma::Col<double>(M2, arma::fill::zeros);

  arma::Col<double> msa1_ref1_mean_distances =
    arma::Col<double>(M1, arma::fill::zeros);
  arma::Col<double> msa2_ref1_mean_distances =
    arma::Col<double>(M2, arma::fill::zeros);

  arma::Col<double> msa1_ref2_mean_distances =
    arma::Col<double>(M1, arma::fill::zeros);
  arma::Col<double> msa2_ref2_mean_distances =
    arma::Col<double>(M2, arma::fill::zeros);

  arma::wall_clock timer;
  timer.tic();
  std::cout << "computing distances between msa 1 and msa 2... " << std::flush;
  {
    int* m1_ptr = nullptr;
    int* m2_ptr = nullptr;
    int count = 0;
    double id = 0;
    for (int m1 = 0; m1 < M1; m1++) {
      m1_ptr = alignment1_T.colptr(m1);
      for (int m2 = 0; m2 < M2; m2++) {
        count = 0;
        m2_ptr = alignment2_T.colptr(m2);
        for (int n = 0; n < N1; n++) {
          if (*(m1_ptr + n) == *(m2_ptr + n)) {
            count++;
          }
        }
        id = (double)count / N1;

        msa1_msa2_mean_distances(m1) += id;
        msa2_msa1_mean_distances(m2) += id;

        if (id > msa1_msa2_max_distances(m1)) {
          msa1_msa2_max_distances(m1) = id;
        }
        if (id > msa2_msa1_max_distances(m2)) {
          msa2_msa1_max_distances(m2) = id;
        }
      }
    }
    msa1_msa2_mean_distances = msa1_msa2_mean_distances / M2;
    msa2_msa1_mean_distances = msa2_msa1_mean_distances / M1;
  }
  std::cout << timer.toc() << " sec" << std::endl;

  if ((reference_idx1 != -1) & (reference_idx1 < M1)) {
    timer.tic();
    std::cout << "computing distances from msa 1 to msa 1 reference sequence "
              << reference_idx1 << " " << std::flush;
    {
      int* m1_ptr = nullptr;
      int* ref_ptr = nullptr;
      int count = 0;
      double id = 0;
      for (int m1 = 0; m1 < M1; m1++) {
        m1_ptr = alignment1_T.colptr(m1);
        count = 0;
        ref_ptr = alignment1_T.colptr(reference_idx1);
        for (int n = 0; n < N1; n++) {
          if (*(m1_ptr + n) == *(ref_ptr + n)) {
            count++;
          }
        }
        id = (double)count / N1;
        msa1_ref1_max_distances(m1) = id;
      }
    }
    std::cout << timer.toc() << " sec" << std::endl;

    timer.tic();
    std::cout << "computing distances from msa 2 to msa 1 reference sequence "
              << reference_idx1 << " " << std::flush;
    {
      int* ref_ptr = nullptr;
      int* m2_ptr = nullptr;
      int count = 0;
      double id = 0;
      for (int m2 = 0; m2 < M2; m2++) {
        m2_ptr = alignment2_T.colptr(m2);
        count = 0;
        ref_ptr = alignment1_T.colptr(reference_idx1);
        for (int n = 0; n < N2; n++) {
          if (*(m2_ptr + n) == *(ref_ptr + n)) {
            count++;
          }
        }
        id = (double)count / N2;
        msa2_ref1_max_distances(m2) = id;
      }
    }
    std::cout << timer.toc() << " sec" << std::endl;
  }

  if ((reference_idx2 != -1) & (reference_idx2 < M2)) {
    timer.tic();
    std::cout << "computing distances from msa 1 to msa 2 reference sequence "
              << reference_idx2 << " " << std::flush;
    {
      int* m1_ptr = nullptr;
      int* ref_ptr = nullptr;
      int count = 0;
      double id = 0;
      for (int m1 = 0; m1 < M1; m1++) {
        m1_ptr = alignment1_T.colptr(m1);
        count = 0;
        ref_ptr = alignment2_T.colptr(reference_idx2);
        for (int n = 0; n < N1; n++) {
          if (*(m1_ptr + n) == *(ref_ptr + n)) {
            count++;
          }
        }
        id = (double)count / N1;
        msa1_ref2_max_distances(m1) = id;
      }
    }
    std::cout << timer.toc() << " sec" << std::endl;

    timer.tic();
    std::cout << "computing distances from msa 2 to msa 2 reference sequence "
              << reference_idx1 << " " << std::flush;
    {
      int* ref_ptr = nullptr;
      int* m2_ptr = nullptr;
      int count = 0;
      double id = 0;
      for (int m2 = 0; m2 < M2; m2++) {
        m2_ptr = alignment2_T.colptr(m2);
        count = 0;
        ref_ptr = alignment2_T.colptr(reference_idx2);
        for (int n = 0; n < N2; n++) {
          if (*(m2_ptr + n) == *(ref_ptr + n)) {
            count++;
          }
        }
        id = (double)count / N2;
        msa2_ref2_max_distances(m2) = id;
      }
    }
    std::cout << timer.toc() << " sec" << std::endl;
  }

  std::cout << "writing output... " << std::flush;
  {
    std::ofstream output_stream("msa1_msa2_max_similarity.txt");
    for (int i = 0; i < M1; i++) {
      output_stream << msa1_msa2_max_distances(i) << std::endl;
    }
    output_stream.close();
  }
  {
    std::ofstream output_stream("msa1_msa2_mean_similarity.txt");
    for (int i = 0; i < M1; i++) {
      output_stream << msa1_msa2_mean_distances(i) << std::endl;
    }
    output_stream.close();
  }
  {
    std::ofstream output_stream("msa2_msa1_max_similarity.txt");
    for (int i = 0; i < M2; i++) {
      output_stream << msa2_msa1_max_distances(i) << std::endl;
    }
    output_stream.close();
  }
  {
    std::ofstream output_stream("msa2_msa1_mean_similarity.txt");
    for (int i = 0; i < M2; i++) {
      output_stream << msa2_msa1_mean_distances(i) << std::endl;
    }
    output_stream.close();
  }

  if ((reference_idx1 != -1) & (reference_idx1 < M1)) {
    {
      std::ofstream output_stream(
        "msa1_msa1-ref" + std::to_string(reference_idx1) + "_max_similarity.txt");
      for (int i = 0; i < M1; i++) {
        output_stream << msa1_ref1_max_distances(i) << std::endl;
      }
      output_stream.close();
    }
    {
      std::ofstream output_stream(
        "msa2_msa1-ref" + std::to_string(reference_idx1) + "_max_similarity.txt");
      for (int i = 0; i < M2; i++) {
        output_stream << msa2_ref1_max_distances(i) << std::endl;
      }
      output_stream.close();
    }
  }
  if ((reference_idx2 != -1) & (reference_idx2 < M2)) {
    {
      std::ofstream output_stream(
        "msa1_msa2-ref" + std::to_string(reference_idx2) + "_max_similarity.txt");
      for (int i = 0; i < M1; i++) {
        output_stream << msa1_ref2_max_distances(i) << std::endl;
      }
      output_stream.close();
    }
    {
      std::ofstream output_stream(
        "msa2_msa2-ref" + std::to_string(reference_idx2) + "_max_similarity.txt");
      for (int i = 0; i < M2; i++) {
        output_stream << msa2_ref2_max_distances(i) << std::endl;
      }
      output_stream.close();
    }
  }
  std::cout << "done" << std::endl;
}
