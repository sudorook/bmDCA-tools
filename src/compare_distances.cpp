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

  char c;
  // while ((c = getopt(argc, argv, "i:n:N:o:")) != -1) {
  while ((c = getopt(argc, argv, "i:I:n:N:")) != -1) {
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
  msa1.computeHammingDistances();
  msa1.writeHammingDistances("msa1_distances.txt");
  std::cout << "done" << std::endl;

  std::cout << "reading msa 2 sequences... " << std::flush;
  MSA msa2 = MSA(infile2, "", false, is_numeric2, 0.8);
  msa2.computeHammingDistances();
  msa2.writeHammingDistances("msa2_distances.txt");
  std::cout << "done" << std::endl;

  int M1 = msa1.M;
  int M2 = msa2.M;
  
  int N1 = msa1.N;
  int N2 = msa2.N;

  if (N1 != N2) {
    std::cerr << "ERROR: alignments are of different length." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::Col<double> msa1_msa2_distances =
    arma::Col<double>(M1, arma::fill::zeros);
  arma::Col<double> msa2_msa1_distances =
    arma::Col<double>(M2, arma::fill::zeros);
  arma::Mat<int> alignment1_T = msa1.alignment.t();
  arma::Mat<int> alignment2_T = msa2.alignment.t();

  arma::wall_clock timer;
  timer.tic();
  std::cout << "computing distances from msa 1 to msa 2... " << std::flush;
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

        if (id > msa1_msa2_distances(m1)) {
          msa1_msa2_distances(m1) = id;
        }
      }
    }
  }
  std::cout << timer.toc() << " sec" << std::endl;

  std::cout << "computing distances from msa 2 to msa 1... " << std::flush;
  {
    int* m1_ptr = nullptr;
    int* m2_ptr = nullptr;
    int count = 0;
    double id = 0;
    for (int m2 = 0; m2 < M2; m2++) {
      m2_ptr = alignment2_T.colptr(m2);
      for (int m1 = 0; m1 < M1; m1++) {
        count = 0;
        m1_ptr = alignment1_T.colptr(m1);
        for (int n = 0; n < N1; n++) {
          if (*(m1_ptr + n) == *(m2_ptr + n)) {
            count++;
          }
        }
        id = (double)count / N1;

        if (id > msa2_msa1_distances(m2)) {
          msa2_msa1_distances(m2) = id;
        }
      }
    }
  }
  std::cout << timer.toc() << " sec" << std::endl;

  std::cout << "writing output... " << std::flush;
  {
    std::ofstream output_stream("msa1_msa2_distances.txt");
    for (int i = 0; i < M1; i++) {
      output_stream << msa1_msa2_distances(i) << std::endl;
    }
    output_stream.close();
  }
  {
    std::ofstream output_stream("msa2_msa1_distances.txt");
    for (int i = 0; i < M2; i++) {
      output_stream << msa2_msa1_distances(i) << std::endl;
    }
    output_stream.close();
  }
  std::cout << "done" << std::endl;
}
