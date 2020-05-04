#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.hpp"
#include "msa_compare.hpp"

int main(int argc, char* argv[]) {
  std::string msa1_file;
  std::string msa2_file;
  std::string msa1_weight_file = "";
  std::string msa2_weight_file = "";
  double threshold = 0.8;
  int bins = 201;

  bool numeric_msa1 = false;
  bool numeric_msa2 = false;
  bool reweight1 = false;
  bool reweight2 = false;

  int positions = 2;

  char c;
  while ((c = getopt(argc, argv, "i:I:n:N:w:W:t:rRb:hp:")) != -1) {
    switch (c) {
      case 'b':
        bins = std::stoi(optarg);
        break;
      case 'p':
        positions = std::stoi(optarg);
        break;
      case 'i':
        msa1_file = optarg;
        break;
      case 'I':
        msa2_file = optarg;
        break;
      case 'n':
        msa1_file = optarg;
        numeric_msa1 = true;
        break;
      case 'N':
        msa2_file = optarg;
        numeric_msa2 = true;
        break;
      case 'w':
        msa1_weight_file = optarg;
        break;
      case 'W':
        msa2_weight_file = optarg;
        break;
      case 'r':
        reweight1 = true;
        break;
      case 'R':
        reweight2 = true;
        break;
      case 't':
        threshold = std::stod(optarg);
        break;
      case 'h':
        std::cerr << "go fuck yourself" << std::endl;
        std::exit(EXIT_FAILURE);
      case '?':
        std::cerr << "what the fuck?" << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  MSA msa = MSA(msa1_file, msa1_weight_file, reweight1, numeric_msa1, threshold);
  MSA mc = MSA(msa2_file, msa2_weight_file, reweight2, numeric_msa2, threshold);

  MSACompare stats = MSACompare(&msa, &mc, bins);

  if (positions >= 1) {
    stats.computeFrequency1p();
    stats.makeFrequency1pHistogram();
  }

  if (positions >= 2) {
    stats.computeFrequency2p();
    // stats.computeCorrelation2p();
    stats.makeFrequency2pHistogram();
    stats.makeCorrelation2pHistogram();
  }

  if (positions == 3) {
    stats.makeEfficient3pHistograms();
  } else if (positions > 3) {
    stats.makeEfficient3pHistograms();
    stats.computeFrequency3p();
    // stats.computeCorrelation3p();
    // stats.makeFrequency3pHistogram();
    // stats.makeCorrelation3pHistogram();
  }

  if (positions == 4) {
    stats.makeEfficient4pHistograms();
  } else if (positions > 4) {
    std::cout << "well aren't you ambitious..." << std::endl;
    stats.makeEfficient4pHistograms();
  }
};
