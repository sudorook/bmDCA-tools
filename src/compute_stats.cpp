#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.hpp"
#include "msa_stats.hpp"

int main(int argc, char* argv[]) {
  std::string infile;
  bool reweight = false;
  bool is_numeric = false;
  double threshold = 0.8;

  char c;
  while ((c = getopt(argc, argv, "i:n:rt:")) != -1) {
    switch (c) {
      case 'i':
        infile = optarg;
        is_numeric = false;
        break;
      case 'n':
        infile = optarg;
        is_numeric = true;
        break;
      case 'r':
        reweight = true;
        break;
      case 't':
        threshold = std::stod(optarg);
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  std::cout << "reading sequences" << std::endl;
  MSA msa = MSA(infile, reweight, is_numeric, threshold);

  int idx = infile.find_last_of("."); 
  std::string prefix = infile.substr(0, idx);
  msa.writeSequenceWeights(prefix + "_weights.txt");
  
  std::cout << "computing stats" << std::endl;
  MSAStats msa_stats = MSAStats(msa);
  
  std::cout << "writing 1p stats" << std::endl;
  msa_stats.writeFrequency1p(prefix + "_freq_1p.bin");
  
  std::cout << "writing 2p stats" << std::endl;
  msa_stats.writeFrequency2p(prefix + "_freq_2p.bin");
  msa_stats.writeCorrelation2p(prefix + "_corr_2p.bin");
  
  // std::cout << "writing 3p stats" << std::endl;
  // msa_stats.writeFrequency3p(prefix + "_freq_3p.bin");
  // msa_stats.writeFrequency3pAscii(prefix + "_freq_3p.txt");
  // msa_stats.writeCorrelation3p(prefix + "_corr_3p.bin");
}
