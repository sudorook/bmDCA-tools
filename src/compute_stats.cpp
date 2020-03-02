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
  std::string outfile;
  std::string dest;
  // bool reweight = true;

  char c;
  while ((c = getopt(argc, argv, "i:o:")) != -1) {
    switch (c) {
      case 'i':
        infile = optarg;
        break;
      case 'o':
        outfile = optarg;
        break;
      // case 'd':
      //   dest = optarg;
      //   break;
      // case 'r':
      //   reweight = true;
      //   break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  std::cout << "reading sequences" << std::endl;
  MSA msa = MSA(infile, true, true, 0.8);

  int idx = infile.find_last_of("."); 
  std::string prefix = infile.substr(0, idx);
  msa.writeSequenceWeights(prefix + "_weights.txt");
  
  std::cout << "computing stats" << std::endl;
  MSAStats msa_stats = MSAStats(msa);
  
  std::cout << "writing 1p stats" << std::endl;
  msa_stats.writeFrequency1p(outfile + "_freq_1p.txt");
  
  std::cout << "writing 2p stats" << std::endl;
  msa_stats.writeCorrelation2p(outfile + "_corr_2p.txt");
  
  std::cout << "writing 3p stats" << std::endl;
  // msa_stats.writeCorrelation3p(outfile + "_corr_3p.txt", 0.01);
  msa_stats.writeCorrelation3p(outfile + "_corr_3p.txt");
}
