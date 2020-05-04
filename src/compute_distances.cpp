#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.hpp"

int main(int argc, char* argv[]) {
  std::string infile;
  std::string outfile;
  std::string dest;
  bool is_numeric = false;
  // bool reweight = true;

  char c;
  while ((c = getopt(argc, argv, "i:n:o:")) != -1) {
    switch (c) {
      case 'i':
        infile = optarg;
        break;
      case 'n':
        infile = optarg;
        is_numeric = true;
        break;
      case 'o':
        outfile = optarg;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  std::cout << "reading sequences" << std::endl;
  MSA msa = MSA(infile, "", false, is_numeric, 0.8);

  msa.writeHammingDistances(outfile);
}
