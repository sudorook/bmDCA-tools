#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa_binary.hpp"

int
main(int argc, char* argv[])
{
  std::string infile;
  bool is_numeric = false;
  double reweighting_threshold = 0.8;
  double similarity_threshold = 1.0;
  double position_gap_threshold = 0.2;
  double sequence_gap_threshold = 0.2;
  std::vector<int> keep_seq;

  std::string outfile = "";
  std::string outfile_weights = "";

  char c;
  while ((c = getopt(argc, argv, "i:n:k:g:G:s:t:o:O:")) != -1) {
    switch (c) {
      case 'i':
        infile = optarg;
        is_numeric = false;
        break;
      case 'n':
        infile = optarg;
        is_numeric = true;
        break;
      case 'k':
        keep_seq.push_back(std::stoi(optarg));
        break;
      case 't':
        reweighting_threshold = std::stod(optarg);
        break;
      case 's':
        similarity_threshold = std::stod(optarg);
        break;
      case 'g':
        position_gap_threshold = std::stod(optarg);
        break;
      case 'G':
        sequence_gap_threshold = std::stod(optarg);
        break;
      case 'o':
        outfile = optarg;
        break;
      case 'O':
        outfile = optarg;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  std::cout << "reading sequences... " << std::flush;
  MSABin msa = MSABin(infile, is_numeric, keep_seq);
  std::cout << "done" << std::endl;
  
  // std::cout << "filtering gaps... " << std::flush;
  // msa.filterGaps(sequence_gap_threshold, position_gap_threshold);
  // std::cout << "done" << std::endl;
  
  if ((position_gap_threshold < 1) & (position_gap_threshold >= 0)) {
    std::cout << "filter gapped positions... " << std::flush;
    msa.filterPositionGaps(position_gap_threshold);
    std::cout << "done" << std::endl;
  }

  if ((sequence_gap_threshold < 1) & (sequence_gap_threshold >= 0)) {
    std::cout << "filter gapped sequences... " << std::flush;
    msa.filterSequenceGaps(sequence_gap_threshold);
    std::cout << "done" << std::endl;
  }
 
  if ((similarity_threshold < 1) & (similarity_threshold >= 0)) {
    std::cout << "filtering similar sequences... " << std::flush;
    msa.filterSimilarSequences(similarity_threshold);
    std::cout << "done" << std::endl;
  }
  
  std::cout << "computing sequence weights... " << std::flush;
  msa.computeSequenceWeights(reweighting_threshold);
  std::cout << "done" << std::endl;

  if (outfile.size() == 0) {
    int idx = infile.find_last_of(".");
    std::string outfile_prefix = infile.substr(0, idx);
    outfile = outfile_prefix + "_processed.txt";
  }
  if (outfile_weights.size() == 0) {
    if (outfile.size() != 0) {
      int idx = outfile.find_last_of(".");
      std::string outfile_prefix = outfile.substr(0, idx);
      outfile_weights = outfile_prefix + "_weights.txt";
    } else {
      int idx = infile.find_last_of(".");
      std::string outfile_prefix = infile.substr(0, idx);
      outfile_weights = outfile_prefix + "_processed_weights.txt";
    }
  }

  msa.writeNumericAlignment(outfile);
  msa.writeSequenceWeights(outfile_weights);

  return 0;
};
