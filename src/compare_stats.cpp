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
  std::string msa_file;
  std::string mc_file;
  double threshold = 0;

  char c;
  while ((c = getopt(argc, argv, "s:c:t:")) != -1) {
    switch (c) {
      case 's':
        msa_file = optarg;
        break;
      case 'c':
        mc_file = optarg;
        break;
      case 't':
        threshold = std::stod(optarg);
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  std::cout << "reading sequences" << std::endl;
  MSA msa = MSA(msa_file, true, true, 0.8);
  MSA mc = MSA(mc_file, true, true, 0.8);

  int idx = msa_file.find_last_of("."); 
  std::string msa_prefix = msa_file.substr(0, idx);
  idx = mc_file.find_last_of(".");
  std::string mc_prefix = mc_file.substr(0, idx);
  
  std::cout << "computing stats" << std::endl;
  MSAStats msa_stats = MSAStats(msa);
  MSAStats mc_stats = MSAStats(mc);
  
  std::cout << "writing 1p stats" << std::endl;
  msa_stats.writeFrequency1p(msa_prefix + "_freq_1p.txt");
  mc_stats.writeFrequency1p(mc_prefix + "_freq_1p.txt");
  
  std::cout << "writing 2p stats" << std::endl;
  msa_stats.writeCorrelation2p(msa_prefix + "_corr_2p.txt");
  mc_stats.writeCorrelation2p(mc_prefix + "_corr_2p.txt");
  
  std::cout << "writing 3p stats" << std::endl;
  {
    std::ofstream msa_stream(msa_prefix + "_corr_3p.txt");
    std::ofstream mc_stream(mc_prefix + "_corr_3p.txt");

    int N = msa.N;
    int Q = msa.Q;

    double msa_tmp = 0;
    double mc_tmp = 0;
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int k = j + 1; k < N; k++) {
          for (int aa1 = 0; aa1 < Q; aa1++) {
            for (int aa2 = 0; aa2 < Q; aa2++) {
              for (int aa3 = 0; aa3 < Q; aa3++) {
                msa_tmp = msa_stats.frequency_3p.at(i, j, k).at(aa1, aa2, aa3) -
                          msa_stats.frequency_2p.at(i, j).at(aa1, aa2) *
                            msa_stats.frequency_1p.at(aa3, k) -
                          msa_stats.frequency_2p.at(i, k).at(aa1, aa3) *
                            msa_stats.frequency_1p.at(aa2, j) -
                          msa_stats.frequency_2p.at(j, k).at(aa2, aa3) *
                            msa_stats.frequency_1p.at(aa1, i) +
                          2.0 * msa_stats.frequency_1p.at(aa1, i) *
                            msa_stats.frequency_1p.at(aa2, j) *
                            msa_stats.frequency_1p.at(aa3, k);
                mc_tmp = mc_stats.frequency_3p.at(i, j, k).at(aa1, aa2, aa3) -
                         mc_stats.frequency_2p.at(i, j).at(aa1, aa2) *
                           mc_stats.frequency_1p.at(aa3, k) -
                         mc_stats.frequency_2p.at(i, k).at(aa1, aa3) *
                           mc_stats.frequency_1p.at(aa2, j) -
                         mc_stats.frequency_2p.at(j, k).at(aa2, aa3) *
                           mc_stats.frequency_1p.at(aa1, i) +
                         2.0 * mc_stats.frequency_1p.at(aa1, i) *
                           mc_stats.frequency_1p.at(aa2, j) *
                           mc_stats.frequency_1p.at(aa3, k);
                if ((fabs(msa_tmp) > threshold) || (fabs(mc_tmp) > threshold)) {
                  msa_stream << msa_tmp << std::endl;
                  mc_stream << mc_tmp << std::endl;
                }
              }
            }
          }
        }
      }
    }
  }
};
