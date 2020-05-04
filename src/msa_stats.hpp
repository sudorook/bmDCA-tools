#ifndef MSA_STATS_HPP
#define MSA_STATS_HPP

#include "msa.hpp"

#include <armadillo>

class MSAStats
{
public:
  MSAStats(MSA);
  double getEffectiveM();
  double getN();
  double getM();
  double getQ();

  void writeFrequency1p(std::string);
  void writeFrequency2p(std::string);
  void writeCorrelation2p(std::string);

  void writeFrequency1pAscii(std::string);
  void writeFrequency2pAscii(std::string);
  void writeCorrelation2pAscii(std::string);

  arma::Col<double> frequency_1p;
  arma::Col<double> frequency_2p;
  arma::Col<double> correlation_2p;

private:
  int M;              // number of sequences
  int N;              // number of positions
  int Q;              // amino acid alphabet size
  double M_effective; // effect number of sequences

  MSA msa;

  void computeFrequency1p(void);
  void computeFrequency2p(void);
  void computeFrequency3p(void);
  void computeCorrelation2p(void);
  void computeCorrelation3p(void);
  // void computeCorrelation4p(void);
};

#endif
