#ifndef MSA_STATS_HPP
#define MSA_STATS_HPP

#include "msa.hpp"

#include <armadillo>

class MSAStats
{
public:
  MSAStats(MSA);
  arma::Mat<double> GetRelEntropyGradient();
  arma::Mat<double> GetFrequency1p();
  double GetEffectiveM();
  double GetL();
  double GetM();
  double GetQ();
  void WriteRelEntropyGradientCompat(std::string);
  void WriteFrequency1pCompat(std::string);
  void WriteFrequency2pCompat(std::string);
  void WriteFrequency3pCompat(std::string);
  
  void WriteCorrelation2pCompat(std::string);
  void WriteCorrelation3pCompat(std::string);

  arma::Mat<double> frequency_1p;
  arma::field<arma::Mat<double>> frequency_2p;
  arma::field<arma::Cube<double>> frequency_3p;
  // arma::field<arma::Mat<double>> correlation_2p;
  // arma::field<arma::Cube<double>> correlation_3p;
  arma::Mat<double> rel_entropy_grad_1p;

private:
  double pseudocount;
  int M;              // number of sequences
  int L;              // number of positions
  int Q;              // amino acid alphabet size
  double M_effective; // effect number of sequences

  arma::Col<double> aa_background_frequencies;
  // arma::Mat<double> frequency_1p;
  // arma::field<arma::Mat<double>> frequency_2p;
  // arma::Mat<double> rel_entropy_grad_1p;
};

#endif
