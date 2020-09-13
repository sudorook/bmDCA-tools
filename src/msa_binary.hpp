#ifndef MSA_BINARY_HPP
#define MSA_BINARY_HPP

#include <armadillo>
#include <string>
#include <vector>

#include "utils.hpp"

class MSABin
{
public:
  arma::Mat<unsigned int> alignment;
  arma::Col<double> sequence_weights;
  int M;                              // number of sequences
  int N;                              // number of positions
  int Q;                              // number of amino acids

  MSABin(std::string, bool, std::vector<int>);
  MSABin(arma::Mat<unsigned int>, int, int, int, std::vector<int>);

  void writeNumericAlignment(std::string);
  void writeSequenceWeights(std::string);

  // void filterGaps(double = 0.2, double = 0.2);
  void filterSequenceGaps(double = 0.2);
  void filterPositionGaps(double = 0.2);
  void filterSimilarSequences(double = 0.8);
  void computeSequenceWeights(double = 0.8);

private:
  std::vector<SeqRecord> seq_records;
  int getSequenceLength(std::string);
  void readInputMSA(std::string);
  void readInputNumericMSA(std::string);
  void makeNumericalMatrix(void);

  std::vector<int> keep_seq;
};

#endif
