#ifndef MSA_HPP
#define MSA_HPP

#include <armadillo>
#include <string>
#include <vector>

class MSA
{
public:
  arma::Mat<int> Alignment;
  arma::Col<double> SequenceWeights;
  int SequenceCount;
  int SequenceLength;

  MSA(std::string, bool = true, double = 0.8);
  void WriteMatrixCompat(std::string);
  void WriteSequenceWeightsCompat(std::string);

private:
  void ReadInputMSA(std::string);
  void ComputeSequenceWeights(double);
};

#endif
