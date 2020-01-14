#include "msa.hpp"

#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

MSA::MSA(std::string msa_file, bool reweight, double threshold)
{
  ReadInputMSA(msa_file);
  // SequenceCount = SeqRecords.size();
  // SequenceLength = GetSequenceLength(SeqRecords.begin()->GetSequence());
  // MakeNumericalMatrix();

  if (reweight) {
    ComputeSequenceWeights(threshold);
  } else {
    SequenceWeights = arma::vec(SequenceCount, arma::fill::ones);
  }
};

void
MSA::ReadInputMSA(std::string msa_file)
{
  std::ifstream input_stream(msa_file);

  if (!input_stream) {
    std::cerr << "what the fuck" << std::endl;
    exit(2);
  }

  /*
   * Read a FASTA-formatted multiple sequence alignment. Each record from the
   * file is stored as a SeqRecord object and appended to the SeqRecords
   * vector.
   */
  
  int AlphabetSize = 0;
  input_stream >> SequenceCount >> SequenceLength >> AlphabetSize;
  // std::cout<<line<<std::endl;
  // sscanf(line.c_str(), "%d %d %d\n", &SequenceCount, &SequenceLength, &AlphabetSize);
  // std::cout<<SequenceCount<<std::endl;
  // std::cout<<SequenceLength<<std::endl;
  // std::cout<<AlphabetSize<<std::endl;
  Alignment = arma::Mat<int>(SequenceCount, SequenceLength);
  
  int counter = 0;
  int i = 0;
  // std::string tmp;
  std::string line;
  std::getline(input_stream, line);
  while(std::getline(input_stream, line)) {
    std::istringstream iss(line);
    int n;
    i = 0;

    while (iss >> n) {
      // std::cout<<counter<<","<<i<<": "<<n<<std::endl;
      Alignment.at(counter, i) = n;
      i++;
    }
    counter++;
  }
  // while (!input_stream.eof()) {
  //   counter = 0;
  //   i = 0;
  //   // input_stream >> tmp;
  //   std::getline(input_stream, tmp);
  //   for (auto aa = tmp.begin(); i < SequenceLength; aa++) {
  //     if (*aa != ' ') {
  //       // Alignment(counter, i) = (int)std::atoi(*aa);
  //       Alignment(counter, i) = (int)(*aa);
  //       i++;
  //     }
  //   }
  //   counter++;
  // };
  input_stream.close();
  // Alignment.print();
}

// void
// MSA::MakeNumericalMatrix(void)
// {
//   Alignment = arma::Mat<int>(SequenceCount, SequenceLength);
// 
//   int row_idx = 0;
//   for (auto seq = SeqRecords.begin(); seq != SeqRecords.end(); seq++) {
//     std::string sequence = seq->GetSequence();
//     int col_idx = 0;
//     for (auto aa = sequence.begin(); aa != sequence.end(); aa++) {
//       switch (*aa) {
//         case '-':
//         case 'B':
//         case 'J':
//         case 'O':
//         case 'U':
//         case 'X':
//         case 'Z':
//           Alignment(row_idx, col_idx) = 0;
//           col_idx++;
//           break;
//         case 'A':
//           Alignment(row_idx, col_idx) = 1;
//           col_idx++;
//           break;
//         case 'C':
//           Alignment(row_idx, col_idx) = 2;
//           col_idx++;
//           break;
//         case 'D':
//           Alignment(row_idx, col_idx) = 3;
//           col_idx++;
//           break;
//         case 'E':
//           Alignment(row_idx, col_idx) = 4;
//           col_idx++;
//           break;
//         case 'F':
//           Alignment(row_idx, col_idx) = 5;
//           col_idx++;
//           break;
//         case 'G':
//           Alignment(row_idx, col_idx) = 6;
//           col_idx++;
//           break;
//         case 'H':
//           Alignment(row_idx, col_idx) = 7;
//           col_idx++;
//           break;
//         case 'I':
//           Alignment(row_idx, col_idx) = 8;
//           col_idx++;
//           break;
//         case 'K':
//           Alignment(row_idx, col_idx) = 9;
//           col_idx++;
//           break;
//         case 'L':
//           Alignment(row_idx, col_idx) = 10;
//           col_idx++;
//           break;
//         case 'M':
//           Alignment(row_idx, col_idx) = 11;
//           col_idx++;
//           break;
//         case 'N':
//           Alignment(row_idx, col_idx) = 12;
//           col_idx++;
//           break;
//         case 'P':
//           Alignment(row_idx, col_idx) = 13;
//           col_idx++;
//           break;
//         case 'Q':
//           Alignment(row_idx, col_idx) = 14;
//           col_idx++;
//           break;
//         case 'R':
//           Alignment(row_idx, col_idx) = 15;
//           col_idx++;
//           break;
//         case 'S':
//           Alignment(row_idx, col_idx) = 16;
//           col_idx++;
//           break;
//         case 'T':
//           Alignment(row_idx, col_idx) = 17;
//           col_idx++;
//           break;
//         case 'V':
//           Alignment(row_idx, col_idx) = 18;
//           col_idx++;
//           break;
//         case 'W':
//           Alignment(row_idx, col_idx) = 19;
//           col_idx++;
//           break;
//         case 'Y':
//           Alignment(row_idx, col_idx) = 20;
//           col_idx++;
//           break;
//       }
//     }
//     row_idx++;
//   }
// }

void
MSA::WriteMatrixCompat(std::string output_file)
{
  std::ofstream output_stream(output_file);
  output_stream << SequenceCount << " " << SequenceLength << " " << 21
                << std::endl;
  // output_stream << SeqRecords.size() << " " << SequenceLength << " " << 21
  //               << std::endl;
  for (int i = 0; i < SequenceCount; i++) {
    for (int j = 0; j < SequenceLength; j++) {
      if (j + 1 == SequenceLength && i + 1 != SequenceCount) {
        output_stream << Alignment.at(i, j) << std::endl;
      } else {
        output_stream << Alignment.at(i, j) << " ";
      }
    }
  }
  output_stream << std::endl;
}

/* void
 * MSA::SaveAlignment1D(std::string output_file)
 * {} */

/* void
 * MSA::SaveAlignment2D(std::string output_file)
 * {} */

// void
// MSA::PrintAlignment(void)
// {
//   for (std::vector<SeqRecord>::iterator it = SeqRecords.begin();
//        it != SeqRecords.end();
//        ++it) {
//     it->PrintRecord();
//   }
// }

// int
// MSA::GetSequenceLength(std::string sequence)
// {
//   int valid_aa_count = 0;
//   for (std::string::iterator it = sequence.begin(); it != sequence.end();
//        ++it) {
//     switch (*it) {
//       case '-':
//       case 'B':
//       case 'J':
//       case 'O':
//       case 'U':
//       case 'X':
//       case 'Z':
//       case 'A':
//       case 'C':
//       case 'D':
//       case 'E':
//       case 'F':
//       case 'G':
//       case 'H':
//       case 'I':
//       case 'K':
//       case 'L':
//       case 'M':
//       case 'N':
//       case 'P':
//       case 'Q':
//       case 'R':
//       case 'S':
//       case 'T':
//       case 'V':
//       case 'W':
//       case 'Y':
//         valid_aa_count += 1;
//         break;
//     }
//   }
//   return valid_aa_count;
// }

void
MSA::ComputeSequenceWeights(double threshold)
{
  SequenceWeights = arma::vec(SequenceCount, arma::fill::zeros);
  arma::Mat<int> Alignment_T = Alignment.t();

  double id;
  int* m1_ptr = NULL;
  int* m2_ptr = NULL;
  for (int m1 = 0; m1 < SequenceCount; ++m1) {
    SequenceWeights.at(m1) += 1;
    m1_ptr = Alignment_T.colptr(m1);
    for (int m2 = m1 + 1; m2 < SequenceCount; ++m2) {
      m2_ptr = Alignment_T.colptr(m2);
      id = 0;
      for (int i = 0; i < SequenceLength; ++i) {
        if (*(m1_ptr + i) == *(m2_ptr + i)) {
          id += 1;
        }
      }
      if (id > threshold * SequenceLength) {
        SequenceWeights.at(m1) += 1;
        SequenceWeights.at(m2) += 1;
      }
    }
  }

  for (int m1 = 0; m1 < SequenceCount; ++m1) {
    SequenceWeights.at(m1) = 1. / SequenceWeights.at(m1);
  }
}

void
MSA::WriteSequenceWeightsCompat(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < SequenceCount; i++) {
    output_stream << SequenceWeights.at(i) << std::endl;
  }
}
