#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.hpp"
#include "utils.hpp"
#include "pcg_random.hpp"

int
find_the_fucker(potts_model* params,
                arma::Mat<double>* bias_field,
                arma::Col<int> start,
                arma::Col<int> stop,
                arma::Mat<int>* sequences,
                arma::Col<int>* mutations_before,
                arma::Col<int>* mutations_after,
                arma::Col<int>* positions,
                arma::Col<double>* energies,
                size_t max_iter,
                size_t m,
                size_t n,
                size_t q,
                long int seed,
                double bias,
                double energy_start,
                double temperature)
{
  pcg32 rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> uniform(0, 1);

  int distance = 0;
  arma::Col<size_t> conf = arma::Col<size_t>(n);
  for (size_t i = 0; i < n; ++i) {
    conf(i) = start(i);
    if (start(i) != stop(i)) {
      distance += 1;
    }
  }

  std::cout << "start distance: " << distance << std::endl;
  double E = energy_start;
  double dE = 0;
  size_t step = 0;
  while ((distance > 0) & (step < max_iter)) {
    size_t i = size_t(n * uniform(rng));
    size_t dq = 1 + size_t((q - 1) * uniform(rng));

    size_t q0 = conf(i);
    size_t q1 = (q0 + dq) % q;

    double e0 = -params->h(q0, i);
    double e1 = -params->h(q1, i);
    for (size_t j = 0; j < n; ++j) {
      if (i > j) {
        e0 -= params->J(j, i)(conf(j), q0);
        e1 -= params->J(j, i)(conf(j), q1);
      } else if (i < j) {
        e0 -= params->J(i, j)(q0, conf(j));
        e1 -= params->J(i, j)(q1, conf(j));
      }
    }
    dE = e1 - e0;
    double de_bias = bias_field->at(q0, i) - bias_field->at(q1, i);
    // double de = dE * (1 - bias) + bias * de_bias;
    double de = dE + bias * de_bias;
    if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
      conf(i) = q1;
      E += dE;
      if (stop(i) == (int)q1) {
        distance -= 1;
      }
      if (stop(i) == (int)q0) {
        distance += 1;
      }
      for (size_t i = 0; i < n; ++i) {
        (*sequences)(step, i) = conf(i);
      }
      (*mutations_before)(step) = q0;
      (*mutations_after)(step) = q1;
      (*energies)(step) = E;
      (*positions)(step) = i+1;
      // distance = 0;
      // for (size_t i = 0; i < n; ++i) {
      //   if ((int)conf(i) != stop(i)) {
      //     distance += 1;
      //   }
      // }
      step++;
    }
    if (step % 10000) {
      std::cout << "distance: " << distance << std::endl;
    }
  }
  std::cout << "found the fucker in " << step << " steps" << std::endl;
  return step;
};

char
convertAA(int n)
{ 
  char aa = '\0';
  switch (n) {
    case 0:
      aa = '-';
      break;
    case 1:
      aa = 'A';
      break;
    case 2:
      aa = 'C';
      break;
    case 3:
      aa = 'D';
      break;
    case 4:
      aa = 'E';
      break;
    case 5:
      aa = 'F';
      break;
    case 6:
      aa = 'G';
      break;
    case 7:
      aa = 'H';
      break;
    case 8:
      aa = 'I';
      break;
    case 9:
      aa = 'K';
      break;
    case 10:
      aa = 'L';
      break;
    case 11:
      aa = 'M';
      break;
    case 12:
      aa = 'N';
      break;
    case 13:
      aa = 'P';
      break;
    case 14:
      aa = 'Q';
      break;
    case 15:
      aa = 'R';
      break;
    case 16:
      aa = 'S';
      break;
    case 17:
      aa = 'T';
      break;
    case 18:
      aa = 'V';
      break;
    case 19:
      aa = 'W';
      break;
    case 20:
      aa = 'Y';
      break;
  }
  return aa;
};
  


int
main(int argc, char* argv[])
{
  std::string msa_start_file;
  std::string msa_stop_file;
  std::string params_file;
  std::string params_J_file;
  double bias = 0;
  bool compat_mode = true;
  bool is_start_numeric = false;
  bool is_stop_numeric = false;

  char c;
  while ((c = getopt(argc, argv, "b:i:I:n:N:p:P:")) != -1) {
    switch (c) {
      case 'i':
        msa_start_file = optarg;
        break;
      case 'n':
        msa_start_file = optarg;
        is_start_numeric = true;
        break;
      case 'I':
        msa_stop_file = optarg;
        break;
      case 'N':
        msa_stop_file = optarg;
        is_stop_numeric = true;
        break;
      case 'p':
        params_file = optarg;
        break;
      case 'b':
        bias = std::stod(optarg);
        break;
      case 'P':
        params_J_file = optarg;
        compat_mode = false;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  std::cout << "reading start sequence... " << std::flush;
  MSA msa_start = MSA(msa_start_file, "", false, is_start_numeric, 0.8);
  std::cout << "done" << std::endl;
  arma::Col<int> seq_start = msa_start.alignment.row(0);
  seq_start.print("start:");

  std::cout << "reading end sequence... " << std::flush;
  MSA msa_stop = MSA(msa_stop_file, "", false, is_stop_numeric, 0.8);
  std::cout << "done" << std::endl;
  arma::Col<int> seq_stop = msa_stop.alignment.row(0);
  seq_stop.print("stop:");

  int M = msa_start.M;
  int Q = msa_start.Q;
  int N = msa_start.N;

  std::cout << "M: " << M << std::endl;
  std::cout << "Q: " << Q << std::endl;
  std::cout << "N: " << N << std::endl;

  std::cout << "reading parameters... " << std::flush;
  potts_model params;
  if (compat_mode) {
    params = loadPottsModelCompat(params_file);
  } else {
    params = loadPottsModel(params_file, params_J_file);
  }
  std::cout << "done" << std::endl;

  std::cout << "computing energies... " << std::flush;
  double E_start = 0;
  double E_stop = 0;
  for (int i = 0; i < N; i++) {
    E_start -= params.h.at(seq_start(i), i);
    E_stop -= params.h.at(seq_stop(i), i);
    for (int j = i + 1; j < N; j++) {
      E_start -= params.J.at(i, j).at(seq_start(i), seq_start(j));
      E_stop -= params.J.at(i, j).at(seq_stop(i), seq_stop(j));
    }
  }
  std::cout << "done" << std::endl;

  std::cout << "starting energy: " << E_start << std::endl;
  std::cout << "stopping energy: " << E_stop << std::endl;

  std::cout << "computing field bias... " << std::flush;
  arma::Mat<double> bias_field = arma::Mat<double>(Q, N, arma::fill::zeros);
  double pseudocount = 0.003;
  for (int i = 0; i < N; i++) {
    double avg = 0;
    for (int aa = 0; aa < Q; aa++) {
      if (seq_stop(i) == aa) {
        avg += log((1. - pseudocount) * 1 + pseudocount * (1. / Q));
      } else {
        avg += log((1. - pseudocount) * 0 + pseudocount * (1. / Q));
      }
    }
    for (int aa = 0; aa < Q; aa++) {
      if (seq_stop(i) == aa) {
        bias_field(aa, i) =
          log((1. - pseudocount) * 1 + pseudocount * (1. / Q)) - avg / Q;
      } else {
        bias_field(aa, i) =
          log((1. - pseudocount) * 0 + pseudocount * (1. / Q)) - avg / Q;
      }
    }
  }
  std::cout << "done" << std::endl;

  double max_iter = 100000;
  arma::Mat<int> sequences = arma::Mat<int>(max_iter, N, arma::fill::zeros);
  arma::Col<int> mutations_before = arma::Col<int>(max_iter, arma::fill::zeros);
  arma::Col<int> mutations_after = arma::Col<int>(max_iter, arma::fill::zeros);
  arma::Col<int> positions = arma::Col<int>(max_iter, arma::fill::zeros);
  arma::Col<double> energies = arma::Col<double>(max_iter, arma::fill::zeros);

  std::cout << "searching..." << std::endl;
  int steps = find_the_fucker(&params,
                              &bias_field,
                              seq_start,
                              seq_stop,
                              &sequences,
                              &mutations_before,
                              &mutations_after,
                              &positions,
                              &energies,
                              max_iter,
                              M,
                              N,
                              Q,
                              10,
                              bias,
                              E_start,
                              1.0);

  std::cout << "writing output... " << std::flush;
  {
    std::ofstream output_stream("path_numerical.txt");
    output_stream << steps << " " << N << " " << Q << std::endl;
    for (int i = 0; i < steps; i++) {
      for (int j = 0; j < N; j++) {
        if (j + 1 == N) {
          output_stream << sequences(i, j) << std::endl;
        } else {
          output_stream << sequences(i, j) << " ";
        }
      }
    }
    output_stream.close();
  }
  {
    std::ofstream output_stream("path_energies.txt");
    for (int i = 0; i < steps; i++) {
      output_stream << energies(i) << std::endl;
    }
    output_stream.close();
  }
  {
    std::ofstream output_stream("path_mutations_before.txt");
    for (int i = 0; i < steps; i++) {
      output_stream << convertAA(mutations_before(i)) << std::endl;
    }
    output_stream.close();
  }
  {
    std::ofstream output_stream("path_mutations_after.txt");
    for (int i = 0; i < steps; i++) {
      output_stream << convertAA(mutations_after(i)) << std::endl;
    }
    output_stream.close();
  }
  {
    std::ofstream output_stream("path_positions.txt");
    for (int i = 0; i < steps; i++) {
      output_stream << positions(i) << std::endl;
    }
    output_stream.close();
  }
  std::cout << "done" << std::endl;
}
