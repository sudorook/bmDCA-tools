#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.hpp"
#include "msa_stats.hpp"

#define BINS 201

typedef struct {
  arma::Mat<unsigned long long int> grid;
  int bins=BINS;
  double bin_width;
  double min=0;
  double max=1;
} histogram;

typedef struct {
  double a;
  double b;
  double R2;
} linear_model;

void
writeHistogram(std::string file, histogram hist)
{
  std::ofstream output_stream(file);
  for (int i = 0; i < BINS; i++) {
    for (int j = 0; j < BINS; j++) {
      output_stream << hist.min + i * hist.bin_width << "\t"
                    << hist.min + j * hist.bin_width << "\t"
                    << hist.grid.at(i, j) << std::endl;
    }
  }
  output_stream.close();
  return;
};

void
writeLinearModel(std::string file, linear_model model)
{
  std::ofstream output_stream(file);
  output_stream << "a"
                << "\t" << model.a << std::endl;
  output_stream << "b"
                << "\t" << model.b << std::endl;
  output_stream << "R2"
                << "\t" << model.R2 << std::endl;
  output_stream.close();
  return;
};

int main(int argc, char* argv[]) {
  std::string msa_file;
  std::string mc_file;
  double threshold = 0.8;

  bool reweight = false;

  char c;
  while ((c = getopt(argc, argv, "s:c:t:r")) != -1) {
    switch (c) {
      case 's':
        msa_file = optarg;
        break;
      case 'c':
        mc_file = optarg;
        break;
      case 'r':
        reweight = true;
        break;
      case 't':
        threshold = std::stod(optarg);
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  std::cout << "reading sequences" << std::endl;
  MSA msa = MSA(msa_file, reweight, true, threshold);
  MSA mc = MSA(mc_file, false, true, threshold);

  int idx = msa_file.find_last_of(".");
  std::string msa_prefix = msa_file.substr(0, idx);
  idx = mc_file.find_last_of(".");
  std::string mc_prefix = mc_file.substr(0, idx);

  std::cout << "computing stats" << std::endl;
  MSAStats msa_stats = MSAStats(msa);
  MSAStats mc_stats = MSAStats(mc);
  std::cout << std::endl;

  int N = msa_stats.getN();
  int Q = msa_stats.getQ();
  if ( (N != mc_stats.getN()) | (Q != mc_stats.getQ()) ) {
    std::cerr << "Mismatched number of positions and/or alphabet size"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int M_msa = msa_stats.getM();
  int M_mc = mc_stats.getM();
  double M_effective_msa = msa_stats.getEffectiveM();
  double M_effective_mc = mc_stats.getEffectiveM();

  arma::wall_clock timer;

  std::cout << "computing 1p frequency regression... " << std::flush;
  timer.tic();
  {
    double f = 1. / (double)Q;
    double xy = arma::accu((msa_stats.frequency_1p - f) %
                           (mc_stats.frequency_1p - f));
    double xx = arma::accu(arma::pow((msa_stats.frequency_1p - f), 2));
    double yy = arma::accu(arma::pow((mc_stats.frequency_1p - f), 2));

    double b = xy / xx;
    double a = f*(1-b);
    double r2 = (2*b*xy - b*b*xx)/yy;
    std::cout << timer.toc() << " sec" << std::endl;
    std::cout << "model: " << b << " x + " << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;

    linear_model model;
    model.a = a;
    model.b = b;
    model.R2 = r2;

    histogram hist;
    hist.grid =
      arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
    double bin_width = 1. / (double)(BINS - 1);
    hist.bin_width = bin_width;
    for (int i = 0; i < Q * N; i++) {
      hist.grid(floor(msa_stats.frequency_1p(i) / bin_width),
           floor(mc_stats.frequency_1p(i) / bin_width))++;
    }
    // hist.grid.print("freq 1p:");
    writeHistogram("freq_1p_hist.tsv", hist);
    writeLinearModel("freq_1p_model.tsv", model);
    std::cout << std::endl;
  }

  std::cout << "computing 2p frequency regression... " << std::flush;
  timer.tic();
  {
    double f = 1. / (double)(Q*Q);
    double xy = arma::accu((msa_stats.frequency_2p - f) %
                           (mc_stats.frequency_2p - f));
    double xx = arma::accu(arma::pow((msa_stats.frequency_2p - f), 2));
    double yy = arma::accu(arma::pow((mc_stats.frequency_2p - f), 2));

    double b = xy / xx;
    double a = f*(1-b);
    double r2 = (2*b*xy - b*b*xx)/yy;
    std::cout << timer.toc() << " sec" << std::endl;
    std::cout << "model: " << b << " x + " << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")"
              << std::endl;

    linear_model model;
    model.a = a;
    model.b = b;
    model.R2 = r2;

    histogram hist;
    hist.grid =
      arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
    double bin_width = 1. / (double)(BINS - 1);
    hist.bin_width = bin_width;
    for (int i = 0; i < N * (N - 1) / 2 * Q * Q; i++) {
      hist.grid(floor(msa_stats.frequency_2p(i) / bin_width),
           floor(mc_stats.frequency_2p(i) / bin_width))++;
    }
    // hist.grid.print("freq 2p:");
    writeHistogram("freq_2p_hist.tsv", hist);
    writeLinearModel("freq_2p_model.tsv", model);
    std::cout << std::endl;
  }

  std::cout << "computing 2p correlation regression... " << std::flush;
  timer.tic();
  {
    double xy = arma::accu((msa_stats.correlation_2p) %
                           (mc_stats.correlation_2p));
    double xx = arma::accu(arma::pow((msa_stats.correlation_2p), 2));
    double yy = arma::accu(arma::pow((mc_stats.correlation_2p), 2));

    double b = xy / xx;
    double a = 0;
    double r2 = (2*b*xy - b*b*xx)/yy;
    std::cout << timer.toc() << " sec" << std::endl;
    std::cout << "model: " << b << " x + " << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;

    linear_model model;
    model.a = a;
    model.b = b;
    model.R2 = r2;

    histogram hist;
    hist.grid =
      arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
    double max = std::max(arma::max(msa_stats.correlation_2p),
                     arma::max(mc_stats.correlation_2p));
    double min = std::min(arma::min(msa_stats.correlation_2p),
                     arma::min(mc_stats.correlation_2p));
    double bin_width = 1.05 * (max - min) / (double)(BINS - 1);
    hist.bin_width = bin_width;
    hist.max = 1.05 * max;
    hist.min = 1.05 * min;
    for (int i = 0; i < N * (N - 1) / 2 * Q * Q; i++) {
      hist.grid(floor((msa_stats.correlation_2p(i) - hist.min) / bin_width + .5),
                floor((mc_stats.correlation_2p(i) - hist.min) / bin_width + .5))++;
    }
    // hist.grid.print("corr 2p:");
    writeHistogram("corr_2p_hist.tsv", hist);
    writeLinearModel("corr_2p_model.tsv", model);
    std::cout << std::endl;
  }

  std::cout << "computing 3p frequency and correlation regressions... "
            << std::flush;
  timer.tic();
  {
    arma::Col<double> freq_xx_vec =
      arma::Col<double>((int)(N * (N - 1) * (N - 2) / 6), arma::fill::zeros);
    arma::Col<double> freq_xy_vec =
      arma::Col<double>((int)(N * (N - 1) * (N - 2) / 6), arma::fill::zeros);
    arma::Col<double> freq_yy_vec =
      arma::Col<double>((int)(N * (N - 1) * (N - 2) / 6), arma::fill::zeros);

    arma::Col<double> corr_xx_vec =
      arma::Col<double>((int)(N * (N - 1) * (N - 2) / 6), arma::fill::zeros);
    arma::Col<double> corr_xy_vec =
      arma::Col<double>((int)(N * (N - 1) * (N - 2) / 6), arma::fill::zeros);
    arma::Col<double> corr_yy_vec =
      arma::Col<double>((int)(N * (N - 1) * (N - 2) / 6), arma::fill::zeros);

    double f = 1. / (double)(Q*Q*Q);

    double freq_bin_width = 1. / (double)(BINS - 1);
    double max = std::max(arma::max(msa_stats.correlation_2p),
                     arma::max(mc_stats.correlation_2p));
    double min = std::min(arma::min(msa_stats.correlation_2p),
                     arma::min(mc_stats.correlation_2p));
    double corr_bin_width = (max - min) / (double)(BINS - 1);

    arma::Cube<unsigned long long int> freq_hist_i =
      arma::Cube<unsigned long long int>(BINS, BINS, N, arma::fill::zeros);
    arma::Cube<unsigned long long int> corr_hist_i =
      arma::Cube<unsigned long long int>(BINS, BINS, N, arma::fill::zeros);

    histogram freq_hist;
    freq_hist.grid =
      arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
    freq_hist.bin_width = freq_bin_width;
    freq_hist.bins = BINS;

    histogram corr_hist;
    corr_hist.grid =
      arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
    corr_hist.bin_width = corr_bin_width;
    corr_hist.bins = BINS;
    corr_hist.max = max;
    corr_hist.min = min;

#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < N; i++) {

        double* msa_weight_ptr = msa.sequence_weights.memptr();
        double* mc_weight_ptr = mc.sequence_weights.memptr();

        for (int j = i + 1; j < N; j++) {
          for (int k = j + 1; k < N; k++) {

            int idx =
              (N * (N - 1) * (N - 2) / 6 -
               (N - i - 1) * (N - i) * (N - i - 2) / 6 + k - j +
               (j - i - 1) * (N - i - 1) - (j - i) * (j - i - 1) / 2 - 1);

            arma::Col<double> msa_frequency_3p_ijk = arma::Col<double>(Q*Q*Q,
                arma::fill::zeros);
            arma::Col<double> msa_correlation_3p_ijk = arma::Col<double>(Q*Q*Q,
                arma::fill::zeros);
            int* msa_align_ptr1 = msa.alignment.colptr(i);
            int* msa_align_ptr2 = msa.alignment.colptr(j);
            int* msa_align_ptr3 = msa.alignment.colptr(k);

            arma::Col<double> mc_frequency_3p_ijk = arma::Col<double>(Q*Q*Q,
                arma::fill::zeros);
            arma::Col<double> mc_correlation_3p_ijk = arma::Col<double>(Q*Q*Q,
                arma::fill::zeros);
            int* mc_align_ptr1 = mc.alignment.colptr(i);
            int* mc_align_ptr2 = mc.alignment.colptr(j);
            int* mc_align_ptr3 = mc.alignment.colptr(k);

            for (int m = 0; m < M_msa; m++) {
              int msa_idx = *(msa_align_ptr3 + m) + *(msa_align_ptr2 + m) * Q +
                            *(msa_align_ptr1 + m) * Q * Q;

              msa_frequency_3p_ijk.at(msa_idx) += *(msa_weight_ptr + m);
            }
            msa_frequency_3p_ijk = msa_frequency_3p_ijk / M_effective_msa;

            for (int m = 0; m < M_mc; m++) {
              int mc_idx = *(mc_align_ptr3 + m) + *(mc_align_ptr2 + m) * Q +
                           *(mc_align_ptr1 + m) * Q * Q;

              mc_frequency_3p_ijk.at(mc_idx) += *(mc_weight_ptr + m);
            }
            mc_frequency_3p_ijk = mc_frequency_3p_ijk / M_effective_mc;

            freq_xx_vec.at(idx) =
              arma::accu(arma::pow(msa_frequency_3p_ijk - f, 2));
            freq_yy_vec.at(idx) =
              arma::accu(arma::pow(mc_frequency_3p_ijk - f, 2));
            freq_xy_vec.at(idx) = arma::accu((msa_frequency_3p_ijk - f) %
                                             (mc_frequency_3p_ijk - f));

            for (int aa1 = 0; aa1 < Q; aa1++) {
              for (int aa2 = 0; aa2 < Q; aa2++) {
                for (int aa3 = 0; aa3 < Q; aa3++) {
                  int idx_aa = aa3 + aa2 * Q + aa1 * Q * Q;
                  int idx2_12 =
                    aa2 + aa1 * Q +
                    Q * Q * (j - i - 1 + i * (N - 1) - i * (i - 1) / 2);
                  int idx2_13 =
                    aa3 + aa1 * Q +
                    Q * Q * (k - i - 1 + i * (N - 1) - i * (i - 1) / 2);
                  int idx2_23 =
                    aa3 + aa2 * Q +
                    Q * Q * (k - j - 1 + j * (N - 1) - j * (j - 1) / 2);

                  msa_correlation_3p_ijk.at(idx_aa) =
                    msa_frequency_3p_ijk.at(idx_aa) -
                    msa_stats.frequency_2p.at(idx2_12) *
                      msa_stats.frequency_1p.at(aa3 + k * Q) -
                    msa_stats.frequency_2p.at(idx2_13) *
                      msa_stats.frequency_1p.at(aa2 + j * Q) -
                    msa_stats.frequency_2p.at(idx2_23) *
                      msa_stats.frequency_1p.at(aa1 + i * Q) +
                    2.0 * msa_stats.frequency_1p.at(aa1 + i * Q) *
                      msa_stats.frequency_1p.at(aa2 + j * Q) *
                      msa_stats.frequency_1p.at(aa3 + k * Q);

                  mc_correlation_3p_ijk.at(idx_aa) =
                    mc_frequency_3p_ijk.at(idx_aa) -
                    mc_stats.frequency_2p.at(idx2_12) *
                      mc_stats.frequency_1p.at(aa3 + k * Q) -
                    mc_stats.frequency_2p.at(idx2_13) *
                      mc_stats.frequency_1p.at(aa2 + j * Q) -
                    mc_stats.frequency_2p.at(idx2_23) *
                      mc_stats.frequency_1p.at(aa1 + i * Q) +
                    2.0 * mc_stats.frequency_1p.at(aa1 + i * Q) *
                      mc_stats.frequency_1p.at(aa2 + j * Q) *
                      mc_stats.frequency_1p.at(aa3 + k * Q);
                }
              }
            }

            {
              arma::Col<double> tmp_msa_frequency_3p_ijk =
                arma::floor(msa_frequency_3p_ijk / freq_bin_width);
              arma::Col<double> tmp_mc_frequency_3p_ijk =
                arma::floor(mc_frequency_3p_ijk / freq_bin_width);

              arma::Col<double> tmp_msa_correlation_3p_ijk = arma::floor(
                (msa_correlation_3p_ijk - min) / corr_bin_width + .5);
              arma::Col<double> tmp_mc_correlation_3p_ijk = arma::floor(
                (mc_correlation_3p_ijk - min) / corr_bin_width + .5);

              for (int q = 0; q < Q * Q * Q; q++) {
                freq_hist_i(tmp_msa_frequency_3p_ijk.at(q),
                            tmp_mc_frequency_3p_ijk.at(q),
                            i)++;
                corr_hist_i(tmp_msa_correlation_3p_ijk.at(q),
                            tmp_mc_correlation_3p_ijk.at(q),
                            i)++;
              }
            }

            corr_xx_vec.at(idx) =
              arma::accu(arma::pow(msa_correlation_3p_ijk, 2));
            corr_yy_vec.at(idx) =
              arma::accu(arma::pow(mc_correlation_3p_ijk, 2));
            corr_xy_vec.at(idx) =
              arma::accu((msa_correlation_3p_ijk) % (mc_correlation_3p_ijk));
          }
        }
      }
    }
    std::cout << timer.toc() << " sec" << std::endl;

    {
      double xx = arma::accu(freq_xx_vec);
      double xy = arma::accu(freq_xy_vec);
      double yy = arma::accu(freq_yy_vec);

      double b = xy / xx;
      double a = f*(1-b);
      double r2 = (2*b*xy - b*b*xx)/yy;
      std::cout << "model: " << b << " x + " << a << ", R^2=" << r2
                << " (r=" << sqrt(r2) << ")" << std::endl;

      linear_model model;
      model.a = a;
      model.b = b;
      model.R2 = r2;

      freq_hist.grid = arma::sum(freq_hist_i, 2);
      // freq_hist.grid.print("freq 3p:");
      writeHistogram("freq_3p_hist.tsv", freq_hist);
      writeLinearModel("freq_3p_model.tsv", model);
      std::cout << std::endl;
    }
    {
      double xx = arma::accu(corr_xx_vec);
      double xy = arma::accu(corr_xy_vec);
      double yy = arma::accu(corr_yy_vec);

      double b = xy / xx;
      double a = 0;
      double r2 = (2*b*xy - b*b*xx)/yy;
      std::cout << "model: " << b << " x + " << a << ", R^2=" << r2
                << " (r=" << sqrt(r2) << ")" << std::endl;

      linear_model model;
      model.a = a;
      model.b = b;
      model.R2 = r2;

      corr_hist.grid = arma::sum(corr_hist_i, 2);
      // corr_hist.grid.print("corr 3p:");
      writeHistogram("corr_3p_hist.tsv", corr_hist);
      writeLinearModel("corr_3p_model.tsv", model);
    }
  }
};
