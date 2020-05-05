#include "msa_compare.hpp"

#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>

MSACompare::MSACompare(MSA* msa, MSA* mc, int bins)
  : msa(msa)
  , mc(mc)
  , bins(bins)
{
  if (msa->Q != mc->Q) {
    std::cerr << "alphabet mismatch: " << msa->Q << " vs " << mc->Q
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (msa->N != mc->N) {
    std::cerr << "position mismatch: " << msa->N << " vs " << mc->N
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Q = msa->Q;
  N = msa->N;
  M_msa = msa->M;
  M_effective_msa = sum(msa->sequence_weights);
  M_mc = mc->M;
  M_effective_mc = sum(mc->sequence_weights);

  std::cout << "reading alignment 1..." << std::endl;
  std::cout << M_msa << " sequences" << std::endl;
  std::cout << N << " positions" << std::endl;
  std::cout << Q << " amino acids (including gaps)" << std::endl;
  std::cout << M_effective_msa << " effective sequences" << std::endl;
  std::cout << std::endl;

  std::cout << "reading alignment 2..." << std::endl;
  std::cout << M_mc << " sequences" << std::endl;
  std::cout << N << " positions" << std::endl;
  std::cout << Q << " amino acids (including gaps)" << std::endl;
  std::cout << M_effective_mc << " effective sequences" << std::endl;
  std::cout << std::endl;
};

void
MSACompare::computeFrequency1p(void)
{
  msa_frequency_1p = arma::Col<double>(Q * N, arma::fill::zeros);
  mc_frequency_1p = arma::Col<double>(Q * N, arma::fill::zeros);

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      int* msa_align_ptr = msa->alignment.colptr(i);
      double* msa_weight_ptr = msa->sequence_weights.memptr();
      for (int m = 0; m < M_msa; m++) {
        msa_frequency_1p(i * Q + *(msa_align_ptr + m)) += *(msa_weight_ptr + m);
      }

      int* mc_align_ptr = mc->alignment.colptr(i);
      double* mc_weight_ptr = mc->sequence_weights.memptr();
      for (int m = 0; m < M_mc; m++) {
        mc_frequency_1p(i * Q + *(mc_align_ptr + m)) += *(mc_weight_ptr + m);
      }
    }
  }
  msa_frequency_1p = msa_frequency_1p / M_effective_msa;
  mc_frequency_1p = mc_frequency_1p / M_effective_mc;

  freq_1p_max =
    std::max(arma::max(msa_frequency_1p), arma::max(mc_frequency_1p));
  freq_1p_min =
    std::min(arma::min(msa_frequency_1p), arma::min(mc_frequency_1p));
  frequency_1p_set = true;
};

void
MSACompare::computeFrequency2p(void)
{

  if (frequency_1p_set == false) {
    std::cerr << "ERROR: compute the 1p frequencies first" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  msa_frequency_2p =
    arma::Col<double>((int)N * (N - 1) / 2 * Q * Q, arma::fill::zeros);
  mc_frequency_2p =
    arma::Col<double>((int)N * (N - 1) / 2 * Q * Q, arma::fill::zeros);

#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        double* msa_weight_ptr = msa->sequence_weights.memptr();
        int* msa_align_ptr1 = msa->alignment.colptr(i);
        int* msa_align_ptr2 = msa->alignment.colptr(j);
        for (int m = 0; m < M_msa; m++) {
          msa_frequency_2p(*(msa_align_ptr2 + m) + *(msa_align_ptr1 + m) * Q +
                           Q * Q *
                             (j - i - 1 + i * (N - 1) - i * (i - 1) / 2)) +=
            *(msa_weight_ptr + m);
        }

        double* mc_weight_ptr = mc->sequence_weights.memptr();
        int* mc_align_ptr1 = mc->alignment.colptr(i);
        int* mc_align_ptr2 = mc->alignment.colptr(j);
        for (int m = 0; m < M_mc; m++) {
          mc_frequency_2p(*(mc_align_ptr2 + m) + *(mc_align_ptr1 + m) * Q +
                          Q * Q *
                            (j - i - 1 + i * (N - 1) - i * (i - 1) / 2)) +=
            *(mc_weight_ptr + m);
        }
      }
    }
  }
  msa_frequency_2p = msa_frequency_2p / M_effective_msa;
  mc_frequency_2p = mc_frequency_2p / M_effective_mc;

  freq_2p_max =
    std::max(arma::max(msa_frequency_2p), arma::max(mc_frequency_2p));
  freq_2p_min =
    std::min(arma::min(msa_frequency_2p), arma::min(mc_frequency_2p));
  frequency_2p_set = true;
};

void
MSACompare::computeFrequency3p(void)
{
  msa_frequency_3p = arma::Col<double>(
    (int)N * (N - 1) * (N - 2) / 6 * Q * Q * Q, arma::fill::zeros);
  mc_frequency_3p = arma::Col<double>(
    (int)N * (N - 1) * (N - 2) / 6 * Q * Q * Q, arma::fill::zeros);

#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int k = j + 1; k < N; k++) {
          double* msa_weight_ptr = msa->sequence_weights.memptr();
          int* msa_align_ptr1 = msa->alignment.colptr(i);
          int* msa_align_ptr2 = msa->alignment.colptr(j);
          int* msa_align_ptr3 = msa->alignment.colptr(k);
          for (int m = 0; m < M_msa; m++) {
            msa_frequency_3p(
              Q * Q * Q *
                (N * (N - 1) * (N - 2) / 6 -
                 (N - i - 1) * (N - i) * (N - i - 2) / 6 + k - j +
                 (j - i - 1) * (N - i - 1) - (j - i) * (j - i - 1) / 2 - 1) +
              *(msa_align_ptr3 + m) + *(msa_align_ptr2 + m) * Q +
              *(msa_align_ptr1 + m) * Q * Q) += *(msa_weight_ptr + m);
          }

          double* mc_weight_ptr = mc->sequence_weights.memptr();
          int* mc_align_ptr1 = mc->alignment.colptr(i);
          int* mc_align_ptr2 = mc->alignment.colptr(j);
          int* mc_align_ptr3 = mc->alignment.colptr(k);
          for (int m = 0; m < M_mc; m++) {
            mc_frequency_3p(
              Q * Q * Q *
                (N * (N - 1) * (N - 2) / 6 -
                 (N - i - 1) * (N - i) * (N - i - 2) / 6 + k - j +
                 (j - i - 1) * (N - i - 1) - (j - i) * (j - i - 1) / 2 - 1) +
              *(mc_align_ptr3 + m) + *(mc_align_ptr2 + m) * Q +
              *(mc_align_ptr1 + m) * Q * Q) += *(mc_weight_ptr + m);
          }
        }
      }
    }
  }
  msa_frequency_3p = msa_frequency_3p / M_effective_msa;
  mc_frequency_3p = mc_frequency_3p / M_effective_mc;

  freq_3p_max =
    std::max(arma::max(msa_frequency_3p), arma::max(mc_frequency_3p));
  freq_3p_min =
    std::min(arma::min(msa_frequency_3p), arma::min(mc_frequency_3p));
  frequency_3p_set = true;
}

void
MSACompare::computeCorrelation2p(void)
{
  if ((frequency_1p_set == false) | (frequency_2p_set == false)) {
    std::cerr << "ERROR: compute the 1p and 2p frequencies first" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  msa_correlation_2p =
    arma::Col<double>((int)N * (N - 1) / 2 * Q * Q, arma::fill::zeros);
  mc_correlation_2p =
    arma::Col<double>((int)N * (N - 1) / 2 * Q * Q, arma::fill::zeros);

#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int aa1 = 0; aa1 < Q; aa1++) {
          for (int aa2 = 0; aa2 < Q; aa2++) {
            int idx = aa2 + aa1 * Q +
                      Q * Q * (j - i - 1 + i * (N - 1) - i * (i - 1) / 2);
            msa_correlation_2p(idx) =
              msa_frequency_2p(idx) -
              msa_frequency_1p(i * Q + aa1) * msa_frequency_1p(j * Q + aa2);
            mc_correlation_2p(idx) =
              mc_frequency_2p(idx) -
              mc_frequency_1p(i * Q + aa1) * mc_frequency_1p(j * Q + aa2);
          }
        }
      }
    }
  }
  corr_2p_max =
    std::max(arma::max(msa_correlation_2p), arma::max(mc_correlation_2p));
  corr_2p_min =
    std::min(arma::min(msa_correlation_2p), arma::min(mc_correlation_2p));
  correlation_2p_set = true;
};

void
MSACompare::computeCorrelation3p(void)
{
  if ((frequency_1p_set == false) | (frequency_2p_set == false) |
      (frequency_3p_set == false)) {
    std::cerr << "ERROR: compute the 1p, 2p, and 3p frequencies first"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  msa_correlation_3p = arma::Col<double>(
    (int)N * (N - 1) * (N - 2) / 6 * Q * Q * Q, arma::fill::zeros);
  mc_correlation_3p = arma::Col<double>(
    (int)N * (N - 1) * (N - 2) / 6 * Q * Q * Q, arma::fill::zeros);

#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int k = j + 1; k < N; k++) {
          for (int aa1 = 0; aa1 < Q; aa1++) {
            for (int aa2 = 0; aa2 < Q; aa2++) {
              for (int aa3 = 0; aa3 < Q; aa3++) {
                int idx = Q * Q * Q *
                            (N * (N - 1) * (N - 2) / 6 -
                             (N - i - 1) * (N - i) * (N - i - 2) / 6 + k - j +
                             (j - i - 1) * (N - i - 1) -
                             (j - i) * (j - i - 1) / 2 - 1) +
                          aa3 + aa2 * Q + aa1 * Q * Q;
                int idx2_12 =
                  aa2 + aa1 * Q +
                  Q * Q * (j - i - 1 + i * (N - 1) - i * (i - 1) / 2);
                int idx2_13 =
                  aa3 + aa1 * Q +
                  Q * Q * (k - i - 1 + i * (N - 1) - i * (i - 1) / 2);
                int idx2_23 =
                  aa3 + aa2 * Q +
                  Q * Q * (k - j - 1 + j * (N - 1) - j * (j - 1) / 2);

                int idx1_1 = aa1 + i * Q;
                int idx1_2 = aa2 + j * Q;
                int idx1_3 = aa3 + k * Q;

                msa_correlation_3p(idx) =
                  msa_frequency_3p(idx) -
                  msa_frequency_2p(idx2_12) * msa_frequency_1p(idx1_3) -
                  msa_frequency_2p(idx2_13) * msa_frequency_1p(idx1_2) -
                  msa_frequency_2p(idx2_23) * msa_frequency_1p(idx1_1) +
                  2.0 * msa_frequency_1p(idx1_1) * msa_frequency_1p(idx1_2) *
                    msa_frequency_1p(idx1_3);

                mc_correlation_3p(idx) =
                  mc_frequency_3p(idx) -
                  mc_frequency_2p(idx2_12) * mc_frequency_1p(idx1_3) -
                  mc_frequency_2p(idx2_13) * mc_frequency_1p(idx1_2) -
                  mc_frequency_2p(idx2_23) * mc_frequency_1p(idx1_1) +
                  2.0 * mc_frequency_1p(idx1_1) * mc_frequency_1p(idx1_2) *
                    mc_frequency_1p(idx1_3);
              }
            }
          }
        }
      }
    }
  }
  correlation_3p_set = true;

  corr_3p_max =
    std::max(arma::max(msa_correlation_3p), arma::max(mc_correlation_3p));
  corr_3p_min =
    std::min(arma::min(msa_correlation_3p), arma::min(mc_correlation_3p));
  corr_3p_max =
    std::max(arma::max(msa_correlation_3p), arma::max(mc_correlation_3p));
  corr_3p_min =
    std::min(arma::min(msa_correlation_3p), arma::min(mc_correlation_3p));
};

void
MSACompare::makeFrequency1pHistogram(void)
{
  if (frequency_1p_set == false) {
    std::cerr << "ERROR: compute the 1p frequencies first" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::wall_clock timer;
  std::cout << "computing 1p frequency regression... " << std::flush;
  timer.tic();

  double f = 1. / (double)Q;
  double xy = arma::accu((msa_frequency_1p - f) % (mc_frequency_1p - f));
  double xx = arma::accu(arma::pow((msa_frequency_1p - f), 2));
  double yy = arma::accu(arma::pow((mc_frequency_1p - f), 2));

  double b = xy / xx;
  double a = f * (1 - b);
  double r2 = (2 * b * xy - b * b * xx) / yy;
  std::cout << timer.toc() << " sec" << std::endl;
  std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
            << " (r=" << sqrt(r2) << ")" << std::endl;

  linear_model model;
  model.a = a;
  model.b = b;
  model.R2 = r2;

  histogram hist;
  hist.grid = arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
  double bin_width = 1. / (double)(BINS - 1);
  hist.bin_width = bin_width;
  for (int i = 0; i < Q * N; i++) {
    hist.grid(floor(msa_frequency_1p(i) / bin_width),
              floor(mc_frequency_1p(i) / bin_width))++;
  }
  writeHistogram("freq_1p_hist.tsv", hist);
  writeLinearModel("freq_1p_model.tsv", model);
  std::cout << std::endl;
};

void
MSACompare::makeFrequency2pHistogram(void)
{
  if ((frequency_1p_set == false) | (frequency_2p_set == false)) {
    std::cerr << "ERROR: compute the 1p and 2p frequencies first" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::wall_clock timer;
  std::cout << "computing 2p frequency regression... " << std::flush;
  timer.tic();

  double f = 1. / (double)(Q * Q);
  double xy = arma::accu((msa_frequency_2p - f) % (mc_frequency_2p - f));
  double xx = arma::accu(arma::pow((msa_frequency_2p - f), 2));
  double yy = arma::accu(arma::pow((mc_frequency_2p - f), 2));

  double b = xy / xx;
  double a = f * (1 - b);
  double r2 = (2 * b * xy - b * b * xx) / yy;
  std::cout << timer.toc() << " sec" << std::endl;
  std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
            << " (r=" << sqrt(r2) << ")" << std::endl;

  linear_model model;
  model.a = a;
  model.b = b;
  model.R2 = r2;

  histogram hist;
  hist.grid = arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
  double bin_width = 1. / (double)(BINS - 1);
  hist.bin_width = bin_width;
  for (int i = 0; i < N * (N - 1) / 2 * Q * Q; i++) {
    hist.grid(floor(msa_frequency_2p(i) / bin_width),
              floor(mc_frequency_2p(i) / bin_width))++;
  }
  // hist.grid.print("freq 2p:");
  writeHistogram("freq_2p_hist.tsv", hist);
  writeLinearModel("freq_2p_model.tsv", model);
  std::cout << std::endl;
};

void
MSACompare::makeCorrelation2pHistogram(void)
{
  if ((frequency_1p_set == false) | (frequency_2p_set == false)) {
    std::cerr << "ERROR: compute the 1p and 2p frequencies first" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  linear_model model;
  histogram hist;

  arma::wall_clock timer;
  std::cout << "computing 2p correlation regression... " << std::flush;
  timer.tic();
  if (correlation_2p_set) {
    double xy = arma::accu((msa_correlation_2p) % (mc_correlation_2p));
    double xx = arma::accu(arma::pow((msa_correlation_2p), 2));
    double yy = arma::accu(arma::pow((mc_correlation_2p), 2));

    double b = xy / xx;
    double a = 0;
    double r2 = (2 * b * xy - b * b * xx) / yy;
    std::cout << timer.toc() << " sec" << std::endl;
    std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;

    model.a = a;
    model.b = b;
    model.R2 = r2;

    hist.grid =
      arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
    double bin_width = 1.05 * (corr_2p_max - corr_2p_min) / (double)(BINS - 1);
    hist.bin_width = bin_width;
    hist.max = 1.05 * corr_2p_max;
    hist.min = 1.05 * corr_2p_min;

    for (int i = 0; i < N * (N - 1) / 2 * Q * Q; i++) {
      hist.grid(floor((msa_correlation_2p(i) - hist.min) / bin_width + .5),
                floor((mc_correlation_2p(i) - hist.min) / bin_width + .5))++;
    }
  } else {
    arma::Col<double> tmp_msa_correlation_2p =
      arma::Col<double>((int)N * (N - 1) / 2 * Q * Q, arma::fill::zeros);
    arma::Col<double> tmp_mc_correlation_2p =
      arma::Col<double>((int)N * (N - 1) / 2 * Q * Q, arma::fill::zeros);
#pragma omp parallel
    {
#pragma omp for schedule(dynamic,1)
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          for (int aa1 = 0; aa1 < Q; aa1++) {
            for (int aa2 = 0; aa2 < Q; aa2++) {
              int idx = aa2 + aa1 * Q +
                        Q * Q * (j - i - 1 + i * (N - 1) - i * (i - 1) / 2);
              tmp_msa_correlation_2p(idx) =
                msa_frequency_2p(idx) -
                msa_frequency_1p(i * Q + aa1) * msa_frequency_1p(j * Q + aa2);
              tmp_mc_correlation_2p(idx) =
                mc_frequency_2p(idx) -
                mc_frequency_1p(i * Q + aa1) * mc_frequency_1p(j * Q + aa2);
            }
          }
        }
      }
    }

    double xy = arma::accu((tmp_msa_correlation_2p) % (tmp_mc_correlation_2p));
    double xx = arma::accu(arma::pow((tmp_msa_correlation_2p), 2));
    double yy = arma::accu(arma::pow((tmp_mc_correlation_2p), 2));

    double b = xy / xx;
    double a = 0;
    double r2 = (2 * b * xy - b * b * xx) / yy;
    std::cout << timer.toc() << " sec" << std::endl;
    std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;
    std::cout << std::endl;

    model.a = a;
    model.b = b;
    model.R2 = r2;

    corr_2p_max = std::max(arma::max(tmp_msa_correlation_2p),
                           arma::max(tmp_mc_correlation_2p));
    corr_2p_min = std::min(arma::min(tmp_msa_correlation_2p),
                           arma::min(tmp_mc_correlation_2p));

    hist.grid =
      arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
    double bin_width = 1.05 * (corr_2p_max - corr_2p_min) / (double)(BINS - 1);
    hist.bin_width = bin_width;
    hist.max = 1.05 * corr_2p_max;
    hist.min = 1.05 * corr_2p_min;

    for (int i = 0; i < N * (N - 1) / 2 * Q * Q; i++) {
      hist.grid(
        floor((tmp_msa_correlation_2p(i) - hist.min) / bin_width + .5),
        floor((tmp_mc_correlation_2p(i) - hist.min) / bin_width + .5))++;
    }
  }
  // hist.grid.print("corr 2p:");
  writeHistogram("corr_2p_hist.tsv", hist);
  writeLinearModel("corr_2p_model.tsv", model);
};

void
MSACompare::makeEfficient3pHistograms(void)
{
  if ((frequency_1p_set == false) | (frequency_2p_set == false)) {
    std::cerr << "ERROR: compute the 1p and 2p frequencies first" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::wall_clock timer;
  std::cout << "computing 3p frequency and correlation regressions... "
            << std::flush;
  timer.tic();

  arma::Col<double> freq_xx_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> freq_xy_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> freq_yy_vec = arma::Col<double>(N, arma::fill::zeros);

  arma::Col<double> corr_xx_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> corr_xy_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> corr_yy_vec = arma::Col<double>(N, arma::fill::zeros);

  arma::Col<double> corr_max_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> corr_min_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> freq_max_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> freq_min_vec = arma::Col<double>(N, arma::fill::zeros);

  double f = 1. / (double)(Q * Q * Q);

  double freq_bin_width = 1. / (double)(BINS - 1);

  arma::Cube<unsigned long long int> freq_hist_i =
    arma::Cube<unsigned long long int>(BINS, BINS, N, arma::fill::zeros);
  arma::Cube<unsigned long long int> corr_hist_i =
    arma::Cube<unsigned long long int>(BINS, BINS, N, arma::fill::zeros);

  histogram freq_hist;
  freq_hist.grid =
    arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
  freq_hist.bin_width = freq_bin_width;
  freq_hist.bins = BINS;

  // double max = 1;
  // double min = -1;
  double max = .75 * corr_2p_max;
  double min = .75 * corr_2p_min;
  double corr_bin_width = (max - min) / (double)(BINS - 1);

  histogram corr_hist;
  corr_hist.grid =
    arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
  corr_hist.bin_width = corr_bin_width;
  corr_hist.bins = BINS;
  corr_hist.max = max;
  corr_hist.min = min;

#pragma omp parallel
  {
    double* msa_weight_ptr = msa->sequence_weights.memptr();
    double* mc_weight_ptr = mc->sequence_weights.memptr();
#pragma omp for schedule(dynamic,1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int k = j + 1; k < N; k++) {

          // double* msa_weight_ptr = msa->sequence_weights.memptr();
          // double* mc_weight_ptr = mc->sequence_weights.memptr();

          // int idx =
          //   (N * (N - 1) * (N - 2) / 6 -
          //    (N - i - 1) * (N - i) * (N - i - 2) / 6 + k - j +
          //    (j - i - 1) * (N - i - 1) - (j - i) * (j - i - 1) / 2 - 1);

          arma::Col<double> msa_frequency_3p_ijk =
            arma::Col<double>(Q * Q * Q, arma::fill::zeros);
          arma::Col<double> msa_correlation_3p_ijk =
            arma::Col<double>(Q * Q * Q, arma::fill::zeros);
          int* msa_align_ptr1 = msa->alignment.colptr(i);
          int* msa_align_ptr2 = msa->alignment.colptr(j);
          int* msa_align_ptr3 = msa->alignment.colptr(k);

          arma::Col<double> mc_frequency_3p_ijk =
            arma::Col<double>(Q * Q * Q, arma::fill::zeros);
          arma::Col<double> mc_correlation_3p_ijk =
            arma::Col<double>(Q * Q * Q, arma::fill::zeros);
          int* mc_align_ptr1 = mc->alignment.colptr(i);
          int* mc_align_ptr2 = mc->alignment.colptr(j);
          int* mc_align_ptr3 = mc->alignment.colptr(k);

          for (int m = 0; m < M_msa; m++) {
            int msa_idx = *(msa_align_ptr3 + m) + *(msa_align_ptr2 + m) * Q +
                          *(msa_align_ptr1 + m) * Q * Q;

            msa_frequency_3p_ijk(msa_idx) += *(msa_weight_ptr + m);
          }
          msa_frequency_3p_ijk = msa_frequency_3p_ijk / M_effective_msa;

          for (int m = 0; m < M_mc; m++) {
            int mc_idx = *(mc_align_ptr3 + m) + *(mc_align_ptr2 + m) * Q +
                         *(mc_align_ptr1 + m) * Q * Q;

            mc_frequency_3p_ijk(mc_idx) += *(mc_weight_ptr + m);
          }
          mc_frequency_3p_ijk = mc_frequency_3p_ijk / M_effective_mc;

          freq_xx_vec(i) += arma::accu(arma::pow(msa_frequency_3p_ijk - f, 2));
          freq_yy_vec(i) += arma::accu(arma::pow(mc_frequency_3p_ijk - f, 2));
          freq_xy_vec(i) +=
            arma::accu((msa_frequency_3p_ijk - f) % (mc_frequency_3p_ijk - f));

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

                int idx1_1 = aa1 + i * Q;
                int idx1_2 = aa2 + j * Q;
                int idx1_3 = aa3 + k * Q;

                msa_correlation_3p_ijk(idx_aa) =
                  msa_frequency_3p_ijk(idx_aa) -
                  msa_frequency_2p(idx2_12) * msa_frequency_1p(idx1_3) -
                  msa_frequency_2p(idx2_13) * msa_frequency_1p(idx1_2) -
                  msa_frequency_2p(idx2_23) * msa_frequency_1p(idx1_1) +
                  2.0 * msa_frequency_1p(idx1_1) * msa_frequency_1p(idx1_2) *
                    msa_frequency_1p(idx1_3);

                mc_correlation_3p_ijk(idx_aa) =
                  mc_frequency_3p_ijk(idx_aa) -
                  mc_frequency_2p(idx2_12) * mc_frequency_1p(idx1_3) -
                  mc_frequency_2p(idx2_13) * mc_frequency_1p(idx1_2) -
                  mc_frequency_2p(idx2_23) * mc_frequency_1p(idx1_1) +
                  2.0 * mc_frequency_1p(idx1_1) * mc_frequency_1p(idx1_2) *
                    mc_frequency_1p(idx1_3);
              }
            }
          }

          corr_xx_vec(i) += arma::accu(arma::pow(msa_correlation_3p_ijk, 2));
          corr_yy_vec(i) += arma::accu(arma::pow(mc_correlation_3p_ijk, 2));
          corr_xy_vec(i) +=
            arma::accu((msa_correlation_3p_ijk) % (mc_correlation_3p_ijk));

          {
            arma::Col<double> tmp_msa_frequency_3p_ijk =
              arma::floor(msa_frequency_3p_ijk / freq_bin_width);
            arma::Col<double> tmp_mc_frequency_3p_ijk =
              arma::floor(mc_frequency_3p_ijk / freq_bin_width);

            arma::Col<double> tmp_msa_correlation_3p_ijk =
              arma::floor((msa_correlation_3p_ijk - min) / corr_bin_width + .5);
            arma::Col<double> tmp_mc_correlation_3p_ijk =
              arma::floor((mc_correlation_3p_ijk - min) / corr_bin_width + .5);

            for (int q = 0; q < Q * Q * Q; q++) {
              freq_hist_i(
                tmp_msa_frequency_3p_ijk(q), tmp_mc_frequency_3p_ijk(q), i)++;
              corr_hist_i(tmp_msa_correlation_3p_ijk(q),
                          tmp_mc_correlation_3p_ijk(q),
                          i)++;
            }
          }

          {
            double tmp_max = std::max(arma::max(msa_frequency_3p_ijk),
                                      arma::max(mc_frequency_3p_ijk));
            double tmp_min = std::min(arma::min(msa_frequency_3p_ijk),
                                      arma::min(mc_frequency_3p_ijk));
            if (tmp_max > freq_max_vec(i)) {
              freq_max_vec(i) = tmp_max;
            }
            if (tmp_min < freq_min_vec(i)) {
              freq_min_vec(i) = tmp_min;
            }
          }

          {
            double tmp_max = std::max(arma::max(msa_correlation_3p_ijk),
                                      arma::max(mc_correlation_3p_ijk));
            double tmp_min = std::min(arma::min(msa_correlation_3p_ijk),
                                      arma::min(mc_correlation_3p_ijk));
            if (tmp_max > corr_max_vec(i)) {
              corr_max_vec(i) = tmp_max;
            }
            if (tmp_min < corr_min_vec(i)) {
              corr_min_vec(i) = tmp_min;
            }
          }
        }
      }
    }
  }
  std::cout << timer.toc() << " sec" << std::endl;

  freq_3p_min = arma::min(freq_min_vec);
  freq_3p_max = arma::max(freq_max_vec);
  corr_3p_min = arma::min(corr_min_vec);
  corr_3p_max = arma::max(corr_max_vec);

  {
    double xx = arma::accu(freq_xx_vec);
    double xy = arma::accu(freq_xy_vec);
    double yy = arma::accu(freq_yy_vec);

    double b = xy / xx;
    double a = f * (1 - b);
    double r2 = (2 * b * xy - b * b * xx) / yy;
    std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;
    // std::cout << std::endl;

    linear_model model;
    model.a = a;
    model.b = b;
    model.R2 = r2;

    freq_hist.grid = arma::sum(freq_hist_i, 2);
    // freq_hist.grid.print("freq 3p:");
    writeHistogram("freq_3p_hist.tsv", freq_hist);
    writeLinearModel("freq_3p_model.tsv", model);
  }
  {
    double xx = arma::accu(corr_xx_vec);
    double xy = arma::accu(corr_xy_vec);
    double yy = arma::accu(corr_yy_vec);

    double b = xy / xx;
    double a = 0;
    double r2 = (2 * b * xy - b * b * xx) / yy;
    std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;
    std::cout << std::endl;

    linear_model model;
    model.a = a;
    model.b = b;
    model.R2 = r2;

    corr_hist.grid = arma::sum(corr_hist_i, 2);
    // corr_hist.grid.print("corr 3p:");
    writeHistogram("corr_3p_hist.tsv", corr_hist);
    writeLinearModel("corr_3p_model.tsv", model);
  }

  // std::cout << "efficient freq 3p: " << freq_3p_min << ", " << freq_3p_max <<
  //   std::endl;
  // std::cout << "efficient corr 3p: " << corr_3p_min << ", " << corr_3p_max <<
  //   std::endl;
};

void
MSACompare::makeFrequency3pHistogram(void)
{
  if ((frequency_1p_set == false) | (frequency_2p_set == false) |
      (frequency_3p_set == false)) {
    std::cerr << "ERROR: compute the 1p, 2p, and 3p frequencies first"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::wall_clock timer;
  std::cout << "computing 3p frequency regression... " << std::flush;
  timer.tic();

  double f = 1. / (double)(Q * Q * Q);
  double xy = arma::accu((msa_frequency_3p - f) % (mc_frequency_3p - f));
  double xx = arma::accu(arma::pow((msa_frequency_3p - f), 2));
  double yy = arma::accu(arma::pow((mc_frequency_3p - f), 2));

  double b = xy / xx;
  double a = f * (1 - b);
  double r2 = (2 * b * xy - b * b * xx) / yy;
  std::cout << timer.toc() << " sec" << std::endl;
  std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
            << " (r=" << sqrt(r2) << ")" << std::endl;
  std::cout << std::endl;

  linear_model model;
  model.a = a;
  model.b = b;
  model.R2 = r2;

  histogram hist;
  hist.grid = arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
  double bin_width = 1. / (double)(BINS - 1);
  hist.bin_width = bin_width;
  for (int i = 0; i < N * (N - 1) * (N - 2) / 6 * Q * Q * Q; i++) {
    hist.grid(floor(msa_frequency_3p(i) / bin_width),
              floor(mc_frequency_3p(i) / bin_width))++;
  }
  // hist.grid.print("freq 3p:");
  writeHistogram("freq_3p_hist.tsv", hist);
  writeLinearModel("freq_3p_model.tsv", model);
};

void
MSACompare::makeCorrelation3pHistogram(void)
{
  if ((frequency_1p_set == false) | (frequency_2p_set == false) |
      (frequency_3p_set == false)) {
    std::cerr << "ERROR: compute the 1p, 2p, and 3p frequencies first"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "3p corr range (start): " << corr_3p_min << ", " << corr_3p_max
            << std::endl;

  linear_model model;
  histogram hist;

  arma::wall_clock timer;
  std::cout << "computing 3p correlation regression... " << std::flush;
  timer.tic();

  if (correlation_3p_set) {
    double xy = arma::accu((msa_correlation_3p) % (mc_correlation_3p));
    double xx = arma::accu(arma::pow((msa_correlation_3p), 2));
    double yy = arma::accu(arma::pow((mc_correlation_3p), 2));

    double b = xy / xx;
    double a = 0;
    double r2 = (2 * b * xy - b * b * xx) / yy;
    std::cout << timer.toc() << " sec" << std::endl;
    std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;
    std::cout << std::endl;

    model.a = a;
    model.b = b;
    model.R2 = r2;

    hist.grid =
      arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
    double bin_width = 1.05 * (corr_3p_max - corr_3p_min) / (double)(BINS - 1);
    hist.bin_width = bin_width;
    hist.max = 1.05 * corr_3p_max;
    hist.min = 1.05 * corr_3p_min;

    for (int i = 0; i < N * (N - 1) * (N - 2) / 6 * Q * Q * Q; i++) {
      hist.grid(floor((msa_correlation_3p(i) - hist.min) / bin_width + .5),
                floor((mc_correlation_3p(i) - hist.min) / bin_width + .5))++;
    }
  } else {
    arma::Col<double> tmp_msa_correlation_3p = arma::Col<double>(
      (int)N * (N - 1) * (N - 2) / 6 * Q * Q * Q, arma::fill::zeros);
    arma::Col<double> tmp_mc_correlation_3p = arma::Col<double>(
      (int)N * (N - 1) * (N - 2) / 6 * Q * Q * Q, arma::fill::zeros);

    arma::Col<double> corr_max_vec = arma::Col<double>(N, arma::fill::zeros);
    arma::Col<double> corr_min_vec = arma::Col<double>(N, arma::fill::zeros);

#pragma omp parallel
    {
#pragma omp for schedule(dynamic,1)
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          for (int k = j + 1; k < N; k++) {
            for (int aa1 = 0; aa1 < Q; aa1++) {
              for (int aa2 = 0; aa2 < Q; aa2++) {
                for (int aa3 = 0; aa3 < Q; aa3++) {
                  int idx = Q * Q * Q *
                              (N * (N - 1) * (N - 2) / 6 -
                               (N - i - 1) * (N - i) * (N - i - 2) / 6 + k - j +
                               (j - i - 1) * (N - i - 1) -
                               (j - i) * (j - i - 1) / 2 - 1) +
                            aa3 + aa2 * Q + aa1 * Q * Q;
                  int idx2_12 =
                    aa2 + aa1 * Q +
                    Q * Q * (j - i - 1 + i * (N - 1) - i * (i - 1) / 2);
                  int idx2_13 =
                    aa3 + aa1 * Q +
                    Q * Q * (k - i - 1 + i * (N - 1) - i * (i - 1) / 2);
                  int idx2_23 =
                    aa3 + aa2 * Q +
                    Q * Q * (k - j - 1 + j * (N - 1) - j * (j - 1) / 2);

                  tmp_msa_correlation_3p(idx) =
                    msa_frequency_3p(idx) -
                    msa_frequency_2p(idx2_12) * msa_frequency_1p(aa3 + k * Q) -
                    msa_frequency_2p(idx2_13) * msa_frequency_1p(aa2 + j * Q) -
                    msa_frequency_2p(idx2_23) * msa_frequency_1p(aa1 + i * Q) +
                    2.0 * msa_frequency_1p(aa1 + i * Q) *
                      msa_frequency_1p(aa2 + j * Q) *
                      msa_frequency_1p(aa3 + k * Q);

                  tmp_mc_correlation_3p(idx) =
                    mc_frequency_3p(idx) -
                    mc_frequency_2p(idx2_12) * mc_frequency_1p(aa3 + k * Q) -
                    mc_frequency_2p(idx2_13) * mc_frequency_1p(aa2 + j * Q) -
                    mc_frequency_2p(idx2_23) * mc_frequency_1p(aa1 + i * Q) +
                    2.0 * mc_frequency_1p(aa1 + i * Q) *
                      mc_frequency_1p(aa2 + j * Q) *
                      mc_frequency_1p(aa3 + k * Q);
                }
              }
            }
          }
          // corr_3p_max = std::max(arma::max(tmp_msa_correlation_3p),
          //                        arma::max(tmp_mc_correlation_3p));
          // corr_3p_min = std::min(arma::min(tmp_msa_correlation_3p),
          //                        arma::min(tmp_mc_correlation_3p));
        }
      }
    }

    double xy = arma::accu((tmp_msa_correlation_3p) % (tmp_mc_correlation_3p));
    double xx = arma::accu(arma::pow((tmp_msa_correlation_3p), 2));
    double yy = arma::accu(arma::pow((tmp_mc_correlation_3p), 2));

    double b = xy / xx;
    double a = 0;
    double r2 = (2 * b * xy - b * b * xx) / yy;
    std::cout << timer.toc() << " sec" << std::endl;
    std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;
    std::cout << std::endl;

    model.a = a;
    model.b = b;
    model.R2 = r2;

    std::cout << "basic corr 3p: " << corr_3p_min << ", " << corr_3p_max
              << std::endl;

    corr_3p_max = std::max(arma::max(tmp_msa_correlation_3p),
                           arma::max(tmp_mc_correlation_3p));
    corr_3p_min = std::min(arma::min(tmp_msa_correlation_3p),
                           arma::min(tmp_mc_correlation_3p));

    std::cout << "basic corr 3p (after): " << corr_3p_min << ", " << corr_3p_max
              << std::endl;

    hist.grid =
      arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
    double bin_width = 1.05 * (corr_3p_max - corr_3p_min) / (double)(BINS - 1);
    hist.bin_width = bin_width;
    hist.max = 1.05 * corr_3p_max;
    hist.min = 1.05 * corr_3p_min;

    for (int i = 0; i < N * (N - 1) * (N - 2) / 6 * Q * Q * Q; i++) {
      hist.grid(
        floor((tmp_msa_correlation_3p(i) - hist.min) / bin_width + .5),
        floor((tmp_mc_correlation_3p(i) - hist.min) / bin_width + .5))++;
    }
  }
  // hist.grid.print("corr 3p:");
  writeHistogram("corr_3p_hist.tsv", hist);
  writeLinearModel("corr_3p_model.tsv", model);
};

void
MSACompare::makeEfficient4pHistograms(void)
{
  if ((frequency_1p_set == false) | (frequency_2p_set == false) |
      (frequency_3p_set == false)) {
    std::cerr << "ERROR: compute the 1p, 2p, and 3p frequencies first"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::wall_clock timer;
  std::cout << "computing 4p frequency and correlation regressions... "
            << std::flush;
  timer.tic();

  // std::cout << "flag 0" << std::endl;
  arma::Col<double> freq_xx_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> freq_xy_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> freq_yy_vec = arma::Col<double>(N, arma::fill::zeros);

  arma::Col<double> corr_xx_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> corr_xy_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> corr_yy_vec = arma::Col<double>(N, arma::fill::zeros);

  arma::Col<double> corr_max_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> corr_min_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> freq_max_vec = arma::Col<double>(N, arma::fill::zeros);
  arma::Col<double> freq_min_vec = arma::Col<double>(N, arma::fill::zeros);

  double f = 1. / (double)(Q * Q * Q * Q);

  double freq_bin_width = 1. / (double)(BINS - 1);

  arma::Cube<unsigned long long int> freq_hist_i =
    arma::Cube<unsigned long long int>(BINS, BINS, N, arma::fill::zeros);
  arma::Cube<unsigned long long int> corr_hist_i =
    arma::Cube<unsigned long long int>(BINS, BINS, N, arma::fill::zeros);

  histogram freq_hist;
  freq_hist.grid =
    arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
  freq_hist.bin_width = freq_bin_width;
  freq_hist.bins = BINS;

  // double max = 1;
  // double min = -1;
  double max = .75 * corr_3p_max;
  double min = .75 * corr_3p_min;
  double corr_bin_width = (max - min) / (double)(BINS - 1);

  histogram corr_hist;
  corr_hist.grid =
    arma::Mat<unsigned long long int>(BINS, BINS, arma::fill::zeros);
  corr_hist.bin_width = corr_bin_width;
  corr_hist.bins = BINS;
  corr_hist.max = max;
  corr_hist.min = min;

  // std::cout << "3p corr range: " << corr_3p_min << ", " << corr_3p_max <<
  //   std::endl;

#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1)
    for (int i = 0; i < N; i++) {
      double* msa_weight_ptr = msa->sequence_weights.memptr();
      double* mc_weight_ptr = mc->sequence_weights.memptr();

      for (int j = i + 1; j < N; j++) {
        for (int k = j + 1; k < N; k++) {
          for (int l = k + 1; l < N; l++) {
            arma::Col<double> msa_frequency_4p_ijkl =
              arma::Col<double>(Q * Q * Q * Q, arma::fill::zeros);
            arma::Col<double> msa_correlation_4p_ijkl =
              arma::Col<double>(Q * Q * Q * Q, arma::fill::zeros);
            int* msa_align_ptr1 = msa->alignment.colptr(i);
            int* msa_align_ptr2 = msa->alignment.colptr(j);
            int* msa_align_ptr3 = msa->alignment.colptr(k);
            int* msa_align_ptr4 = msa->alignment.colptr(l);

            arma::Col<double> mc_frequency_4p_ijkl =
              arma::Col<double>(Q * Q * Q * Q, arma::fill::zeros);
            arma::Col<double> mc_correlation_4p_ijkl =
              arma::Col<double>(Q * Q * Q * Q, arma::fill::zeros);
            int* mc_align_ptr1 = mc->alignment.colptr(i);
            int* mc_align_ptr2 = mc->alignment.colptr(j);
            int* mc_align_ptr3 = mc->alignment.colptr(k);
            int* mc_align_ptr4 = mc->alignment.colptr(l);

            for (int m = 0; m < M_msa; m++) {
              int msa_idx = *(msa_align_ptr4 + m) + *(msa_align_ptr3 + m) * Q +
                            *(msa_align_ptr2 + m) * Q * Q +
                            *(msa_align_ptr1 + m) * Q * Q * Q;

              msa_frequency_4p_ijkl(msa_idx) += *(msa_weight_ptr + m);
            }
            msa_frequency_4p_ijkl = msa_frequency_4p_ijkl / M_effective_msa;

            for (int m = 0; m < M_mc; m++) {
              int mc_idx = *(mc_align_ptr4 + m) + *(mc_align_ptr3 + m) * Q +
                           *(mc_align_ptr2 + m) * Q * Q +
                           *(mc_align_ptr1 + m) * Q * Q * Q;

              mc_frequency_4p_ijkl(mc_idx) += *(mc_weight_ptr + m);
            }
            mc_frequency_4p_ijkl = mc_frequency_4p_ijkl / M_effective_mc;

            freq_xx_vec(i) +=
              arma::accu(arma::pow(msa_frequency_4p_ijkl - f, 2));
            freq_yy_vec(i) +=
              arma::accu(arma::pow(mc_frequency_4p_ijkl - f, 2));
            freq_xy_vec(i) += arma::accu((msa_frequency_4p_ijkl - f) %
                                         (mc_frequency_4p_ijkl - f));

            for (int aa1 = 0; aa1 < Q; aa1++) {
              for (int aa2 = 0; aa2 < Q; aa2++) {
                for (int aa3 = 0; aa3 < Q; aa3++) {
                  for (int aa4 = 0; aa4 < Q; aa4++) {
                    int idx_aa = aa4 + aa3 * Q + aa2 * Q * Q + aa1 * Q * Q * Q;

                    int idx2_12 =
                      aa2 + aa1 * Q +
                      Q * Q * (j - i - 1 + i * (N - 1) - i * (i - 1) / 2);
                    int idx2_13 =
                      aa3 + aa1 * Q +
                      Q * Q * (k - i - 1 + i * (N - 1) - i * (i - 1) / 2);
                    int idx2_14 =
                      aa4 + aa1 * Q +
                      Q * Q * (l - i - 1 + i * (N - 1) - i * (i - 1) / 2);
                    int idx2_23 =
                      aa3 + aa2 * Q +
                      Q * Q * (k - j - 1 + j * (N - 1) - j * (j - 1) / 2);
                    int idx2_24 =
                      aa4 + aa2 * Q +
                      Q * Q * (l - j - 1 + j * (N - 1) - j * (j - 1) / 2);
                    int idx2_34 =
                      aa4 + aa3 * Q +
                      Q * Q * (l - k - 1 + k * (N - 1) - k * (k - 1) / 2);

                    int idx3_123 = aa1 * Q * Q + aa2 * Q + aa3 +
                                   Q * Q * Q *
                                     (N * (N - 1) * (N - 2) / 6 -
                                      (N - i - 1) * (N - i) * (N - i - 2) / 6 +
                                      k - j + (j - i - 1) * (N - i - 1) -
                                      (j - i) * (j - i - 1) / 2 - 1);
                    int idx3_124 = aa1 * Q * Q + aa2 * Q + aa4 +
                                   Q * Q * Q *
                                     (N * (N - 1) * (N - 2) / 6 -
                                      (N - i - 1) * (N - i) * (N - i - 2) / 6 +
                                      l - j + (j - i - 1) * (N - i - 1) -
                                      (j - i) * (j - i - 1) / 2 - 1);
                    int idx3_134 = aa1 * Q * Q + aa3 * Q + aa4 +
                                   Q * Q * Q *
                                     (N * (N - 1) * (N - 2) / 6 -
                                      (N - i - 1) * (N - i) * (N - i - 2) / 6 +
                                      l - k + (k - i - 1) * (N - i - 1) -
                                      (k - i) * (k - i - 1) / 2 - 1);
                    int idx3_234 = aa2 * Q * Q + aa3 * Q + aa4 +
                                   Q * Q * Q *
                                     (N * (N - 1) * (N - 2) / 6 -
                                      (N - j - 1) * (N - j) * (N - j - 2) / 6 +
                                      l - k + (k - j - 1) * (N - j - 1) -
                                      (k - j) * (k - j - 1) / 2 - 1);

                    int idx1_1 = aa1 + i * Q;
                    int idx1_2 = aa2 + j * Q;
                    int idx1_3 = aa3 + k * Q;
                    int idx1_4 = aa4 + l * Q;

                    msa_correlation_4p_ijkl(idx_aa) =
                      msa_frequency_4p_ijkl(idx_aa) -
                      msa_frequency_1p(idx1_1) * msa_frequency_3p(idx3_234) -
                      msa_frequency_1p(idx1_2) * msa_frequency_3p(idx3_134) -
                      msa_frequency_1p(idx1_3) * msa_frequency_3p(idx3_124) -
                      msa_frequency_1p(idx1_4) * msa_frequency_3p(idx3_123) -
                      msa_frequency_2p(idx2_12) * msa_frequency_2p(idx2_34) -
                      msa_frequency_2p(idx2_13) * msa_frequency_2p(idx2_24) -
                      msa_frequency_2p(idx2_14) * msa_frequency_2p(idx2_23) +
                      2 * msa_frequency_2p(idx2_12) * msa_frequency_1p(idx1_3) *
                        msa_frequency_1p(idx1_4) +
                      2 * msa_frequency_2p(idx2_13) * msa_frequency_1p(idx1_2) *
                        msa_frequency_1p(idx1_4) +
                      2 * msa_frequency_2p(idx2_14) * msa_frequency_1p(idx1_2) *
                        msa_frequency_1p(idx1_3) +
                      2 * msa_frequency_2p(idx2_23) * msa_frequency_1p(idx1_1) *
                        msa_frequency_1p(idx1_4) +
                      2 * msa_frequency_2p(idx2_24) * msa_frequency_1p(idx1_1) *
                        msa_frequency_1p(idx1_3) +
                      2 * msa_frequency_2p(idx2_34) * msa_frequency_1p(idx1_1) *
                        msa_frequency_1p(idx1_2) -
                      6 * msa_frequency_1p(idx1_1) * msa_frequency_1p(idx1_2) *
                        msa_frequency_1p(idx1_3) * msa_frequency_1p(idx1_4);

                    mc_correlation_4p_ijkl(idx_aa) =
                      mc_frequency_4p_ijkl(idx_aa) -
                      mc_frequency_1p(idx1_1) * mc_frequency_3p(idx3_234) -
                      mc_frequency_1p(idx1_2) * mc_frequency_3p(idx3_134) -
                      mc_frequency_1p(idx1_3) * mc_frequency_3p(idx3_124) -
                      mc_frequency_1p(idx1_4) * mc_frequency_3p(idx3_123) -
                      mc_frequency_2p(idx2_12) * mc_frequency_2p(idx2_34) -
                      mc_frequency_2p(idx2_13) * mc_frequency_2p(idx2_24) -
                      mc_frequency_2p(idx2_14) * mc_frequency_2p(idx2_23) +
                      2 * mc_frequency_2p(idx2_12) * mc_frequency_1p(idx1_3) *
                        mc_frequency_1p(idx1_4) +
                      2 * mc_frequency_2p(idx2_13) * mc_frequency_1p(idx1_2) *
                        mc_frequency_1p(idx1_4) +
                      2 * mc_frequency_2p(idx2_14) * mc_frequency_1p(idx1_2) *
                        mc_frequency_1p(idx1_3) +
                      2 * mc_frequency_2p(idx2_23) * mc_frequency_1p(idx1_1) *
                        mc_frequency_1p(idx1_4) +
                      2 * mc_frequency_2p(idx2_24) * mc_frequency_1p(idx1_1) *
                        mc_frequency_1p(idx1_3) +
                      2 * mc_frequency_2p(idx2_34) * mc_frequency_1p(idx1_1) *
                        mc_frequency_1p(idx1_2) -
                      6 * mc_frequency_1p(idx1_1) * mc_frequency_1p(idx1_2) *
                        mc_frequency_1p(idx1_3) * mc_frequency_1p(idx1_4);
                  }
                }
              }
            }

            corr_xx_vec(i) += arma::accu(arma::pow(msa_correlation_4p_ijkl, 2));
            corr_yy_vec(i) += arma::accu(arma::pow(mc_correlation_4p_ijkl, 2));
            corr_xy_vec(i) +=
              arma::accu((msa_correlation_4p_ijkl) % (mc_correlation_4p_ijkl));

            {
              arma::Col<double> tmp_msa_frequency_4p_ijkl =
                arma::floor(msa_frequency_4p_ijkl / freq_bin_width);
              arma::Col<double> tmp_mc_frequency_4p_ijkl =
                arma::floor(mc_frequency_4p_ijkl / freq_bin_width);

              arma::Col<double> tmp_msa_correlation_4p_ijkl = arma::floor(
                (msa_correlation_4p_ijkl - min) / corr_bin_width + .5);
              arma::Col<double> tmp_mc_correlation_4p_ijkl = arma::floor(
                (mc_correlation_4p_ijkl - min) / corr_bin_width + .5);

              for (int q = 0; q < Q * Q * Q * Q; q++) {
                freq_hist_i(tmp_msa_frequency_4p_ijkl(q),
                            tmp_mc_frequency_4p_ijkl(q),
                            i)++;
                corr_hist_i(tmp_msa_correlation_4p_ijkl(q),
                            tmp_mc_correlation_4p_ijkl(q),
                            i)++;
              }
            }

            {
              double tmp_max = std::max(arma::max(msa_frequency_4p_ijkl),
                                        arma::max(mc_frequency_4p_ijkl));
              double tmp_min = std::min(arma::min(msa_frequency_4p_ijkl),
                                        arma::min(mc_frequency_4p_ijkl));
              if (tmp_max > freq_max_vec(i)) {
                freq_max_vec(i) = tmp_max;
              }
              if (tmp_min < freq_min_vec(i)) {
                freq_min_vec(i) = tmp_min;
              }
            }

            {
              double tmp_max = std::max(arma::max(msa_correlation_4p_ijkl),
                                        arma::max(mc_correlation_4p_ijkl));
              double tmp_min = std::min(arma::min(msa_correlation_4p_ijkl),
                                        arma::min(mc_correlation_4p_ijkl));
              if (tmp_max > corr_max_vec(i)) {
                corr_max_vec(i) = tmp_max;
              }
              if (tmp_min < corr_min_vec(i)) {
                corr_min_vec(i) = tmp_min;
              }
            }
          }
        }
      }
    }
  }
  std::cout << timer.toc() << " sec" << std::endl;

  freq_4p_min = arma::min(freq_min_vec);
  freq_4p_max = arma::max(freq_max_vec);
  corr_4p_min = arma::min(corr_min_vec);
  corr_4p_max = arma::max(corr_max_vec);

  {
    double xx = arma::accu(freq_xx_vec);
    double xy = arma::accu(freq_xy_vec);
    double yy = arma::accu(freq_yy_vec);

    double b = xy / xx;
    double a = f * (1 - b);
    double r2 = (2 * b * xy - b * b * xx) / yy;
    std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;
    // std::cout << std::endl;

    linear_model model;
    model.a = a;
    model.b = b;
    model.R2 = r2;

    freq_hist.grid = arma::sum(freq_hist_i, 2);
    // freq_hist.grid.print("freq 4p:");
    writeHistogram("freq_4p_hist.tsv", freq_hist);
    writeLinearModel("freq_4p_model.tsv", model);
  }
  {
    double xx = arma::accu(corr_xx_vec);
    double xy = arma::accu(corr_xy_vec);
    double yy = arma::accu(corr_yy_vec);

    double b = xy / xx;
    double a = 0;
    double r2 = (2 * b * xy - b * b * xx) / yy;
    std::cout << "model: y=" << b << "x+" << a << ", R^2=" << r2
              << " (r=" << sqrt(r2) << ")" << std::endl;
    std::cout << std::endl;

    linear_model model;
    model.a = a;
    model.b = b;
    model.R2 = r2;

    corr_hist.grid = arma::sum(corr_hist_i, 2);
    // corr_hist.grid.print("corr 4p:");
    writeHistogram("corr_4p_hist.tsv", corr_hist);
    writeLinearModel("corr_4p_model.tsv", model);
  }

  // std::cout << "4p freq range: " << freq_4p_min << ", " << freq_4p_max <<
  //   std::endl;
  // std::cout << "4p corr range: " << corr_4p_min << ", " << corr_4p_max <<
  //   std::endl;
};

void
MSACompare::writeHistogram(std::string file, histogram hist)
{
  std::ofstream output_stream(file);
  for (int i = 0; i < BINS; i++) {
    for (int j = 0; j < BINS; j++) {
      output_stream << hist.min + i * hist.bin_width << "\t"
                    << hist.min + j * hist.bin_width << "\t" << hist.grid(i, j)
                    << std::endl;
    }
  }
  output_stream.close();
  return;
};

void
MSACompare::writeLinearModel(std::string file, linear_model model)
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

void
MSACompare::writeFrequency1p(void)
{
  std::ofstream msa_output_stream("msa_freq_1p.txt");
  std::ofstream mc_output_stream("mc_freq_1p.txt");

  for (int i = 0; i < N; i++) {
    msa_output_stream << i;
    mc_output_stream << i;
    for (int aa = 0; aa < Q; aa++) {
      msa_output_stream << " " << msa_frequency_1p(aa + i * Q);
      mc_output_stream << " " << mc_frequency_1p(aa + i * Q);
    }
    msa_output_stream << std::endl;
    mc_output_stream << std::endl;
  }
  msa_output_stream.close();
  mc_output_stream.close();
};

void
MSACompare::writeFrequency2p(void)
{
  std::ofstream msa_output_stream("msa_freq_2p.txt");
  std::ofstream mc_output_stream("mc_freq_2p.txt");

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      msa_output_stream << i << " " << j;
      mc_output_stream << i << " " << j;
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          int idx =
            aa2 + aa1 * Q + Q * Q * (j - i - 1 + i * (N - 1) - i * (i - 1) / 2);
          msa_output_stream << " " << msa_frequency_2p(idx);
          mc_output_stream << " " << mc_frequency_2p(idx);
        }
      }
      msa_output_stream << std::endl;
      mc_output_stream << std::endl;
    }
  }
  msa_output_stream.close();
  mc_output_stream.close();
};
