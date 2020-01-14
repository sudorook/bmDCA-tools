#! /bin/bash
set -eu

./plot_stats.py \
  -s msa_freq_1p.txt \
  -c mcmc_new_freq_1p.txt \
  -t "1st order statistics (new)" -o "mcmc_new_1p"
./plot_stats.py \
  -s msa_corr_2p.txt \
  -c mcmc_new_corr_2p.txt \
  -t "2nd order statistics (new)" -o "mcmc_new_2p"
./plot_stats.py \
  -s msa_corr_3p.txt \
  -c mcmc_new_corr_3p.txt \
  -t "3rd order statistics (new)" -o "mcmc_new_3p"

./plot_stats.py \
  -s msa_freq_1p.txt \
  -c mcmc_old_freq_1p.txt \
  -t "1st order statistics (old)" -o "mcmc_old_1p"
./plot_stats.py \
  -s msa_corr_2p.txt \
  -c mcmc_old_corr_2p.txt \
  -t "2nd order statistics (old)" -o "mcmc_old_2p"
./plot_stats.py \
  -s msa_corr_3p.txt \
  -c mcmc_old_corr_3p.txt \
  -t "3rd order statistics (old)" -o "mcmc_old_3p"

./plot_stats.py \
  -s msa_freq_1p.txt \
  -c mcmc_reg_freq_1p.txt \
  -t "1st order statistics (reg)" -o "mcmc_reg_1p"
./plot_stats.py \
  -s msa_corr_2p.txt \
  -c mcmc_reg_corr_2p.txt \
  -t "2nd order statistics (reg)" -o "mcmc_reg_2p"
./plot_stats.py \
  -s msa_corr_3p.txt \
  -c mcmc_reg_corr_3p.txt \
  -t "3rd order statistics (reg)" -o "mcmc_reg_3p"
