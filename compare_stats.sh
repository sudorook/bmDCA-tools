#! /bin/bash
set -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

MSA="$1"
MCMC="$2"
DESCRIPTION="$3"
THRESHOLD="$4"

SLUG=$(echo "$DESCRIPTION" | tr '[:upper:]' '[:lower:]' | sed -e "s/ /_/g")

"${SCRIPT_DIR}/compare_stats" \
  -s "$MSA" \
  -c "$MCMC" \
  -t "$THRESHOLD"

"${SCRIPT_DIR}/plot_stats.py" \
  -s "${MSA%.*}_freq_1p.txt" \
  -c "${MCMC%.*}_freq_1p.txt" \
  -t "1p: MSA vs MCMC (${DESCRIPTION})" \
  -o "msa_mcmc_${SLUG}_1p"
"${SCRIPT_DIR}/plot_stats.py" \
  -s "${MSA%.*}_corr_2p.txt" \
  -c "${MCMC%.*}_corr_2p.txt" \
  -t "2p: MSA vs MCMC (${DESCRIPTION})" \
  -o "msa_mcmc_${SLUG}_2p"
"${SCRIPT_DIR}/plot_stats.py" \
  -s "${MSA%.*}_corr_3p.txt" \
  -c "${MCMC%.*}_corr_3p.txt" \
  -t "3p: MSA vs MCMC (${DESCRIPTION})" \
  -o "msa_mcmc_${SLUG}_3p"

rm -v \
  "${MSA%.*}_freq_1p.txt" \
  "${MCMC%.*}_freq_1p.txt" \
  "${MSA%.*}_corr_2p.txt" \
  "${MCMC%.*}_corr_2p.txt" \
  "${MSA%.*}_corr_3p.txt" \
  "${MCMC%.*}_corr_3p.txt"
