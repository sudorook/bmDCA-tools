#! /bin/bash
set -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

OPTIONS=s:c:d:t:
LONGOPTIONS=msa:,mcmc:,description:,threshold:
PARSED=$(getopt -o ${OPTIONS} --long ${LONGOPTIONS} -n "$0" -- "$@")
eval set -- "$PARSED"

while [ $# -ge 1 ]; do
  case "$1" in
    -s|--msa)
      MSA="$2"
      shift 2
      ;;
    -c|--mcmc)
      MCMC="$2"
      shift 2
      ;;
    -d|--description)
      DESCRIPTION="$2"
      shift 2
      ;;
    -t|--threshold)
      THRESHOLD="$2"
      shift 2
      ;;
    --)
      shift
      break
      ;;
    *)
      echo "DANGER, WILL ROBINSON!"
      exit 3
      ;;
  esac
done

SLUG=$(echo "$DESCRIPTION" | tr '[:upper:]' '[:lower:]' | sed -e "s/ /_/g" | \
       sed -e "s/,//g" -e "s/(//g" -e "s/)//g")

"${SCRIPT_DIR}/compare_stats" \
  -s "$MSA" \
  -c "$MCMC" \
  -t "$THRESHOLD"

echo "plotting 1p frequencies"
"${SCRIPT_DIR}/plot_stats.py" \
  -s "${MSA%.*}_freq_1p.txt" \
  -c "${MCMC%.*}_freq_1p.txt" \
  -t "1p: MSA vs MCMC (${DESCRIPTION})" \
  -o "msa_mcmc_${SLUG}_1p"

echo "plotting 2p correlations"
"${SCRIPT_DIR}/plot_stats.py" \
  -s "${MSA%.*}_corr_2p.txt" \
  -c "${MCMC%.*}_corr_2p.txt" \
  -t "2p: MSA vs MCMC (${DESCRIPTION})" \
  -o "msa_mcmc_${SLUG}_2p"

echo "plotting 3p correlations"
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
