#! /bin/bash
set -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
REWEIGHT=false

OPTIONS=s:c:d:t:r
LONGOPTIONS=msa:,mcmc:,description:,threshold:,reweight
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
    -r|--reweight)
      REWEIGHT=true
      shift
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

if [ "${REWEIGHT}" == true ]; then
  "${SCRIPT_DIR}/compare_stats" -r \
    -s "$MSA" \
    -c "$MCMC" \
    -t "$THRESHOLD"
else
  "${SCRIPT_DIR}/compare_stats" \
    -s "$MSA" \
    -c "$MCMC" \
    -t "$THRESHOLD"
fi

echo "plotting 1p frequencies"
"${SCRIPT_DIR}/plot_stats.py" \
  -s "${MSA%.*}_freq_1p.txt" \
  -c "${MCMC%.*}_freq_1p.txt" \
  -t "1p frequencies (${DESCRIPTION})" \
  -l \
  -o "${SLUG}_freq_1p"

echo "plotting 2p frequencies"
"${SCRIPT_DIR}/plot_stats.py" \
  -s "${MSA%.*}_freq_2p.txt" \
  -c "${MCMC%.*}_freq_2p.txt" \
  -t "2p frequencies (${DESCRIPTION})" \
  -l \
  -o "${SLUG}_freq_2p"

echo "plotting 2p correlations"
"${SCRIPT_DIR}/plot_stats.py" \
  -s "${MSA%.*}_corr_2p.txt" \
  -c "${MCMC%.*}_corr_2p.txt" \
  -t "2p correlations (${DESCRIPTION})" \
  -l \
  -o "${SLUG}_corr_2p"

if (( $(echo "$THRESHOLD == 0" | bc -l) )); then
  echo "plotting 3p frequencies"
  "${SCRIPT_DIR}/plot_stats.py" \
    -s "${MSA%.*}_freq_3p.txt" \
    -c "${MCMC%.*}_freq_3p.txt" \
    -t "3p frequencies (${DESCRIPTION})" \
    -l \
    -o "${SLUG}_freq_3p"

  echo "plotting 3p correlations"
  "${SCRIPT_DIR}/plot_stats.py" \
    -s "${MSA%.*}_corr_3p.txt" \
    -c "${MCMC%.*}_corr_3p.txt" \
    -t "3p correlations (${DESCRIPTION})" \
    -l \
    -o "${SLUG}_corr_3p"
else
  echo "plotting 3p frequencies"
  "${SCRIPT_DIR}/plot_stats.py" \
    -s "${MSA%.*}_freq_3p.txt" \
    -c "${MCMC%.*}_freq_3p.txt" \
    -t "3p frequencies (${DESCRIPTION})" \
    -o "${SLUG}_freq_3p"

  echo "plotting 3p correlations"
  "${SCRIPT_DIR}/plot_stats.py" \
    -s "${MSA%.*}_corr_3p.txt" \
    -c "${MCMC%.*}_corr_3p.txt" \
    -t "3p correlations (${DESCRIPTION})" \
    -o "${SLUG}_corr_3p"
fi

rm -v \
  "${MSA%.*}_freq_1p.txt" \
  "${MCMC%.*}_freq_1p.txt" \
  "${MSA%.*}_freq_2p.txt" \
  "${MCMC%.*}_freq_2p.txt" \
  "${MSA%.*}_corr_2p.txt" \
  "${MCMC%.*}_corr_2p.txt" \
  "${MSA%.*}_freq_3p.txt" \
  "${MCMC%.*}_freq_3p.txt" \
  "${MSA%.*}_corr_3p.txt" \
  "${MCMC%.*}_corr_3p.txt"
