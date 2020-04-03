#! /bin/bash
set -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

OPTIONS=s:c:S:C:d:
LONGOPTIONS=msa:,msa_label:,mcmc:,mcmc_label:,description:
PARSED=$(getopt -o ${OPTIONS} --long ${LONGOPTIONS} -n "$0" -- "$@")
eval set -- "$PARSED"

while [ $# -ge 1 ]; do
  case "$1" in
    -s|--msa)
      MSA="$2"
      shift 2
      ;;
    -S|--msa_label)
      MSA_LABEL="$2"
      shift 2
      ;;
    -c|--mcmc)
      MCMC="$2"
      shift 2
      ;;
    -C|--mcmc_label)
      MCMC_LABEL="$2"
      shift 2
      ;;
    -d|--description)
      DESCRIPTION="$2"
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

SLUG=$(echo "$MSA_LABEL $MCMC_LABEL $DESCRIPTION" | \
       tr '[:upper:]' '[:lower:]' | \
       sed -e "s/ /_/g" | \
       sed -e "s/,//g" -e "s/(//g" -e "s/)//g")

echo "computing pairwise distances for $MSA"
if [[ "${MSA}" =~ "numeric" ]]; then
  "${SCRIPT_DIR}/compute_distances" \
    -n "$MSA" \
    -o "${MSA%.*}_distances.txt"
else
  "${SCRIPT_DIR}/compute_distances" \
    -i "$MSA" \
    -o "${MSA%.*}_distances.txt"
fi

echo "computing pairwise distances for $MCMC"
if [[ "${MCMC}" =~ "numeric" ]]; then
  "${SCRIPT_DIR}/compute_distances" \
    -n "$MCMC" \
    -o "${MCMC%.*}_distances.txt"
else
  "${SCRIPT_DIR}/compute_distances" \
    -i "$MCMC" \
    -o "${MCMC%.*}_distances.txt"
fi

echo "plotting distances"
"${SCRIPT_DIR}/plot_distances.py" \
  -s "${MSA%.*}_distances.txt" \
  -c "${MCMC%.*}_distances.txt" \
  -S "${MSA_LABEL}" \
  -C "${MCMC_LABEL}" \
  -t "MSA vs MCMC (${DESCRIPTION})" \
  -o "${SLUG}.svg"

# clean up
rm -vf "${MSA%.*}_distances.txt"
rm -vf "${MCMC%.*}_distances.txt"
