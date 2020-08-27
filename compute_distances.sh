#! /bin/bash
set -eu

OPTIONS=i:n:l:I:N:L:d:
LONGOPTIONS=msa:,msa_label:,mcmc:,mcmc_label:,description:
PARSED=$(getopt -o ${OPTIONS} --long ${LONGOPTIONS} -n "$0" -- "$@")
eval set -- "$PARSED"

while [ $# -ge 1 ]; do
  case "$1" in
    -i|--msa)
      MSA="$2"
      shift 2
      ;;
    -n|--msa_numeric)
      MSA="$2"
      MSA_NUMERIC=true
      shift 2
      ;;
    -l|--msa_label)
      MSA_LABEL="$2"
      shift 2
      ;;
    -I|--mcmc)
      MCMC="$2"
      shift 2
      ;;
    -N|--mcmc_numeric)
      MCMC="$2"
      MCMC_NUMERIC=true
      shift 2
      ;;
    -L|--mcmc_label)
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
if ! [[ "${MSA_NUMERIC}" == true ]]; then
  compute_distances \
    -i "$MSA" \
    -o "${MSA%.*}_distances.txt"
else
  compute_distances \
    -n "$MSA" \
    -o "${MSA%.*}_distances.txt"
fi

echo "computing pairwise distances for $MCMC"
if ! [[ "${MCMC_NUMERIC}" == true ]]; then
  compute_distances \
    -i "$MCMC" \
    -o "${MCMC%.*}_distances.txt"
else
  compute_distances \
    -n "$MCMC" \
    -o "${MCMC%.*}_distances.txt"
fi

echo "plotting distances"
plot_distances.py \
  -s "${MSA%.*}_distances.txt" \
  -c "${MCMC%.*}_distances.txt" \
  -S "${MSA_LABEL}" \
  -C "${MCMC_LABEL}" \
  -t "MSA vs MCMC (${DESCRIPTION})" \
  -o "${SLUG}.svg"

# clean up
rm -vf "${MSA%.*}_distances.txt"
rm -vf "${MCMC%.*}_distances.txt"
