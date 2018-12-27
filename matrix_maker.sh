#!/bin/bash
set -euo pipefail

gtf=$1
peaks=$2
#prefix=$3
prefix=`basename "$1"`_`basename "$2"`
output=$3
counts=$4

cd "$(dirname "$0")"

#module load bedtools
#module load R

bash scripts/overlap_calculations.sh $gtf $peaks $prefix
Rscript scripts/overlap_calculations_to_matrix.R $prefix $output $counts
