#!/bin/bash
set -euo pipefail

gtf=$1
peaks=$2
prefix=$3
output=$4
counts=$5

cd "$(dirname "$0")"

module load bedtools
module load R

bash scripts/overlap_calculations.sh $gtf $peaks $prefix
Rscript scripts/overlap_calculations_to_matrix.R $prefix $output $counts
