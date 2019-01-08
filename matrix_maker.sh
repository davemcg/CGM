#!/bin/bash
set -euo pipefail

gtf=$1
peaks=$2
output=$3
counts=$4

cd "$(dirname "$0")"

module load bedtools
module load R

prefix=$1_$2

bash scripts/overlap_calculations.sh $gtf $peaks $prefix
Rscript scripts/overlap_calculations_to_matrix.R $prefix $output $counts
