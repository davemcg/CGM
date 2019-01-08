#!/bin/bash
set -euo pipefail

gtf=$1
peaks=$2
output=$3
counts=$4

cd "$(dirname "$0")"

prefix=`basename "$1"`_`basename "$2"`

bash scripts/overlap_calculations.sh $gtf $peaks $prefix
Rscript scripts/overlap_calculations_to_matrix.R $prefix $output $counts
