#!/bin/bash
set -euo pipefail

gtf=$1
peaks=$2
output=$3
score_column=$4
#if [ ! -z $4]
#then
#	counts=$4
#else
#	count=1
#fi

#cd "$(dirname "$0")"

prefix=`basename "$1"`_`basename "$2"`

bash ~/git/CGM/scripts/overlap_calculations.sh $gtf $peaks $prefix
Rscript ~/git/CGM/scripts/overlap_calculations_to_matrix.R /tmp/$prefix $output $score_column 
