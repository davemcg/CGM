#!/bin/bash
set -euo pipefail

gtf=$1
peaks=$2
prefix=`basename "$1"`_`basename "$2"`

# cut down gtf to only include gene gtf that are protein coding
gzcat $gtf | \
	awk '$3 == "gene" {print $0}' | \
	grep "protein_coding" | \
	bedtools sort -i - > /tmp/$prefix.gene.gtf

# return peaks overlapping each gtf gene
bedtools intersect -a /tmp/$prefix.gene.gtf -b $peaks -F 0.5 -wa -wb > /tmp/$prefix.gene_peak_overlaps.data
# return peaks overlapping each gtf gene exon (can infer intron overlap by subtracting gene_peak - exon_peak
gzcat $gtf | \
	awk '$3 == "exon" {print $0}' | \
	grep "protein_coding" | \
	bedtools sort -i - | \
	bedtools merge -s -c 9 -o "first" -i - | 
	bedtools intersect -a - -b $peaks -F 0.5 -wa -wb | \
	awk '$NF != 0 {print $0}' > /tmp/$prefix.exon_peak_overlaps.data

# find 10 closest genes to each peak 
# remove peak <-> gene where distance > 1 mb
# bedtools closest -a <( bedtools intersect -a $peaks -b /tmp//tmp/$prefix.gene.gtf -v ) \ # (removing peaks that overlap genes)
bedtools closest -a $peaks -b /tmp/$prefix.gene.gtf -k 5  -D "a" -io | \
	awk '$NF > -1000001 && $NF < 1000001 {print $0}' > /tmp/$prefix.peaks_closest_genes.data


