#!/bin/bash

gtf=$1
peaks=$2
prefix=$3

# cut down gtf to only include gene gtf that are protein coding
gzcat $gtf | \
	awk '$3 == "gene" {print $0}' | \
	grep "protein_coding" | \
	bedtools sort -i - > /tmp/$prefix.gene.gtf

# count peaks overlapping each gtf gene
bedtools intersect -a /tmp/$prefix.gene.gtf -b $peaks -F 0.5 -c > $prefix.gene_peak_overlaps.data
# count peaks overlapping each gtf gene exon (can infer intron overlap by subtracting gene_peak - exon_peak
gzcat $gtf | \
	awk '$3 == "exon" {print $0}' | \
	grep "protein_coding" | \
	bedtools merge -i - | \
	bedtools intersect -a - -b $peaks -F 0.5 -c | \
	awk '$NF != 0 {print $0}' > $prefix.exon_peak_overlaps.data

# find 10 closest genes to each peak (removing peaks that overlap genes)
# remove peak <-> gene where distance > 1 mb
bedtools closest -a <( bedtools intersect -a $peaks -b /tmp/$prefix.gene.gtf -v ) \
	-b /tmp/$prefix.gene.gtf -k 5  -D "a" -io | \
	awk '$NF > -1000001 && $NF < 1000001 {print $0}' > $prefix.peaks_closest_genes.data

# remove gtf gene file
rm /tmp/$prefix.gene.gtf
