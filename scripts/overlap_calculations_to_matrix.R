#!/usr/local/bin Rscript
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
output_file <- args[2]
score_col <- args[3] 
# this is the column number in the peak bed file with a 'score'
# 1 - based number
# for example chr1 100 200 45 would be 4 
# this should be some kind of 'score' for 
# the peak where a higher number is a 'stronger' peak
# like -log10(p-value) or reads overlapping
# if not given, then just use peak counts
if (missing(score_col)){
	score_col <- NA
} else {score_col <- as.numeric(args[3])}

#count_file <- args[3]
# count file not required
#if (missing(count_file)){
#  count_file <- NA
#} else {count_file <- args[3]}

gtf <- read_tsv(paste0(prefix,'.gene.gtf'), col_names = F)
gene <- read_tsv(paste0(prefix,'.gene_peak_overlaps.data'), col_names = F)
exon <- read_tsv(paste0(prefix,'.exon_peak_overlaps.data'), col_names = F)
peak <- read_tsv(paste0(prefix, '.peaks_closest_genes.data'), col_names = F)

processor <- function(df, score_col = NA){
  info_column_name <- paste0('X', grep('gene_id', df[1,]))
  if (!is.na(score_col)){
	if (score_col == 'last') {
		score_column_name <- paste0('X', ncol(df))
	} else {
  		score_column_name <- paste0('X', grep('gene_id', df[1,]) + score_col) 
	}
  }
  #count_column_name <- paste0('X', ncol(df))
  if (is.na(score_col)){
  df_out <- df %>% 
    group_by(!!as.name(info_column_name)) %>% 
    summarise(Value = n()) %>% 
    ungroup() %>% 
    rowwise() %>% 
    ## extract ensgene and gene name
    mutate(ENSGene_Name = gsub(pattern = "gene_id||\\s+||\"", "", 
                               grep('gene_id', str_split(!!as.name(info_column_name), ';')[[1]], 
                                    value = T)),
           Gene_Name = gsub(pattern = "gene_name||\\s+||\"", "", 
                            grep('gene_name', str_split(!!as.name(info_column_name), ';')[[1]], 
                                 value = T))) %>%
 	ungroup() %>%  
    select(ENSGene_Name, Gene_Name, Value) 
  } else if (!is.na(score_col) & score_col != 'last') {
  df_out <- df %>% 
    group_by(!!as.name(info_column_name)) %>% 
    summarise(Value = sum(!!as.name(score_column_name))) %>% 
    ungroup() %>% 
    rowwise() %>% 
    ## extract ensgene and gene name
    mutate(ENSGene_Name = gsub(pattern = "gene_id||\\s+||\"", "", 
                               grep('gene_id', str_split(!!as.name(info_column_name), ';')[[1]], 
                                    value = T)),
           Gene_Name = gsub(pattern = "gene_name||\\s+||\"", "", 
                            grep('gene_name', str_split(!!as.name(info_column_name), ';')[[1]], 
                                 value = T))) %>% 
	ungroup() %>% 
    select(ENSGene_Name, Gene_Name, Value) 
  } else {
    df_out <- df %>%
      rowwise() %>%
      ## extract ensgene and gene name
      mutate(ENSGene_Name = gsub(pattern = "gene_id||\\s+||\"", "", 
                                 grep('gene_id', str_split(!!as.name(info_column_name), ';')[[1]], 
                                      value = T)),
             Gene_Name = gsub(pattern = "gene_name||\\s+||\"", "", 
                              grep('gene_name', str_split(!!as.name(info_column_name), ';')[[1]], 
                                   value = T))) %>% 
      select(ENSGene_Name, Gene_Name, Value = !!as.name(score_column_name))
  }
  df_out
}
# extract all gene info
print('gtf')
gtf_df <- processor(gtf, NA)
# process gene based info
print('gene')
gene_df <- processor(gene, score_col) 
# process exon based info
print('exon')
exon_df <- processor(exon, score_col) %>% 
	group_by(ENSGene_Name) %>% 
	summarise(Gene_Name = max(Gene_Name), Value = sum(Value)) %>% 
	ungroup()
# join together, remove NA, count intron peaks
out_df <- left_join(gene_df %>% select(ENSGene_Name,
                                       Gene_Name,
                                       Peak_Value = Value), 
                    exon_df %>% select(ENSGene_Name, 
                                       Exon_Value = Value), 
                    by='ENSGene_Name') %>% unique()
out_df[is.na(out_df)] <- 0
out_df <- out_df %>% rowwise() %>% mutate(Intron_Value = max(0, Peak_Value - Exon_Value))

# process intergenic peaks
# find gene info column first
gene_col <- colnames(peak)[grepl('gene_id', peak[1,])]
print('peak')
intergenic <- processor(peak %>% filter(grepl('gene_id', !!as.name(gene_col))), 'last') %>% 
  select(ENSGene_Name, Distance = Value)
# calculate peaks in each window, upstream and downstream of gene
if (is.na(score_col)){
out_df <- intergenic %>% 
  ungroup() %>% 
  group_by(ENSGene_Name) %>% 
  summarise(`1e6` = sum(Distance > 5e5 & Distance <= 1e6),
            `5e5` = sum(Distance > 1e5 & Distance <= 5e5),
            `1e5` = sum(Distance > 5e4 & Distance <= 1e5),
            `5e4` = sum(Distance > 1e4 & Distance <= 5e4),
            `1e4` = sum(Distance > 5000 & Distance <= 1e4),
            `5000` = sum(Distance > 1000 & Distance <= 5000),
            `1000` = sum(Distance > 0 & Distance <= 1000),
            `-1000` = sum(Distance < 0 & Distance >= -1000),
            `-5000` = sum(Distance < -1000 & Distance >= -5000),
            `-1e4` = sum(Distance < -5000 & Distance >= -1e4),
            `-5e4` = sum(Distance < -1e4 & Distance >= -5e4),
            `-1e5` = sum(Distance < -5e4 & Distance >= -1e5),
            `-5e5` = sum(Distance < -1e5 & Distance >= -5e5),
            `-1e6` = sum(Distance < -5e5 & Distance >= -1e6)) %>% 
  left_join(out_df, .) %>% 
  select(-Peak_Value)
} else {
print('grouping')
intergenic <- cbind(intergenic, peak %>% filter(grepl('gene_id', !!as.name(gene_col))) %>% select(Value = X4))
out_df <- intergenic %>% 
  ungroup() %>% 
  group_by(ENSGene_Name) %>% 
  summarise(`1e6` = sum(Value[Distance > 5e5 & Distance <= 1e6]),
            `5e5` = sum(Value[Distance > 1e5 & Distance <= 5e5]),
            `1e5` = sum(Value[Distance > 5e4 & Distance <= 1e5]),
            `5e4` = sum(Value[Distance > 1e4 & Distance <= 5e4]),
            `1e4` = sum(Value[Distance > 5000 & Distance <= 1e4]),
            `5000` = sum(Value[Distance > 1000 & Distance <= 5000]),
            `1000` = sum(Value[Distance > 0 & Distance <= 1000]),
            `-1000` = sum(Value[Distance < 0 & Distance >= -1000]),
            `-5000` = sum(Value[Distance < -1000 & Distance >= -5000]),
            `-1e4` = sum(Value[Distance < -5000 & Distance >= -1e4]),
            `-5e4` = sum(Value[Distance < -1e4 & Distance >= -5e4]),
            `-1e5` = sum(Value[Distance < -5e4 & Distance >= -1e5]),
            `-5e5` = sum(Value[Distance < -1e5 & Distance >= -5e5]),
            `-1e6` = sum(Value[Distance < -5e5 & Distance >= -1e6])) %>% 
  left_join(out_df, .) %>% 
  select(-Peak_Value)

}
out_df[is.na(out_df)] <- 0
# add gene length (TSS to TES)
print('gene length')
gene$length <- abs(gene$X5-gene$X4)
gene$ENSGene_Name <- sapply(gene$X9, function(x)strsplit(x,'"') %>% unlist %>% .[[2]])

# if count file provided, then attach
#if (!is.na(count_file)){
#  counts <- read_csv(count_file,col_names = c("Gene_Name","Line","lsTPM" , "log2(lsTPM)" ,"Rank","TF" ))
#  out_df <- left_join(out_df, gene[,c("ENSGene_Name",'length')],"ENSGene_Name")%>%
#    left_join(counts[,c('Gene_Name','lsTPM')],by='Gene_Name')
#  out_df$lsTPM[is.na(out_df$lsTPM)] <- 0
#}

# finally add in missing genes (present in gtf, not in input, would be genes with 0 peaks near)
print('missing')
out_df <- left_join(gtf_df %>% select(ENSGene_Name, Gene_Name), out_df ) %>% 
  left_join(., gene %>% select(ENSGene_Name, Gene_length = length)) %>%  # add length info
  arrange(Gene_Name)
out_df[is.na(out_df)] <- 0
out_df <- out_df %>% unique()
write_tsv(out_df, path = output_file )
