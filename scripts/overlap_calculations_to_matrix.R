#!/usr/local/bin Rscript
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
output_file <- args[2]
count_file <- args[3]
# count file not required
if (missing(count_file)){
  count_file <- NA
} else {count_file <- args[3]}

gtf <- read_tsv(paste0(prefix,'.gene.gtf'), col_names = F)
gene <- read_tsv(paste0(prefix,'.gene_peak_overlaps.data'), col_names = F)
exon <- read_tsv(paste0(prefix,'.exon_peak_overlaps.data'), col_names = F)
peak <- read_tsv(paste0(prefix, '.peaks_closest_genes.data'), col_names = F)

processor <- function(df){
  info_column_name <- paste0('X', grep('gene_id', df[1,]))
  count_column_name <- paste0('X', ncol(df))
  df_out <- df %>% 
    rowwise() %>% 
    ## extract ensgene and gene name
    mutate(ENSGene_Name = gsub(pattern = "gene_id||\\s+||\"", "", 
                               grep('gene_id', str_split(!!as.name(info_column_name), ';')[[1]], 
                                    value = T)),
           Gene_Name = gsub(pattern = "gene_name||\\s+||\"", "", 
                            grep('gene_name', str_split(!!as.name(info_column_name), ';')[[1]], 
                                 value = T))) %>% 
    select(ENSGene_Name, Gene_Name, bedtools_value = !!as.name(count_column_name)) 
  df_out
}
# extract all gene info
gtf_df <- processor(gtf)
# process gene based info
gene_df <- processor(gene) 
# process exon based info
exon_df <- processor(exon)
# join together, remove NA, count intron peaks
out_df <- left_join(gene_df %>% select(ENSGene_Name,
                                       Gene_Name,
                                       Peak_Count = bedtools_value), 
                    exon_df %>% select(ENSGene_Name, 
                                       Exon_Count = bedtools_value), 
                    by='ENSGene_Name') %>% unique()
out_df[is.na(out_df)] <- 0
out_df <- out_df %>% rowwise() %>% mutate(Intron_Count = max(0, Peak_Count - Exon_Count))

# process intergenic peaks
intergenic <- processor(peak %>% filter(grepl('gene_id', X13))) %>% 
  select(ENSGene_Name, Distance = bedtools_value)
# calculate peaks in each window, upstream and downstream of gene
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
  select(-Peak_Count)
out_df[is.na(out_df)] <- 0
# add gene length (TSS to TES)
gene$length <- abs(gene$X5-gene$X4)
gene$ENSGene_Name <- sapply(gene$X9, function(x)strsplit(x,'"') %>% unlist %>% .[[2]])

# if count file provided, then attach
if (!is.na(count_file)){
  counts <- read_csv(count_file,col_names = c("Gene_Name","Line","lsTPM" , "log2(lsTPM)" ,"Rank","TF" ))
  out_df <- left_join(out_df, gene[,c("ENSGene_Name",'length')],"ENSGene_Name")%>%
    left_join(counts[,c('Gene_Name','lsTPM')],by='Gene_Name')
  out_df$lsTPM[is.na(out_df$lsTPM)] <- 0
}

# finally add in missing genes (present in gtf, not in input, would be genes with 0 peaks near)
out_df <- left_join(gtf_df %>% select(ENSGene_Name, Gene_Name), out_df ) %>% arrange(Gene_Name)
out_df[is.na(out_df)] <- 0
write_tsv(out_df, path = output_file )
