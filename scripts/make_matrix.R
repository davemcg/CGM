library(dplyr)
library(GenomicRanges)
library(gread)
library(parallel)
args=commandArgs(trailingOnly = T)
# gene='OR4F5'
#targ_peaks <- 'data/IPSC_ATAC.common_peaks.blackListed.bed'
targ_peaks <- args[1]
target_peaks=read.table(targ_peaks, stringsAsFactors = F,sep = '\t')
colnames(target_peaks) <- c('p_chr','p_start','p_end','peak')
# gread library lets you add introns into a gtf easily
gtf <- read_format('data/gencode.v27.basic.annotation.gtf')
gtf_df <- as.data.frame(gtf)
all_genes <- unique(gtf_df$gene_name)
all_genes_plus <- filter(gtf_df, strand=='+')[,'gene_name']%>%unique
getOption("mc.cores", 2L)
make_matrix <- function(gtf,gene,target_peaks){
    #a=Sys.time()
    gtf_gene <- gtf[elementMetadata(gtf)[,'gene_name']==gene]
    gtf_gene_introns <- construct_introns(gtf_gene, update=TRUE)%>%as.data.frame
    
    if(any(grepl('appris', gtf_gene_introns[,'tag']))){
        # there is an appris tag for this gene
        app_ids <- gtf_gene_introns[,'tag']%>%grep('appris',.)%>% gtf_gene_introns[.,'tag']%>% unique
        max=strsplit(app_ids,'_')%>%lapply(function(x)x[[3]])%>%unlist%>%which.max
        gtf_target_tx <- filter(gtf_gene_introns, tag==app_ids[max])
    }else if(any(gtf_gene_introns[,'transcript_type']%in%c('protein_coding'))){
        #pick the first protein coding transcript
        tx=filter(gtf_gene_introns,transcript_type=='protein_coding')[,'transcript_id']%>%unique%>%na.omit%>%.[[1]]
        gtf_target_tx <- filter(gtf_gene_introns, transcript_id==tx)
    }else{
       #pick the first transcript
        tx <- unique(gtf_gene_introns[,'transcript_id'])%>%na.omit%>%.[[1]]
        gtf_target_tx <- filter(gtf_gene_introns, transcript_id==tx)
    }
    
    ex_int_bed <- filter(gtf_target_tx, feature%in%c('exon','intron'))%>%select( seqnames, start, end, feature,score,strand)
    #now lets make the distance bins
    # struct-1000000 -100000 -10000 -1000 TSS]  [TES +1000 +10000 +100000 +1000000
    #b=Sys.time()
    #c=Sys.time()
    make_bins <- function(TSS, TES, locs, strand,chr){
        last=0
        upres=list()
        downres=list()
        i=1
        for( loc in locs){
            downres[[i]] <- data.frame(seqnames=chr,start=TES+last,end=TES+loc,feature=paste0('ds',loc/1000,'k'),score=NA,strand=strand)
            upres[[i]] <- data.frame(seqnames=chr,start=TSS-loc,end= TSS-last,feature=paste0('us',loc/1000,'k'),score=NA,strand=strand)
            last=loc
            i <- i+1
        }
        #account for genes at the begnining of a chrom; currently nor accounting for gens at the end of a chrom
        upres_valid <- (lapply(upres, function(x) x[,'end']>0)%>%unlist)%>%upres[.]
        if(length(upres_valid)<length(upres)) upres_valid[[length(upres_valid)]][,'start'] <- 0
        return(do.call(rbind,c(upres_valid,downres)))
    }
    TSS=filter(gtf_target_tx, feature=='transcript')[,'start']
    TES=filter(gtf_target_tx, feature=='transcript')[,'end']
    strand=filter(gtf_target_tx, feature=='transcript')[,'strand']%>%as.character()
    chr=filter(gtf_target_tx, feature=='transcript')[,'seqnames']%>%as.character()
    locs <- c(1000, 10000, 100000,1000000)
    all_features_bed <- rbind(ex_int_bed, make_bins(TSS,TES,locs,strand,chr),stringsAsFactors=F)
    
    all_features_range=GRanges(seqnames=all_features_bed$seqnames, ranges=IRanges(start=all_features_bed$start,end=all_features_bed$end),      
                               strand = all_features_bed$strand, feature=all_features_bed$feature)
    
    target_peaks_chr <- filter(target_peaks,p_chr%in%chr )
    if(nrow(target_peaks_chr)==0){
        print(paste0('chrom for ',gene,' not found in peak file'))
        loc_names <- c(paste0('ds',locs/1000,'k'),paste0('us',locs/1000,'k'),'exon','intron','expression')
        res <- matrix(data=0, nrow = 1, ncol = length(loc_names))
        rownames(res) <- gene
        colnames(res) <- loc_names
        return(res)
        
    }
    #d=Sys.time()
    #e=Sys.time()
    peaks_range <- GRanges(seqnames=target_peaks_chr$p_chr, ranges = IRanges(start=target_peaks_chr$p_start,end = target_peaks_chr$p_end),peaks=target_peaks_chr$peak)    
    hits <- findOverlaps(all_features_range,peaks_range)
    overlaps <- pintersect(all_features_range[queryHits(hits)], peaks_range[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width( peaks_range[subjectHits(hits)])
    hits <- hits[percentOverlap > 0.5]
    joined_hits <- cbind(all_features_bed[from(hits),], target_peaks[to(hits),])
    if(nrow(joined_hits)==0){
        # why are you like this 
        print(c(gene,' found no peaks'))
        loc_names <- c(paste0('ds',locs/1000,'k'),paste0('us',locs/1000,'k'),'exon','intron','expression')
        res <- matrix(data=0, nrow = 1, ncol = length(loc_names))
        rownames(res) <- gene
        colnames(res) <- loc_names
        return(res)
    }
    combined_peaks <- aggregate(joined_hits[,"peak"], by=list(joined_hits$feature), sum)
    join_frame <- data.frame(t(combined_peaks$x))
    colnames(join_frame) <- combined_peaks$Group.1
    loc_names <- c(paste0('ds',locs/1000,'k'),paste0('us',locs/1000,'k'),'exon','intron','expression')
    res <- matrix(data=0, nrow = 1, ncol = length(loc_names))
    colnames(res) <- loc_names
    res[,colnames(join_frame)] <- join_frame[1,]%>%as.numeric 
    rownames(res) <- gene
    #f=Sys.time()
    return(res)
}


k=Sys.time()
all_plus_genes <- mclapply(all_genes,function(x) make_matrix(gtf,x,target_peaks),mc.cores = args[2])
m=Sys.time()
m-k

#current runtime for 


