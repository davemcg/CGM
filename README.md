# CGM
Chromatin Gene Matrix

Collapse chromatin information (e.g. ChIP-seq data) into a well-behaving matrix centered around each gene. 

Input is a bed file for each condition/sample/etc (scores optional) along with a genome reference matched (we use gencode) GTF.

Normalized RNA-seq TPM/RPKM quantification for each gene in the GTF can be given as optional input. 

Output is a matrix where each row is a gene with columns for chromatin signal in windows upstream, downstream, gene expression (optional), exonic, and intronic. 
