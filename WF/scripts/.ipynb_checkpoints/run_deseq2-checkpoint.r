args = commandArgs(trailingOnly=TRUE)

# test if there are all the necessary arguments
if (length(args)!=5) {
  stop("5 arguments must be supplied: counts_file_path, metadata_file_path, up_value, outdir_path, targeted_gene_id", call.=FALSE)
}

require(dplyr)
require(stringr)
require(data.table)
require(DESeq2)
require(apeglm)

count_file_path = args[1]
metadata_file_path = args[2]
up_value = args[3]
outdir_path = args[4]
targeted_gene_id = args[5]

# helper function to parse input paths
trim_quotes_from_path<-function(path_string){
path_string = gsub("'", '', path_string)
path_string = gsub('"', '', path_string)
return(path_string)
}

# loading count matrix
cts <- read.csv(trim_quotes_from_path(count_file_path), sep='\t',header = TRUE,stringsAsFactors=FALSE,check.names = FALSE)
id_column <- colnames(cts)[1]
cts <- cts[!duplicated(cts[ , c(id_column)]),]

# loading metadata
metadata <- read.csv(trim_quotes_from_path(metadata_file_path), sep='\t',header = TRUE,stringsAsFactors=FALSE,check.names = FALSE)
colnames(metadata) <- c('sample','condition')
rownames(metadata) <- metadata$sample
metadata$condition <- (metadata$condition==up_value)

# prepare count matrix for DESEQ format
rownames(cts) <- cts[,id_column]
cts <- select(cts, rownames(metadata))

# create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design = ~ condition)

# VST COUNT TRANSFORMATION
vsd <- vst(dds, blind = FALSE)
vsd_table = assay(vsd)
vsd_table = data.frame(vsd_table)
cols <- colnames(vsd_table)
vsd_table <- cbind(rownames(vsd_table), vsd_table)
colnames(vsd_table) <- c(id_column,cols)
# write the VST transformed counts
write.table(vsd_table,file=paste(trim_quotes_from_path(outdir_path),"/VST_counts.tsv",sep = ""),row.names = FALSE, quote=FALSE, sep='\t')

# do not run Deseq2 if there is only one replicate per condition
if(nrow(metadata)==2){
    res <- vsd_table
    a = metadata[metadata$condition==TRUE,]['index'][[1]]
    b = metadata[metadata$condition==FALSE,]['index'][[1]]
    res$log2FoldChange <- log2(res[,a])-log2(res[,b])
    res$pvalue <- 1
    res$padj <- 1
    res$gene_name = rownames(res)
    res <- res[,c('gene_name','log2FoldChange','pvalue','padj')]
}else{
    dds <- DESeq(dds)
    res <- lfcShrink(dds, coef="conditionTRUE", type="apeglm")
    # res <- lfcShrink(dds, coef="conditionTRUE", apeAdapt=FALSE,type="apeglm")
    # res <- lfcShrink(dds, coef="conditionTRUE", type="normal")
    res = data.frame(res[order(res$padj),])
    res$gene_name = rownames(res)
    res <- res[,c('gene_name','baseMean','log2FoldChange','pvalue','padj')]
}

# write the log2FCs
write.table(res,file=paste(trim_quotes_from_path(outdir_path),"/DE_table.tsv",sep = ""),row.names = FALSE, quote=FALSE, sep='\t')
