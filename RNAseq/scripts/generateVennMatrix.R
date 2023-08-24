# Porthmeus
# 09.06.20

# create a venn matrix for the result tables

# logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(DESeq2)

# construc matrix
cln <- gsub(".RData","",basename(unlist(snakemake@input)))
fl_df <- load(snakemake@input[[1]])
df <- get(fl_df)
vennMat <- matrix(0, ncol = length(snakemake@input), nrow = nrow(df), dimnames = list(contig = rownames(df), comparisons = cln))

# create a matrix were the 1's in the matrix is replaced by the log2FC so that one can see the direction of regulation and the strength
vennDirMat <- vennMat

# read data and fill matrix
for(fl in unlist(snakemake@input)){
    fl_df <- load(fl)
    df <- get(fl_df)
    cln <- gsub(".RData","",basename(fl))
    sel <- df$padj <= 0.05
    sel[is.na(sel)] <- FALSE # remove the NAs introduced by statistical preselection
    vennMat[,cln] <- sel
    vennDirMat[sel,cln] <- df$log2FoldChange[sel]
}

# save the matrix
write.csv(file = snakemake@output[["venn"]], vennMat)
write.csv(file = snakemake@output[["vennDir"]], vennDirMat)
