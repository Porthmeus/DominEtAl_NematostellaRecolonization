# Porthmeus
# 24.10.19

log <- file(snakemake@log[[1]], "wt")
sink(log)
sink(log, type = "message")

require(DESeq2)
require(data.table)
require(pheatmap)
require(RColorBrewer)


# read data
load(snakemake@input[["deseq"]])
res <- fread(snakemake@input[["res"]])

# get the top genes
sigres <- res[padj <= 0.05,]
if(nrow(sigres) >= 40){
    upper <- 40
} else {
    upper <- nrow(sigres)
}

top40 <- sigres[order(abs(log2FoldChange), decreasing = T)[1:upper],V1]

if(nrow(sigres[log2FoldChange > 0,]) >= 20){
    upper <- 20
} else {
    upper <- nrow(sigres[log2FoldChange > 0,])
}

top20up <- sigres[order(log2FoldChange, decreasing = T)[1:upper],V1]

if(nrow(sigres[log2FoldChange < 0,]) >= 20){
    upper <- 20
} else {
    upper <- nrow(sigres[log2FoldChange < 0,])
}

top20down <- sigres[order(log2FoldChange, decreasing = F)[1:upper],V1]



# plot a heatmap of the genes
allgenes <- unique(c(top40,top20up,top20down))
rld <- counts(deseq, normalized=TRUE)[allgenes,] # use counts normalized for sequencing depth
#rld <- rlog(assay(deseq)[allgenes,])
colnames(rld) <- paste(colData(deseq)[["Condition"]], colData(deseq)[["Plate"]])
rownames(rld) <- allgenes


cols <- colorRampPalette(c("blue","white","red"))(255)
svg(file = snakemake@output[["top"]], width = 16/2.54, height = 16/2.54)
pheatmap(rld[top40,], col = cols, scale = "row" )
dev.off()
svg(file = snakemake@output[["up"]], width = 16/2.54, height = 16/2.54)
pheatmap(rld[top20up,], col = cols, scale = "row")
dev.off()
svg(file = snakemake@output[["down"]], width = 16/2.54, height = 16/2.54)
pheatmap(rld[top20down,], col = cols, scale = "row")
dev.off()
