# Porthmeus
# 08.10.19
# redirect messages to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(DESeq2)
require(limma)
require(RColorBrewer)
require(pheatmap)


plotCorHeatmap <- function(mat, main = NA, names = NULL){
    # helper function to plot the pheatmap
    sampleDists <- dist(t(mat))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- names
    colnames(sampleDistMatrix) <- names
    
    pheatmap(sampleDistMatrix,
        col=colors,
        main = main,
        trace = F)
}
# load the DESeq2 matrix
load(snakemake@input[["deseq"]])

# define colors
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

names <- colData(deseq)[["Replicate"]]

# plot the plots... no grid possible, as to it would become quite complicated to achieve it.
for(i in c(0,1,3,10)){
    sel <- apply(counts(deseq), 1, median) > i
    rld <- rlog(counts(deseq)[sel,])
    for(batch in c("wb","wob")){
        if(batch == "wb"){
            svg(file = snakemake@output[[paste0(batch, i)]], width = 16/2.54, height =16/2.54)
            plotCorHeatmap(rld, names = names, main = paste("with batch, thrld =",i))
            dev.off()
        } else if(batch == "wob"){
            rldbr <- removeBatchEffect(rld,
                                       batch = colData(deseq)[["Batch"]],
                                       batch2 = colData(deseq)[["Batch2"]])
            svg(file = snakemake@output[[paste0(batch, i)]], width = 16/2.54, height =16/2.54)
            plotCorHeatmap(rldbr, names = names, main = paste("w/o batch, thrld =",i))
            dev.off()
        }
    }
}

