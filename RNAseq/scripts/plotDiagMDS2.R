# Porthmeus
# 24.10.19
# redirect messages to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(DESeq2)
require(limma)
require(data.table)
require(ggplot2)
require(MASS)

# load the DESeq2 matrix
getDeseq <- load(snakemake@input[["deseq"]])
deseq <- get(getDeseq)

# transform data
rld <- rlog(deseq)

# remove the batch effect
rldBR <- rld
assay(rldBR) <- removeBatchEffect(assay(rld), batch = colData(rld)[["Batch"]], batch2 = colData(rld)[["Batch2"]])

# define the thresholds
thr  <- c(0,1,3,10)

# for each threshold get the most variable genes, compute the most variable genes, transform the data and calculate the distance and iterative multidimensional scaling
df <- data.frame(character(), character(), integer(), numeric(), numeric())
stress_df <- data.frame(character(), integer(), numeric())
for(th in thr){
    # filter for the different thresholds
    sel <- apply(assay(deseq), 1, median) > th
    
    #extract 500 most variable genes
    selBR <- names(sort(apply(assay(rldBR[sel,]), 1, var), decreasing = TRUE)[1:500])
    sel <- names(sort(apply(assay(rld[sel,]), 1, var), decreasing = TRUE)[1:500])
    
    # calculate distances
    mat <- assay(rld)[sel,]
    colnames(mat) <- colData(rld)[["Replicate"]]
    dst <- dist(t(mat))
    
    # also for the matrix without batch
    matBR <- assay(rldBR)[selBR,]
    colnames(matBR) <- colData(rldBR)[["Replicate"]]
    dstBR <- dist(t(matBR))
    
    # calculate the iterative scaling
    mds <- isoMDS(dst, k =2)
    mdsBR <- isoMDS(dstBR, k =2)
    
    # combine everything into a data frame
    df <- rbind(df, data.frame(Replicate = c(row.names(mds$points), 
                                             row.names(mdsBR$points)), 
                               Batch = c(rep("with Batch",ncol(deseq)),
                                         rep("w/o Batch",ncol(deseq))), 
                               Threshold = rep(th, 2*ncol(deseq)), 
                               PC1 = c(mds$points[,1],
                                       mdsBR$points[,1]),
                               PC2 = c(mds$points[,2],
                                       mdsBR$points[,2]),
                               Stress = c(rep(mds$stress, ncol(deseq)),
                                          rep(mdsBR$stress, ncol(deseq)))
    ))
    
    
    
    # save the stress level in another data frame for later evaluation
    stress_df <- rbind(stress_df, data.frame(Batch = c("with Batch","w/o Batch"),
                                             Threshold = c(th,th),
                                             Stress = c(mds$stress,mdsBR$stress)))
}
# merge the data with the meta data
df <- merge(data.table(df), data.table(as.data.frame(colData(deseq))), by ="Replicate")

# generate a mds plot from the results
pMDS <- ggplot(df, aes(x=PC1, y=PC2, color = Condition, label = Replicate)) +
    geom_text() +
    facet_grid(Batch.x~Threshold)
ggsave(pMDS, 
       filename = snakemake@output[["mds"]],
       width = 32,
       height = 16,
       units = "cm",
       device ="svg")

# generate barplot from the stress values
pStress <- ggplot(stress_df, aes(x=as.factor(Threshold), y=Stress)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Batch)
ggsave(pStress, 
       filename = snakemake@output[["stress"]],
       width = 32,
       height = 16,
       units = "cm",
       device ="svg")

# save the data frame for future analysis
write.csv(df, file = snakemake@output[["csv"]])