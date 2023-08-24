# Porthmeus
# 24.10.19

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(DESeq2)
require(ggplot2)
require(data.table)
require(limma)
require(svglite)

calcPCATable <- function(object, ntop=500){
    # takes an object from bioconducter which returns a matrix upon assay()
    # @object the data object
    # @ntop how many of the most variable genes should be used to calculate the variances

    #require(SummarizedExperiment)
    #rv <- apply(assay(object), 1, var)
    rv <- apply(object, 1, var)
    
    if(length(rv) < ntop){
        n <- length(rv)
    } else {
        n <- ntop
    }
    
    sel<- order(rv, decreasing=TRUE)[1:n]
    pca <- prcomp(t(object[sel,]), scale = TRUE)
    pca[["percent"]] <- round(pca$sdev^2/sum(pca$sdev^2),3)*100
    return(pca)

}


# load the DESeq2 matrix
load(snakemake@input[["deseq"]])

# transform the count matrix to log
rld <- rlog(deseq)
rldbe <- rld

# generate a second data set without batch effect
assay(rldbe) <- removeBatchEffect(assay(rld), batch = colData(rld)[["Batch"]], batch2 = colData(rld)[["Batch2"]])

# plot PCAs, MDAs and correlation plots
# calulate the PCA plot for the whole data set
pca <- calcPCATable(assay(rld))
df <- data.table(pca[["x"]])
df[,threshold := "thrld = 0"]
df[,batch:= "with batch"]
df <- as.data.frame(cbind(df, colData(rld)))

# do the same on the data without batch effect
pca <- calcPCATable(assay(rldbe))
df2 <- data.table(pca[["x"]])
df2[,threshold := "thrld = 0"]
df2[,batch:= "w/o batch"]
df2 <- as.data.frame(cbind(df2, colData(rld)))

# combine both data sets
df <- rbind(df, df2)

# apply different filter to low read count genes to check their influence on the plot
for(i in c(1,3,10)){
    sel <- apply(counts(deseq), 1, median) > i
    rldThr <- rlog(counts(deseq)[sel,])
    rldThrBe <- removeBatchEffect(rldThr, batch = colData(rld)[["Batch"]], batch2 = colData(rld)[["Batch2"]])
    # calc PCA for uncorrected data
    pca <- calcPCATable(rldThr)
    df2 <- data.table(pca[["x"]])
    df2[,threshold := paste("thrld =",i)]
    df2[,batch:= "with batch"]
    df2 <- as.data.frame(cbind(df2, colData(rld)))
    df <- rbind(df,df2)
    
    # calc PCA for batch corrected data
    pca <- calcPCATable(rldThrBe)
    df2 <- data.table(pca[["x"]])
    df2[,threshold := paste("thrld =",i)]
    df2[,batch:= "w/o batch"]
    df2 <- as.data.frame(cbind(df2, colData(rld)))
    df <- rbind(df,df2)
}

# create the plot
df[["threshold"]] <- factor(df[["threshold"]], level = paste("thrld =", c(0,1,3,10)))
pp <- ggplot(df, aes(x=PC1, y=PC2, color = Condition, label = Replicate)) + 
    geom_text() +
    ylab(paste("PC2 (", pca[["percent"]][2], "%)")) +
    xlab(paste("PC1 (", pca[["percent"]][1], "%)")) +
    facet_grid(batch~threshold)

ggsave(file = snakemake@output[["PCA"]], width = 32, height = 16, units = "cm", pp)

# save the data for later evaluation
write.csv(df, file = snakemake@output[["csv"]])
