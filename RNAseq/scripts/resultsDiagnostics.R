# Porthmeus
# 24.10.2019

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(DESeq2)

# load results tables
load(snakemake@input[[1]])

# plot vulcano plot
jpeg(file = snakemake@output[["MA"]], width = 16/2.54, height = 16/2.54, units = "in", res=300)
plotMA(df, ylim =c(-1,1))
dev.off()

# plot histogram of p-values, baseMeans and foldChanges
svg(file = snakemake@output[["Hists"]], width = 32/2.54, height = 16/2.54)
par(mfcol = c(1,3))
hist(df[["padj"]], main = "Histogram of adjusted p-values", xlab = "p-value")
hist(abs(df[["log2FoldChange"]][df[["padj"]] <= 0.5]), main = "Histogram of significant log2(FoldChanges)", xlab = "log2(FC)")
hist(log10(df[["baseMean"]][df[["padj"]] <= 0.5]), main = "Histogram of mean expression", xlab = "log10(baseMean)")
dev.off()


# plot the ratios of rejected H0 in dependence of mean normalized counts per gene
########
# code snippet adapted from http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2 (10.10.19)
# create bins using the quantile function
qs <- c( 0, quantile( df[["baseMean"]][df[["baseMean"]] > 0], 0:10/10 ) )
# "cut" the genes into the bins
bins <- cut(df$baseMean, qs)
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(0.5*qs[-1] + 0.5*qs[-length(qs)]))
# calculate the ratio of p values less than .01 for each bin
ratios <- tapply( df$pvalue, bins, function(p) mean( p <= 0.05, na.rm=TRUE ) )
# plot these ratios
svg(file = snakemake@output[["SmallP"]], width = 16/2.54, height = 16/2.54)
barplot(ratios, xlab="mean normalized count", ylab="ratio of p-values <= 0.05")
dev.off()

#######


