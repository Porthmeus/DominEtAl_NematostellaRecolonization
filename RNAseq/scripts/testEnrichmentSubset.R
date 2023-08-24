# Porthmeus
# 19.03.20

# the script calculates the enrichment of gene groupings by RNAseq results. It
# takes a grouping of genes (this may be GO or KEGG annotation, or any other
# grouping, like protein families etc.) and their expression between two
# biological states (e. g. differential gene expression between two treatments
# or conditions). The script is tailored to use DESeq2 results output to
# calculate the enrichement. The enrichment itself is calculated by two
# different methods: first is GSEA [https://doi.org/10.1073%2Fpnas.0506580102]
# and the second is a hypergeometric test (basically a bunch of fisher's exact
# tests, with p-value correction for multiple testing).
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(clusterProfiler)
require(DESeq2)
require(data.table)
require(cowplot)
require(ggridges)
require(ggplot2)

msg <- "No significant cluster found"
# read the RNAseq data to define the universe
mat <- load(snakemake@input[["deseq"]])
deseq <- get(mat)
DEGs <- read.table(snakemake@input[["results"]])[[1]]
# read the annotation
anno <- fread(snakemake@input[["cluster"]])[,1:2, with =FALSE]
colnames(anno) <- c("cluster","gene")
anno[, cluster := substr(cluster, 1, 35)]
annoL <- anno[,.(ClusterSize = length(gene)),by=cluster]

# extract the DEGs from that RNAseq for hypergeometrical testing
DEGs <- DEGs[DEGs %in% anno[[2]]]
universe <- rownames(deseq)
# calculate the hypergeometrical test
hyGeo <- enricher(gene = DEGs, universe = universe, TERM2GENE = anno)
hyGeoDF <- data.frame(hyGeo)
if(nrow(hyGeoDF)>0){
    hyGeoDF <- merge(hyGeoDF, annoL, by.x="ID",by.y="cluster")
    hyGeoDF <- data.table(hyGeoDF)
    hyGeoDF[,ID := factor(ID, levels= ID[order(Count/ClusterSize)])]
    # save the data
    write.csv(hyGeoDF, file = snakemake@output[["hyGeoCSV"]])
    save(hyGeo, file = snakemake@output[["hyGeoRData"]])
    
    # create some plots for evaluation
    # define fontsize in case of many clusters being calculated
    if(nrow(hyGeoDF) > 25){
        fontSize <- 25/nrow(hyGeoDF)*11
    } else {
        fontSize <- 11
    }
    # a dotplot
    pDot <- ggplot(hyGeoDF, aes(x=Count/ClusterSize, y=ID, col=p.adjust, size = Count)) +
        geom_point() +
        scale_color_gradient(high="blue", low="red") +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank())
        
    # a barplot
    pBar <- ggplot(hyGeoDF, aes(x = Count, y=ID, fill=p.adjust)) +
        geom_bar(stat="identity") +
        scale_fill_gradient(high="blue", low = "red") +
        theme(legend.position = "none",
              axis.text.y = element_text(size = fontSize))
    
    pGrid <- plot_grid(pBar,pDot, rel_widths = c(2,1))
} else {
    pGrid <- ggplot(data.frame(1,1))+ ggtitle(msg)
    cat("ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue,geneID,Count,len,ClusterSize", file = snakemake@output[["hyGeoCSV"]])
    cat("ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue,geneID,Count,len,ClusterSize", file = snakemake@output[["hyGeoRData"]])
}

ggsave(pGrid, filename = snakemake@output[["hyGeoPlot"]], width =8, height=4, device = "svg")
