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

# read the RNAseq data
mat <- load(snakemake@input[["results"]])
results <- get(mat)

# read the annotation
anno <- fread(snakemake@input[["cluster"]])[,1:2, with =FALSE]
colnames(anno) <- c("cluster","gene")
anno[, cluster := substr(cluster, 1, 35)]
annoL <- anno[,.(ClusterSize = length(gene)),by=cluster]

# extract the DEGs from that RNAseq for hypergeometrical testing
DEGs <- rownames(results[(results$padj <= 0.05 & !is.na(results$padj)),])
DEGs <- DEGs[DEGs %in% anno[[2]]]
universe <- unique(rownames(results))
# calculate the hypergeometrical test
hyGeo <- enricher(gene = DEGs, universe = universe,TERM2GENE = anno)
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

# now calculate the same for the GSEA
FC <- results[["log2FoldChange"]]
names(FC) <- rownames(results)
FC <- sort(FC, decreasing =TRUE)

# I reduce fold changes slightly to unbind ties in duplicated values, otherwise, GSEA throws a warning message
while(any(duplicated(FC))){
    FC[duplicated(FC)] <- FC[duplicated(FC)]-1E-10
}

gsea <- GSEA(FC, TERM2GENE = anno)

# the data format is kind of annoying and needs a lot of reformation
gseaDF <- data.table(data.frame(gsea))
if(nrow(gseaDF)>0){
    # get the number of genes in each cluster
    gseaDF <- merge(gseaDF, annoL, by.x="ID", by.y="cluster") 
    # get the number of genes in the core enrichment
    gseaDF[, Count := sapply(gregexpr("/", core_enrichment),length)+1]
    gseaDF[, CountClusterRatio := Count/ClusterSize]
    gseaDF[,ID := factor(ID, levels = ID[order(CountClusterRatio)])]
    
    # create some plots:
    # adjust the font size depending on the number of clusters enriched
    if(nrow(gseaDF) > 25){
        fontSize <- 25/nrow(gseaDF)*11
    } else {
        fontSize <- 11
    }
    # a dotplot
    gseaDot <- ggplot(gseaDF, aes(x=Count/ClusterSize, y = ID, size = Count, color = p.adjust)) +
        geom_point() +
        scale_color_gradient(high="blue", low="red") +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank())
    
    # to create nice looking ridgeplot some data needs to be tangeled again...
    pf <- c()
    nve <- c()
    for(i in 1:nrow(gseaDF)){
        n <- strsplit(gseaDF[i,core_enrichment], split = "/")[[1]]
        pf <- c(pf, rep(as.character(gseaDF[i,ID]), length(n)))
        nve <- c(nve,n)
    }
    ridgeDF <- data.frame(ID = factor(pf, levels = levels(gseaDF[,ID])),
                          gene = nve,
                          FC = FC[nve])
    ridgeDF <- merge(ridgeDF, gseaDF[,.(ID,p.adjust)])
    
    # generate the ridgeplot
    gseaRidge <- ggplot(ridgeDF, aes(x = FC, y = ID, fill = p.adjust)) + 
        geom_density_ridges() +
        scale_fill_gradient(high="blue", low = "red") +
        theme(legend.position = "none",
              axis.text.y = element_text(size = fontSize))
    
    # combine plots and save
    pGSEA <- plot_grid(gseaRidge, gseaDot, rel_widths = c(2,1))
    # save the data of gsea
    write.csv(gseaDF, file = snakemake@output[["gseaCSV"]])
    save(gsea, file = snakemake@output[["gseaRData"]])
} else {
    pGSEA <- ggplot(data.frame(1,1)) +ggtitle(msg)
    cat( "ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue,geneID,Count,len,ClusterSize", file = snakemake@output[["gseaCSV"]])
    cat( "ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue,geneID,Count,len,ClusterSize", file = snakemake@output[["gseaRData"]])
}
ggsave(pGSEA, 
       filename = snakemake@output[["gseaPlot"]],
       width = 8,
       height = 4,
       device = "svg")
