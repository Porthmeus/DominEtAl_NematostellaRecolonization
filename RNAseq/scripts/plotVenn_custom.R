# Porthmeus
# 09.06.20

# plot a Venn diagram for the DEGs between GF and bJ, bA, bL respectively

#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

require(DESeq2)
require(limma)
require(data.table)

# load table
venn <- read.csv(snakemake@input[["venn"]], row.names =1)
vennAnno <- fread(snakemake@input[["vennAnno"]])
# create the plot
sel <- grep("GF_vs", colnames(venn), value = TRUE)[-1]
print(sel)
svg(snakemake@output[["venn"]], width =4, height=4)
vennDiagram(venn[,sel], circle.col = c(2,4,7), cex=0.7, names = gsub("Condition_GF_vs_","",sel))
dev.off()

# save the genes from the middle
sel2 <- apply(venn[,sel], 1, sum)==length(sel)
cat(file = snakemake@output[["geneList"]], paste(c("# overlap",rownames(venn)[sel2],"\n"), collapse="\n"))

sel2 <- rownames(venn)[sel2]
outtbl <- vennAnno[Contig %in% sel2,]
outtbl[, vennLocation := "overlap"]


sel2 <- apply(venn[,sel], 1, sum) == 1
for(i in 1:length(sel)){
    cln <- sel[i]
    sel3 <- venn[,cln] == 1
    contigs <- rownames(venn)[sel2 & sel3]
    
    tbl <- vennAnno[Contig %in% contigs,]
    outtbl <- rbind(outtbl,tbl[,vennLocation :=cln])
    cat(file = snakemake@output[["geneList"]], paste(c(paste("#",cln), contigs,"\n"), collapse = "\n"), append = TRUE)
}
write.csv(outtbl, file = snakemake@output[["geneTable"]])
