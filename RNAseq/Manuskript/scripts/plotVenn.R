# Porthmeus
# 05.02.21

# plot a Venn diagram for more insight of DEG expression

library(limma)
library(data.table)

# load data 
venn <- read.csv("../DESeq2/VennMatrix.csv", row.names=1)
sel <- colnames(venn)[c(2:4,9)]
clrs <- fread("Colors.csv")

venn_sel <- venn[,sel]
colnames(venn_sel) <- gsub("bE","bL",gsub("Condition_|_vs_|GF|\\.A","",sel))

setkey(clrs,"Condition")
clrs_sel <- clrs[colnames(venn_sel),]
clrs_sel[is.na(Color), Color := "black"]
pdf("plots/vennDiagram.pdf", width=7, height = 7)
vennDiagram(venn_sel,circle.col = clrs_sel[, Color], main = "DEGs compared to GF", bty="n")
dev.off()
