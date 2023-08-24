# Porthmeus
# 12.06.20

require(DESeq2)
require(ggplot2)

load("../DESeq2/DeseqMat.RData")

cht <- "NVE14301" # contig of the Chitin synthase
cts <- data.frame(counts = counts(deseq)[cht,], Condition = colData(deseq)$Condition)
cts[["Condition"]] <- factor(cts[["Condition"]], levels= c("GF","wt","bE.A","bJ.A","bA.A"))
p <- ggplot(cts, aes(y=counts,x=Condition)) +
    geom_boxplot(outlier.shape = NULL)+
    geom_point(shape= 21, aes(fill=Condition)) +
    ylab("Normalized read counts") +
    theme_bw()
ggsave(file = "ChitinSynthaseCounts.svg", height=3, width=4, p)

p2 <- ggplot(cts[cts[["Condition"]] != "wt",], aes(y=counts,x=Condition)) +
    geom_boxplot(outlier.shape = NULL)+
    geom_point(shape= 21, aes(fill=Condition)) +
    ylab("Normalized read counts") +
    theme_bw()

ggsave(file = "ChitinSynthaseCounts2.svg", height=3, width=4, p2)



