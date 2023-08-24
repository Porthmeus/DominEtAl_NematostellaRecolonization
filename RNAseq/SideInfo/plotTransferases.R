# Porthmeus
# 12.06.20

require(DESeq2)
require(ggplot2)
require(reshape2)

load("../DESeq2/DeseqMat.RData")

cht <- c("NVE10195, NVE13506, NVE2689, NVE11295, NVE21891, NVE8010") # contig of the Chitin synthase
cht <- unlist(strsplit(cht, split = ", "))

cts <- as.data.frame(t(counts(deseq)[cht,]))
cts <- cbind(cts, Condition = gsub(".A$","",colData(deseq)$Condition))
cts[["Condition"]] <- factor(cts[["Condition"]], levels= c("GF","wt","bE","bJ","bA"))

cts <- melt(cts,id.vars="Condition")
p <- ggplot(cts, aes(y=value,x=Condition)) +
    geom_boxplot(outlier.shape = NULL)+
    geom_point(shape= 21, aes(fill=Condition)) +
    ylab("Normalized read counts") +
    facet_wrap(~variable, scale = "free")+
    theme_bw()
ggsave(file = "TransferasesCounts.pdf", height=3, width=6, p)

p2 <- ggplot(cts[cts[["Condition"]] != "wt",], aes(y=value,x=Condition)) +
    geom_boxplot(outlier.shape = NULL)+
    geom_point(shape= 21, aes(fill=Condition)) +
    ylab("Normalized read counts") +
    facet_wrap(~variable, scale = "free")+
    theme_bw()

ggsave(file = "TransferasesCounts2.pdf", height=3, width=6, p2)



