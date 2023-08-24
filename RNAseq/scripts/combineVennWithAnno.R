# Porthmeus
# 10.06.20

# logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(data.table)
# read the tables and merge, save
anno <- fread(snakemake@input[["anno"]])
venn <- fread(snakemake@input[["venn"]])
colnames(venn)[1] <- "Contig"

all <-merge(venn, anno, by ="Contig")
write.csv(file = snakemake@output[["vennAnno"]], all, row.names=FALSE)