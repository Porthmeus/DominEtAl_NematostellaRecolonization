# Porthmeus
# 24.10.19

# redirect messages to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


require(tximport)
require(readr)
require(DESeq2)
require(data.table)


# read the SampleDescription
desc <- fread(snakemake@input[["meta"]], stringsAsFactors = TRUE)
desc[, Batch := factor(Batch)]

# import file into tximport object
files <- file.path(snakemake@params[["files"]],paste(desc[["file_ID"]],"quant", sep = "_"), "quant.sf")
print(files)
txi <- tximport(files, type = "salmon", txOut = TRUE)

# create the count matrices
dds <- DESeqDataSetFromTximport(txi, desc, ~ Batch2 + Batch + Condition)
deseq <- DESeq(dds, betaPrior=T)
rownames(deseq) <- limma::strsplit2(rownames(deseq), split = "_")[,1]
save(file = snakemake@output[[1]], deseq)

