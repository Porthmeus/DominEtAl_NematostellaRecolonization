# Porthmeus
# 24.10.2019

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

reconstructContrast <- function(string){
    # @string - a file name parsed by snakemake containing the information for contrasting
    string  <- gsub(".csv","",basename(string))
    cnd <- strsplit(string, split = "_")[[1]][1]
    contrast <- strsplit(string, split = "_vs_")[[1]]
    c1 <- strsplit(contrast[1], split = "_")[[1]][-1]
    c2 <- strsplit(contrast[2], split = "_")[[1]]
    ctrst <- list(paste0(cnd,c1), paste0(cnd,c2))
    return(ctrst)
}



require(DESeq2)
require(data.table)
# load deseq data and annotation
mat <- load(snakemake@input[["deseq"]])
deseq <- get(mat)
anno <- fread(snakemake@input[["anno"]])
anno[anno==""] <- NA

# create the results tables and save them
cntrst <- reconstructContrast(snakemake@output[["results"]])
print(cntrst)

# save the RData object
df <- results(deseq, contrast = cntrst)
save(file = snakemake@output[["resultsR"]], df)

# create an data frame with all significant results and add all known annotations
sel <- !is.na(df[["padj"]])
df <- df[sel,]
df <- as.data.frame(df)
df[["cluster"]] <- rownames(df)
df <- merge(df,anno[,.(cluster,KO,Description)], all.x =TRUE)
rownames(df) <- df[["cluster"]]
write.csv(df[df[["padj"]] <= 0.05,-1], file = snakemake@output[["results"]])

