# Porthmeus
# 06.02.20

# Takes the all the results files and creates a KEGG-mapping file which can be used to display DE genes in the different KEGG pathways
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(data.table)

# load data and extract the KO numbers
file.create(snakemake@output[["kegg"]])
for(res in snakemake@input[["results"]]){
    name <- paste("# ", gsub("_"," ",gsub(".csv","",res)), "\n", sep ="")
    df <- fread(res)
    df <- df[!is.na(KO),]
    cat(file = snakemake@output[["kegg"]], name, append = TRUE)
    for(i in 1:nrow(df)){
            cat(file = snakemake@output[["kegg"]], paste(paste(df[i,.(V1,KO)], collapse = "\t"), "\n", sep = ""), append = TRUE)
    }
}
