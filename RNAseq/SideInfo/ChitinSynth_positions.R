# Porthmeus
# 26.03.20

# Script for Hanna, to get the position of the Chitin-Synthases among the differential genes.
synth<-"NVE14301"

path <-"../DESeq2/"
files <- list.files(path = path, pattern = "Condition.*.csv")
positions <- c()
for(f in files){
    df <- read.csv(file.path(path,f), row.names = 1)
    if(synth %in% rownames(df)){
        df <- df[order(abs(df$log2FoldChange), decreasing =TRUE),]
        positions <- c(positions,grep(synth, rownames(df)))
    }
}
out <- data.frame(Condition = gsub(".csv","",files), NVE14301=positions)
write.csv(file = "ChitinSynth_positions.csv",out)
    
