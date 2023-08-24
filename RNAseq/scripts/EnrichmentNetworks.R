# Porthmeus
# 24.03.20

require(clusterProfiler)
require(ggplot2)
require(data.table)

# load data
D_T <- load("KEGG/Enrichment/dataCondition_bE.A_bJ.A_bA.A_vs_GF_KEGG_hyGeo.RData")
data <- get(D_T)

anno <- fread("Annotation/Cluster2KeggLong.csv")
anno[, shortName := substr(Name, 1, regexpr(",",Name)-1)]
annoS <- anno[, .(cluster, shortName)]
annoS <- annoS[!duplicated(cluster),]

