# Porthmeus
# 05.02.21


blcklist <- read.csv("./blacklistKegg2.csv")
cluster2Kegg <- read.csv("./Cluster2KeggLong.csv")
sel <- cluster2Kegg[["Pathway"]] %in% blcklist[[1]]
cluster2Kegg <- cluster2Kegg[!sel,]
write.csv(file="./Cluster2KeggLongWhite.csv", cluster2Kegg, row.names=FALSE)
