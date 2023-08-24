# Porthmeus
# 10.06.20

# aggregate annotation to contig

require(data.table)
require(limma)

# write a small aggregation function
agr <- function(x){
    # simply cat non-duplicated values to together by a ';'
    paste(x[!duplicated(x)], collapse = ";")
}

intPro <- fread("Nv_annotation2.csv")

# generate a vector for the contig numbers
intPro[,Contig := strsplit2(V1, split = "::")[,2]]
intPro <- intPro[,.(Contig,V4,V5,V6)]
colnames(intPro) <- c("Contig","DB","ID","Desc")
dt_intPro <- dcast(intPro, Contig~DB, value.var = c("ID","Desc"), fun = agr)

# KEGG
kegg <- fread("Cluster2Kegg.csv")
colnames(kegg) <- c("Contig","ID_KEGG","Desc_KEGG")

# GO
GO <- fread("contig2go.csv")
GO <- GO[-1,2:3]
colnames(GO) <- c("Contig","ID_GO")
dt_GO <- dcast(GO, Contig~"ID_GO", fun = agr, value.var="ID_GO")
colnames(dt_GO) <- c("Contig","ID_GO")

# TM and SP
TMSP <- fread("Nv_SummaryAnno.csv")
dt_TMSP <- dcast(TMSP, transcript~"id", value.var = c("SP","No.TM"), fun = agr)
colnames(dt_TMSP)[1] <- "Contig"

# merge the tables
all <- merge(dt_intPro, kegg, by="Contig", all =TRUE)
all <- merge(all, dt_GO, by = "Contig", all = TRUE)
all <- merge(all, dt_TMSP, by = "Contig", all = TRUE)

# save it to disk
write.csv(file = "Contig2AnnoAggregated.csv", all, row.names = FALSE)
