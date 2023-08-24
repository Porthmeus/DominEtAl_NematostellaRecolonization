# Porthmeus
# 09.06.20

require(DESeq2)
require(limma)
require(ggplot2)
require(Rtsne)
require(umap)

load("../DESeq2/DeseqMat.RData")
meta <- colData(deseq)
clrs <- read.csv("Colors.csv")

# transform
rld <- rlog(deseq)

# remove duplicates
rld <- assay(rld)
rld <- rld[!duplicated(rld),]

# remove batch effects
rld_ba <- removeBatchEffect(rld, batch = meta$Batch, batch2=meta$Batch2)
 meta[["Condition"]] <- gsub("bE","bL", gsub("\\.A","",meta[["Condition"]]))
# remove the wild types
rld_ba_AB <- rld_ba[,meta$Condition != "wt"]


# create a tsne plot
set.seed(5526)
tsne <- Rtsne(t(rld_ba_AB), perplexity =6)
tsney <- cbind(as.data.frame(tsne$Y),meta[meta$Condition != "wt",])
p_tsne <- ggplot(tsney, aes(x=V1, y=V2, color = Condition)) +
                     geom_point(size = 2) +
                     labs(x = "Dim1",
                          y= "Dim2",
                          title = "t-SNE") +
                     scale_color_manual(breaks = clrs[["Condition"]],
                                        values = clrs[["Color"]]) +
                     theme_bw()
p_tsne

ggsave(p_tsne, filename = "plots/tSNE.pdf", width = 4, height = 3)


#umap <- umap(t(rld))
#lyot <- cbind(umap[["layout"]], meta)
umap <- umap(t(rld_ba_AB))
lyot <- cbind(umap[["layout"]], meta[meta$Condition != "wt",])

set.seed(-1)
p_umap <- ggplot(as.data.frame(lyot), aes(x=V1, y=V2, color = Condition)) +
                     geom_point(size = 2) +
                     labs(x = "UMAP1",
                          y= "UMAP2",
                          title = "UMAP") +
                     scale_color_manual(breaks = clrs[["Condition"]],
                                        values = clrs[["Color"]]) +
                     theme_bw()
p_umap
ggsave(p_umap, filename = "plots/UMAP.pdf", width = 4, height = 3)
