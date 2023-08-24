# Porthmeus
# 09.06.20

require(DESeq2)
require(limma)
require(ggplot2)
require(Rtsne)
require(plotly)

load("../DESeq2/DeseqMat.RData")
meta <- colData(deseq)

# transform
rld <- rlog(deseq)

# remove duplicates
rld <- assay(rld)
rld <- rld[!duplicated(rld),]

# remove batch effects
rld_ba <- removeBatchEffect(rld, batch = meta$Batch, batch2=meta$Batch2)

# remove the wild types
rld_ba_AB <- rld_ba[,meta$Condition != "wt"]

# get 500 most variable genes
sel <- names(sort(decreasing = TRUE, apply(rld_ba_AB, 1,var))[1:500])

# caclulate eigenvalues
pca <- prcomp(t(rld_ba_AB[sel,]))
pcax <- as.data.frame(pca$x)
pcax <- cbind(pcax, meta[meta$Condition != "wt",])

# create a plot
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, color = Condition)) +
    geom_point(size=2) +
    labs(x = paste0("PC1 (", round(pca$sdev[1]/sum(pca$sdev), digits=3)*100, "%)"),
         y = paste0("PC2 (", round(pca$sdev[2]/sum(pca$sdev), digits=3)*100, "%)"),
         title = "PCA")

pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, color = Condition)) +
    geom_point(size=2) +
    labs(x = paste0("PC1 (", round(pca$sdev[1]/sum(pca$sdev), digits=3)*100, "%)"),
         y = paste0("PC3 (", round(pca$sdev[3]/sum(pca$sdev), digits=3)*100, "%)"),
         title = "PCA")

pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, color = Condition)) +
    geom_point(size=2) +
    labs(x = paste0("PC2 (", round(pca$sdev[2]/sum(pca$sdev), digits=3)*100, "%)"),
         y = paste0("PC3 (", round(pca$sdev[3]/sum(pca$sdev), digits=3)*100, "%)"),
         title = "PCA")


# create a tsne plot
tsne <- Rtsne(t(rld_ba_AB), perplexity =6)
tsney <- cbind(as.data.frame(tsne$Y),meta[meta$Condition != "wt",])
p_tsne <- ggplot(tsney, aes(x=V1, y=V2, color = Condition)) +
                     geom_point(size = 2) +
                     labs(x = "Dim1",
                          y= "Dim2",
                          title = "t-SNE")


p_pca <- cowplot::plot_grid(pc12,pc13,pc23,p_tsne, nrow=2)
ggsave(p_pca, filename = "PCA_tSNE_plots.svg", width = 9, height = 6)

# create 3D-PCA plot
fig <- plot_ly(pcax, x=~PC1, y=~PC2, z=~PC3, color = ~Condition,
               hoverinfo = "text",
               text = paste("Replicate:",pcax$Replicate, "\n", "Sample:", pcax$Sample))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = paste0("PC1 (",
                                                   round(pca$sdev[1]/sum(pca$sdev), digits=3)*100,
                                                   "%)")),
                       yaxis = list(title = paste0("PC2 (",
                                                   round(pca$sdev[2]/sum(pca$sdev), digits=3)*100,
                                                   "%)")),
                       zaxis = list(title = paste0("PC3 (",
                                                   round(pca$sdev[3]/sum(pca$sdev), digits=3)*100,
                                                   "%)")))
                       ) 
                                    
saveWidget(fig, "3D_PCAplot.html")
