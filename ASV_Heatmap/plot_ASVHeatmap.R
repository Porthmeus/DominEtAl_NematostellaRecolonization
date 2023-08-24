## ----load_data, message = FALSE------------------------------------------
require(pheatmap)
require(data.table)
require(RColorBrewer)
require(ggplot2)
require(factoextra)
require(cowplot)

dataDir <- "data"

meta <- fread(file.path(dataDir,"mapping_recolonization.tsv"))
mat <- fread(file.path(dataDir,"table_rel.tsv"))

colnames(mat)[1] <- "OTU"
mat2  <- as.matrix(mat[,-1])
rownames(mat2) <- mat[,OTU]


## ----rawfilter_data------------------------------------------------------
mat2 <- mat2[,meta[,sampleID]]
mat2 <- mat2[apply(mat2,1,sum)>0,]


## ----reorder_data--------------------------------------------------------
# create a vector of the order, how the samples should appear in the data
cdt_order <- c("bL","bJ","bA",
               "bL_A_2d","bL_A_7d","bL_A_14d","bL_A_28d",
               "bJ_A_2d","bJ_A_7d","bJ_A_14d","bJ_A_28d",
               "bA_A_2d","bA_A_7d","bA_A_14d","bA_A_28d")

# make the meta column factors and level it according to the above order. Use
# the order as well to sort the samples in the matrix and the meta table itself.
# Rename the samples in the matrix so that they are more human readable.
meta[,sample := factor(sample, levels=cdt_order)]
mat2 <- mat2[,order(meta[,sample])]
setkey(meta,sample)
colnames(mat2) <- meta[,name]

# create a vector to introduce gaps in the columns of the heatmap and a
# different coloration
gaps <- c(15,35,55) # each condition was sequenced 5 times
colscale <- colorRampPalette(c("black","orange","red"))(255)


## ----plot_rawHeatmap, fig.width=8, fig.height=4, fig.path="plots/"-------
pheatmap(mat2,
         show_rownames = FALSE,
         cluster_cols=FALSE,
         cluster_rows = FALSE,
         gaps_col = gaps,
         color = colscale,
         border_color = NA
)


## ----plot_ECDF, echo= FALSE, warning=FALSE, fig.height=4, fig.width=8, fig.path = "plots/"----
# melt the matrix into a data.table and merge this with the meta data to have
# everything in one table
df <- melt(
        data.table(mat2, keep.rownames = TRUE),
        id.vars = "rn", 
        variable.name ="name" )
df <- merge(df, meta, by = "name")

# extract the minmal value for each recol_inoc
minDf <- df[identity != "inocula",  .(minAb = min(value), name), by =.(dpr,rn )]

# generate a plot of the data
ecdf_p <- ggplot(minDf, aes(minAb)) + 
    stat_ecdf(geom = "step") +
    xlab("Min abundance for diff. inocular for recol.") +
    ylab("Fraction of total data")
# zoom a little into the plot to better understand the data 
ecdfZoom_p <-    ecdf_p + ylim(0.99,1) + xlim(0,0.001)
plot_grid(ecdf_p, ecdfZoom_p, align = "h", labels = "auto")


## ----plot_filteredHeatmap, fig.height = 4, fig.width = 8, fig.path = "plots/"----
# filter the matrix for only those bacteria with at least some reads in one of
# the conditions
sel <- unique(minDf[minAb >0, rn])
pheatmap(mat2[sel,],
         show_rownames = FALSE,
         cluster_cols=FALSE,
         cluster_rows = T,
         gaps_col = gaps,
         color = colscale,
         border_color = NA
)


## ----plot_filteredScaledHeatmap, fig.width=8, fig.height=4,fig.path="plots/"----
# helper function
minMaxNorm<-function(x){
    # normalizes a vector of numeric values to min=0 and max=1
    # returns 0.5 if x is a single value
    if(length(x)==1){
        return(0.5)
    }else{
        res<-c()
        mi<-min(x)
        ma<-max(x)
        for(i in x){
            res<-c(res,(i-mi)/(ma-mi))
        }
    }
    
    res[is.na(res)] <- 0
    
    return(res)
    
}
# construct a vector to indicate the order for normalization of the OTUs
plotMeta <- meta[,treatment]
inocL <- nrow(meta[identity == "inoculate",])
plotMeta[1:inocL] <- meta[identity == "inoculate",identity]
meta[["plotMeta"]] <- plotMeta

# do the actual normalization
matMinMax <- mat2[sel,]
for(i in unique(meta[,plotMeta])){
    smpl <- meta[plotMeta == i,name]
    matMinMax[,smpl] <- t(apply(mat2[sel,smpl], 1, minMaxNorm))
}

# do the clustering on the inoculate samples only
clustered_samples <- meta[identity == "inoculate", name]

value_order <- hclust(
    
    get_dist(mat2[sel,clustered_samples],
             method = "spearman",
             stand = FALSE),
    
    method = "median")[["order"]]

# after manual resorting use this vector:
ASV_order <- fread(file.path("data","ASV_order.csv"), header = FALSE)[[1]]
ASV2Number <- fread(file.path("data","ASV2Number.csv"))
setkey(ASV2Number, "ASV_ID")
rownames(matMinMax) <- ASV2Number[rownames(matMinMax), ASV]
# create the plot
p_minMax <- pheatmap(matMinMax[ASV_order,],
                  show_rownames = TRUE,
                  cluster_cols=FALSE,
                  cluster_rows = FALSE,
                  gaps_col = gaps,
                  color = colscale,
                  border_color = NA)
ggsave(p_minMax, height=8, width = 8,
       file = file.path("..","figures_raw","Fig2_ASVHeatmap.pdf"))

