# Porthmeus
# 05.02.21

# plot the chitin synthasis as two seperate box jitter plots

require(ggplot2)
require(DESeq2)
require(data.table)

dds <- load("../DESeq2/DeseqMat.RData")
deseq <- get(dds)

chitins <- c("NVE14301","NVE8515","NVE2200","NVE25290","NVE9076","NVE1769","NVE12930")
clrs <- fread("Colors.csv") 
setkey(clrs, "Condition")

chtn_RC <- counts(deseq[chitins,], normalized=TRUE)
colnames(chtn_RC) <-gsub("bE","bL", 
                         gsub("\\.A","",
                              colData(deseq)[,"Condition"]
                              )
                        )
chtn_pl <- data.table(reshape2::melt(chtn_RC, value.name="RC"))
colnames(chtn_pl) <- c("Gene","Condition","RC")

clrs_sel <- clrs[unique(as.character(chtn_pl[,Condition])),]
setkey(clrs_sel, "Color")

chtn_pl[,Condition := factor(Condition, levels = c("wt","GF","bL","bJ","bA"))]

# EDIT 17.10.22 - remove the wt samples and set ylim to 0
chtn_pl <- chtn_pl[Condition != "wt",]

boxJitter1 <- ggplot(chtn_pl[Gene == chitins[1],], aes(x=Condition, y = RC)) +
                     geom_boxplot() +
                     geom_point(aes(fill = Condition), shape =21)+
                     scale_fill_manual(breaks = clrs_sel[,Condition],
                                        values = clrs_sel[,Color]) +
                     theme_bw() +
                     ylim(0,max(chtn_pl[Gene == chitins[1],RC])) +
                     ylab("Median ratio normalized read counts") +
                     ggtitle(chitins[1]) + 
                     theme(legend.position = "None")
boxJitter1

boxJitter2 <- ggplot(chtn_pl[Gene == chitins[2],], aes(x=Condition, y = RC)) +
                     geom_boxplot() +
                     geom_point(aes(fill = Condition), shape =21)+
                     scale_fill_manual(breaks = clrs_sel[,Condition],
                                        values = clrs_sel[,Color]) +
                     theme_bw() +
                     ylim(0,max(chtn_pl[Gene == chitins[2],RC])) +
                     ylab("Median ratio normalized read counts") +
                     ggtitle(chitins[2])+
                     theme(legend.position = "None")
boxJitter2

ggsave(boxJitter1, filename = paste0("plots/boxJitter_",chitins[1],".pdf"), height = 3, width =3)
ggsave(boxJitter2, filename = paste0("plots/boxJitter_",chitins[2],".pdf"), height = 3, width =3)


pp <- ggplot(chtn_pl[Gene %in% chitins[3:length(chitins)],], aes(x=Condition, y = RC)) +
                         geom_boxplot() +
                         geom_point(aes(fill = Condition), shape =21)+
                         scale_fill_manual(breaks = clrs_sel[,Condition],
                                            values = clrs_sel[,Color]) +
                         facet_wrap(~Gene, scale = "free")+
                         theme_bw() +
                         ylab("Median ratio normalized read counts") +
                         ggtitle("Bona fide mucins") + 
                         theme(legend.position = "None")
                     pp

ggsave(pp, filename = "plots/boxJitter_Mucins.pdf",
       width =5,
       height =5)


