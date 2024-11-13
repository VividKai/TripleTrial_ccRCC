#########Set work dir
setwd("/rsrch5/scratch/genomic_med/kyu3/projects/02.Pavlos_Tri_Nivo_Ipil_Sitra_RCC/bin")

#########libraries
source("/rsrch5/home/genomic_med/kyu3/scAnaly/function.R")
library(Seurat)
library(harmony)
library(PCAtools)
library(dplyr)
library(data.table)
library(limma)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(Nebulosa)
library(cetcolor)
library(doBy)
library(smplot2)
library(ggSCvis)
library(ggbreak)
library(immcp)
library(presto)
library(reshape2)
library(viridis)
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(dplyr)
library(plyr)
library(scRNAtoolVis)
library(clusterProfiler)
library(ComplexHeatmap)



############################Figure 3############################
##############Figure 3a##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)

combined$ModCluster <- factor(combined$ModCluster, c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast"))

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$ModCluster))]
pdf("../Results/Overview2311/231221.overview.tme.pdf", width = 4.5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "ModCluster", cols = color_v)+ggtitle(""))
dev.off()

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])

write.table(combined@meta.data[, c("ID", "ModCluster", "UMAP1", "UMAP2")], "../Results/RawData/figure3a_c.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Figure 3b##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)

combined$ModCluster <- factor(combined$ModCluster, rev(c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast")))

cellmarker <- fread("../../../datasets/01.Markers/human.major.cluster.txt", header = T, data.table = F)
features = cellmarker$Gene

excludemarkers <- c("PTPRC", "KRT5", "KRT8", "CLDN4", "CLDN7", "EPCAM", 
                    "LAMP3", "CLEC9A", "IDO1", "CD1C", "CLEC10A", "LILRA4",
                    "CLEC4C", "FCGR3B", "SELE", "MKI67", "NCAM1", "FCGR3A", "ACTG2",
                    "CDK1", "MSLN", "CALB2")
features <- setdiff(features, excludemarkers)

Idents(combined) <- combined$ModCluster
pdf ("../Results/Overview2311/231221.overview.tme.markers.pdf", height = 3, width = 8)
print(DotPlot(combined, features = unique(features))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "white", mid = "#7C93C3", high = "#750E21", midpoint = 1.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

tmpres <- DotPlot(combined, features = unique(features))
write.table(tmpres$data, "../Results/RawData/figure3b.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Figure 3c##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)

combined$ModCluster <- factor(combined$ModCluster, c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast"))

combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

tempObj <- combined
Idents(tempObj) <- tempObj$ModCluster
coord = Embeddings(object = tempObj, reduction = "umap")
coord = coord[,c(1,2)]
colnames(coord) = c("UMAP_1", "UMAP_2")
coord = data.frame(ID = rownames(coord), coord)
meta = tempObj@meta.data
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')
print(table(meta$Class))
minNum <- min(table(meta$Class))
print(minNum)

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(meta$ModCluster))]

meta <- meta %>%
  group_by(Class) %>%
  slice_sample(n = minNum)

g <- ggplot(data = meta, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = F,bins = 5) +
  #geom_density_2d_filled()
  geom_point(color = 'white', size = .0005, alpha = 0.1) +
  scale_fill_viridis(option="A") +
  facet_wrap(~Class,ncol = 4) +
  theme_black()
pdf("../Results/Overview2311/231221.tme.umap_density.pdf", width = 12, height = 3)
print(g)
dev.off()

g <- ggplot(data = meta, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = ModCluster), size = 0.2) +
  scale_color_manual(values = color_v) +
  facet_wrap(~Class,ncol = 4) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.position = "none")
pdf("../Results/Overview2311/231221.tme.umap.pdf", width = 14, height = 3)
print(g)
dev.off()

##############Figure 3d##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)

combined$ModCluster <- factor(combined$ModCluster, c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast"))

combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

cellratio <- prop.table(table(combined$ModCluster, as.factor(as.character(combined$Class))), margin = 2)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

write.table(cellratio, "../Results/RawData/figure3d.txt", row.names = F, col.names = T, quote = F, sep = "\t")

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$ModCluster))]

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var2, y= Freq, fill = Var1), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

pdf("../Results/Overview2311/231221.tme.cellfraction.barplot.pdf", width = 4, height = 3.5)
p
dev.off()

cellfrac <- fread("../Results/Overview2311/231221.tme.cellfraction.txt", header = T, data.table = F, stringsAsFactors = F)
cellfrac$Var1 <- factor(cellfrac$Var1, c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast"))
cellfrac$Var2 <- factor(cellfrac$Var2, c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

cellfrac <- cellfrac %>% group_by(Var1) %>% mutate(rescale = scale(Freq))
cellfrac$Freq <- round(cellfrac$Freq, 2)

p1 <- ggplot(cellfrac, aes(x = Var2, y = Var1, fill = rescale))+
  geom_tile()+
  geom_text(aes(label = Freq)) +
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Overview2311/231221.tme.diff.frac.pdf", p1, width = 5, height = 4)

##############Figure 3e##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)
combined$ModCluster <- factor(combined$ModCluster, c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast"))

tmp.combined <- combined
tmp.combined <- subset(tmp.combined, subset = Collection.point %in% c("Baseline"))
cellratio <- prop.table(table(tmp.combined$ModCluster, as.character(tmp.combined$MiratiID)), margin = 2)
cellratio <- as.data.frame(cellratio)

saminfo <- fread("../datasets/sample.qc.doublets.100623.txt", header = T, data.table = F)
saminfo$MiratiID <- paste("RCC", saminfo$`Mirati ID`, sep = "")

cellratio <- merge(cellratio, saminfo, by.x = "Var2", by.y = "MiratiID")
cellratio$Var1 <- as.character(cellratio$Var1)

allcors <- c()
allpva <- c()
for (cell in unique(cellratio$Var1)){
  tmp.cellratio <- unique(cellratio[cellratio$Var1 == cell, ])
  x <- tmp.cellratio$Freq
  y <- tmp.cellratio$`ORR (%)`
  diam.lm <- lm(y~x, tmp.cellratio)
  allcors <- c(allcors, (summary(diam.lm)$coefficients[2, 1] / abs(summary(diam.lm)$coefficients[2, 1]))*summary(diam.lm)$r.squared)
  allpva <- c(allpva, summary(diam.lm)$coefficients[2, 4])
  
}


res <- data.frame(ID = unique(cellratio$Var1),
                  R2 = allcors,
                  Pva = allpva)
res <- res[order(res$R2), ]
res$ID <- factor(res$ID, as.character(sort(unique(combined$ModCluster))))

write.table(res, "../Results/RawData/figure3e.txt", row.names = F, col.names = T, quote = F, sep = "\t")

res$logp <- -log10(res$Pva)
allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$ModCluster))]

p <- ggplot(res, aes(x = ID, y = R2))+
  geom_hline(yintercept = 0, color = "black", size = 0.5)+
  geom_segment(aes(x = ID, xend = ID, y = 0, yend = R2), color = "gray", linetype = "dashed")+
  theme_classic2()+
  geom_point(aes(color = ID, size = logp))+
  geom_text_repel(aes(x = ID, y = R2, label = ID))+
  ylab("Response Index")+
  scale_x_discrete(limit = as.character(res$ID))+
  scale_color_manual(values = color_v)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Overview2311/231221.tme.cellfraction.lm.r2.pva.pdf", p, width = 5, height = 4)




