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


############################Figure 4############################
##############Figure 4a##############
combined <- readRDS("../3.Annotation/MinorCluster/01.epi/harmony.Merged.pca.major.minorcluster.rds")
saminfo <- fread("../datasets/sample.qc.doublets.100623.txt", header = T, data.table = F)
saminfo$MiratiID <- paste("RCC", saminfo$`Mirati ID`, sep = "")
cellinfo <- data.frame(ID = rownames(combined@meta.data),
                       HRR = combined@meta.data$GSEID)
cellinfo <- merge(cellinfo, saminfo, by.x = "HRR", by.y = "Sample ID")
rownames(cellinfo) <- cellinfo$ID
combined <- AddMetaData(combined, cellinfo)
combined$MiratiID <- paste("RCC", combined$`Mirati ID`, sep = "")

nums <- paste("C", c(0:(length(unique(combined$MinorCluster)) -1)), sep = "")
Idents(combined) <- combined$MinorCluster
majorcluster <- nums
names(majorcluster)<-levels(combined)
combined<-RenameIdents(combined,majorcluster)
combined$CCluster<- Idents(combined)

combined$clusters <- as.numeric(gsub("C", "", Idents(combined)))

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$CCluster))]

p1 <- DimPlot(combined, label = T) +
  ## scale_color_manual(values = hue_pal( h.start = 0, direction = 1, l = 80)(15)) +
  ## scale_color_manual(values = umapColor) +
  scale_color_manual(values = color_v) +
  theme_void() + theme(text = element_text(size = 10),
                       legend.position = "right")

allcolors <- c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Accent"))
color_v2=allcolors[1:length(unique(combined$Patient))]
cellcount1 <- as.data.frame(table(combined$CCluster, combined$MiratiID))
names(cellcount1) <- c("Cluster", "MiratiID", "Freq")
cellcount1$Cluster <- factor(cellcount1$Cluster, sort(unique(cellcount1$Cluster), decreasing = T))
p2 <- ggplot(cellcount1) +
  geom_bar(aes(x =Freq, y= Cluster, fill = MiratiID), stat = "sum", width = 0.7, linewidth = 0.25,colour = '#222222')+
  theme_classic() +
  labs(x='',y ='')+
  #coord_flip()+
  scale_fill_manual(values = color_v2)+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v3=allcolors[1:length(unique(combined$Biopsy_organ_type))]

cellcount1 <- as.data.frame(table(combined$CCluster, combined$Biopsy_organ_type))
names(cellcount1) <- c("Cluster", "Biopsy_organ_type", "Freq")
cellcount1$Cluster <- factor(cellcount1$Cluster, sort(unique(cellcount1$Cluster), decreasing = T))
p3 <- ggplot(cellcount1) +
  geom_bar(aes(x =Freq, y= Cluster, fill = Biopsy_organ_type), stat = "sum", width = 0.7, linewidth = 0.25,colour = '#222222')+
  theme_classic() +
  labs(x='',y ='')+
  #coord_flip()+
  scale_fill_manual(values = color_v3)+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3


pdf("../Results/Overview2311/231221.overview.epi.harmony.cluster.patient.biopsy.pdf", width = 12, height = 3)
p1|p2|p3
dev.off()


combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])

write.table(combined@meta.data[, c("ID", "CCluster", "Biopsy_organ_type", "MiratiID", "UMAP1", "UMAP2")], "../Results/RawData/figure4a_b.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Figure 4b##############
combined <- readRDS("../3.Annotation/MinorCluster/01.epi/harmony.Merged.pca.major.minorcluster.rds")
saminfo <- fread("../datasets/sample.qc.doublets.100623.txt", header = T, data.table = F)
saminfo$MiratiID <- paste("RCC", saminfo$`Mirati ID`, sep = "")
cellinfo <- data.frame(ID = rownames(combined@meta.data),
                       HRR = combined@meta.data$GSEID)
cellinfo <- merge(cellinfo, saminfo, by.x = "HRR", by.y = "Sample ID")
rownames(cellinfo) <- cellinfo$ID
combined <- AddMetaData(combined, cellinfo)
combined$MiratiID <- paste("RCC", combined$`Mirati ID`, sep = "")
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

nums <- paste("C", c(0:(length(unique(combined$MinorCluster)) -1)), sep = "")
Idents(combined) <- combined$MinorCluster
majorcluster <- nums
names(majorcluster)<-levels(combined)
combined<-RenameIdents(combined,majorcluster)
combined$CCluster<- Idents(combined)

tempObj <- combined
Idents(tempObj) <- tempObj$CCluster
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
color_v=allcolors[1:length(unique(meta$CCluster))]

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
pdf("../Results/Overview2311/231221.epi.umap_density.pdf", width = 12, height = 3)
print(g)
dev.off()

g <- ggplot(data = meta, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = CCluster), size = 0.2) +
  scale_color_manual(values = color_v) +
  facet_wrap(~Class,ncol = 4) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.position = "none")
pdf("../Results/Overview2311/231221.epi.umap.pdf", width = 14, height = 3)
print(g)
dev.off()



##############Figure 4c##############
combined <- readRDS("../3.Annotation/MinorCluster/01.epi/harmony.Merged.pca.major.minorcluster.rds")
saminfo <- fread("../datasets/sample.qc.doublets.100623.txt", header = T, data.table = F)
saminfo$MiratiID <- paste("RCC", saminfo$`Mirati ID`, sep = "")
cellinfo <- data.frame(ID = rownames(combined@meta.data),
                       HRR = combined@meta.data$GSEID)
cellinfo <- merge(cellinfo, saminfo, by.x = "HRR", by.y = "Sample ID")
rownames(cellinfo) <- cellinfo$ID
combined <- AddMetaData(combined, cellinfo)
combined$MiratiID <- paste("RCC", combined$`Mirati ID`, sep = "")
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))


nums <- paste("C", c(0:(length(unique(combined$MinorCluster)) -1)), sep = "")
Idents(combined) <- combined$MinorCluster
majorcluster <- nums
names(majorcluster)<-levels(combined)
combined<-RenameIdents(combined,majorcluster)
combined$CCluster<- Idents(combined)

cellratio <- prop.table(table(combined$CCluster, as.factor(as.character(combined$Class))), margin = 2)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

write.table(cellratio, "../Results/RawData/figure4c.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


cellfrac <- fread("../Results/Overview2311/231221.epi.cellfraction.txt", header = T, data.table = F, stringsAsFactors = F)
cellfrac$Var1 <- factor(cellfrac$Var1, rev(unique(cellfrac$Var1)))
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


ggsave("../Results/Overview2311/231221.epi.diff.frac.corre.pdf", p1, width = 5, height = 4)

##############Figure 4d##############
combined <- readRDS("../3.Annotation/MinorCluster/01.epi/harmony.Merged.pca.major.minorcluster.rds")
saminfo <- fread("../datasets/sample.qc.doublets.100623.txt", header = T, data.table = F)
saminfo$MiratiID <- paste("RCC", saminfo$`Mirati ID`, sep = "")
cellinfo <- data.frame(ID = rownames(combined@meta.data),
                       HRR = combined@meta.data$GSEID)
cellinfo <- merge(cellinfo, saminfo, by.x = "HRR", by.y = "Sample ID")
rownames(cellinfo) <- cellinfo$ID
combined <- AddMetaData(combined, cellinfo)
combined$MiratiID <- paste("RCC", combined$`Mirati ID`, sep = "")
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

nums <- paste("C", c(0:(length(unique(combined$MinorCluster)) -1)), sep = "")
Idents(combined) <- combined$MinorCluster
majorcluster <- nums
names(majorcluster)<-levels(combined)
combined<-RenameIdents(combined,majorcluster)
combined$CCluster<- Idents(combined)

tmp.combined <- combined
tmp.combined <- subset(tmp.combined, subset = Collection.point %in% c("Baseline"))
cellratio <- prop.table(table(tmp.combined$CCluster, as.character(tmp.combined$MiratiID)), margin = 2)
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
res$ID <- factor(res$ID, as.character(sort(unique(combined$CCluster))))

write.table(res, "../Results/RawData/figure4d.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

res$logp <- -log10(res$Pva)

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$CCluster))]

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

ggsave("../Results/Overview2311/231221.epi.cellfraction.lm.r2.pva.pdf", p, width = 5, height = 4)





##############Figure 4e##############
path <- fread("../Results/Overview2401/2401.c8.enriched.path.txt", header = F, stringsAsFactors = F, data.table = F)
names(path) <- c("ID", "BG", "Descrip", "Insize", "Ratio", "Pva", "FDR")

path$logp <- -log10(path$FDR)
path$ID <- factor(path$ID, rev(path$ID))

p <- ggplot(path) +
  geom_col(aes(x =logp, y= ID, fill = logp), width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  #scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0))

p

pdf("../Results/Overview2401/2401.epi.c8.path.pdf", width = 8, height = 3)
p
dev.off()

##############Figure 4f##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$CCluster))]

metapros <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cancer.mps.txt", header = T, stringsAsFactors = F, data.table = F)
metapros$ID <- "ID"
metapros <- melt(metapros, id.vars = "ID")

allscores <- split(metapros$value, metapros$variable)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)


selgenes <- names(allscores)[c(12:16)]

resdat <- combined@meta.data[, c("CCluster", selgenes)]
resdat <- melt(resdat, id.vars = "CCluster")
names(resdat) <- c("CCluster", "MPs", "Score")

resdat$CCluster <- factor(resdat$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))

p <- ggplot(resdat, aes(x = CCluster, y = Score, fill = CCluster))+
  geom_boxplot(outlier.size = 0)+
  scale_fill_manual(values = color_v) +
  theme_classic()
q <- p + facet_wrap(vars(MPs), scales = "free_y", ncol = 3)

ggsave(q, width = 12, height = 4, filename = "../Results/Overview2401/2401.epi.cluster.selmetaprograms.pdf")

write.table(resdat, "../Results/RawData/figure4f.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


##############Figure 4g##############
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

scenicLoomPath <- "../4.Analy/02.epi/scenic/scenic.sce_regulon_AUC.loom"

loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)

combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
cellClusters <- combined@meta.data
close_loom(loom)

selectedResolution <- "CCluster" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pdf("../Results/Overview2401/2401.epi.scenic.tfs.heatmap.ccluster.selected.pdf", width = 6, height = 2)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[c("TWIST1(+),", "ZEB1(+),"), ], name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=6))
dev.off()

write.table(regulonActivity_byCellType_Scaled, "../Results/RawData/figure4g.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Figure 4h##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
Idents(combined) <- combined$CCluster

combined <- subset(combined, subset = CCluster %in% c("C2", "C8"))

selmarkers <- c("B2M", "CD63", "TGFB1", "TGFBI", "SERPINE2", "H19", 
                "TIMP1", "CAV1", "HGF")

Idents(combined) <- combined$ModCluster
pdf ("../Results/Overview2401/2401.epi.markers.pdf", height = 2, width = 5)
print(DotPlot(combined, features = unique(selmarkers))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "#3652AD", mid = "white", high = "#FF004D", midpoint = 0.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

tmpres <- DotPlot(combined, features = unique(selmarkers))
write.table(tmpres$data, "../Results/RawData/figure4h.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Figure 4i##############
topmarkers <- fread("../4.Analy/02.epi/harmony.epi.cluster.markers.top50.txt", header = T, stringsAsFactors = F, data.table = F)

library(GSVA)
expdat <- fread("../../../datasets/03.TCGA_GDC/gene_exp/TCGA-KIRC.htseq_fpkm.gene.uniq.tsv", header = T, data.table = F)
rownames(expdat) <- expdat$gene
expdat <- expdat[, -1]

pro_up.gnas <- data.frame(gene_symbol = topmarkers[topmarkers$cluster == "C8", ]$gene)
pro_up.gnas$gs_name <- "C8_signature"
tfsgmt <- split(pro_up.gnas$gene_symbol, pro_up.gnas$gs_name)
res_tfs <- gsva(as.matrix(expdat), tfsgmt, kcdf="Gaussian",method = "ssgsea",parallel.sz=10)
res_tfs <- as.data.frame(t(res_tfs))
res_tfs$sample = rownames(res_tfs)
surv_info <- fread("../../../datasets/03.TCGA_GDC/survival/TCGA-KIRC.survival.tsv", header = T)
surv_info$NewID <- substr(surv_info$sample, 1, 15)
dss_surv <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/03.TCGA/survival/survival%2FKIRC_survival.txt", header = T)
surv_info <- merge(surv_info, dss_surv, by.x = "NewID", by.y = "sample")

surv_info <- merge(surv_info, res_tfs, by.x = "sample", by.y = "sample")
sum(is.na(surv_info$DSS))
surv_info <- surv_info[is.na(surv_info$DSS) == FALSE, ]
surv_info$OS <- surv_info$OS.y
surv_info$OS.time <- surv_info$OS.time.y

surv_info$Group <- ifelse(surv_info$C8_signature > median(surv_info$C8_signature), "1.High", "0.Low")

surv_info$DSS.time <- surv_info$DSS.time / 365
library(survival)
library(survminer)
survplot <- function(dat = mydata, type = "DSS", fit = fit, pval = pval){
  p <- ggsurvplot(fit,
                  linetype = 1,
                  #censor.shape=45,
                  data = dat,
                  size = 1, # change line size
                  palette = c("#6495EDFF", "#FF4500FF"),# custom color palettes
                  #conf.int = TRUE, # Add confidence interval
                  pval = paste('p = ', pval), # Add p-value
                  risk.table = T, # Add risk table
                  #tables.theme = theme_survminer(font.main = 10),
                  #risk.table.col = "strata",# Risk table color by groups
                  legend = "right",
                  #legend.labs = c("G1 (n = 7)", "G2 (n = 26)", "G3 (n = 24)"), # Change legend labels
                  risk.table.height = 0.25, # Useful to change when you have multiple groups
                  ggtheme = theme_classic2(), # Change ggplot2 theme
                  xlab = "Time (days)",
                  ylab = paste0("Probability of ", type))
  return(p)
}
tmp <- summary(coxph((Surv(DSS.time, DSS)) ~ Group, data = surv_info))
fit <- survfit(Surv(DSS.time, DSS) ~ Group, data = surv_info)
pfs <- survplot(surv_info, type = "DSS-KIRC", fit = fit, pval = tmp$logtest[3])
pdf("../Results/Overview2401/2401.epi.signature.C8.survival.pdf", width = 6.5, height = 5)
print(pfs)
dev.off()

write.table(surv_info, "../Results/RawData/figure4i.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

