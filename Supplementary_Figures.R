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


############################Supplementary codes############################
############################Figure S2############################
##############overview tme cells##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)

combined$ModCluster <- factor(combined$ModCluster, c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast"))

allcolors <- c(RColorBrewer::brewer.pal(8, "Set2"))
color_v=allcolors[1:length(unique(combined$Biopsy.organ.type))]
pdf("../Results/Overview2401/2401.overview.tme.biopsy.site.pdf", width = 5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "Biopsy.organ.type", cols = color_v)+ggtitle(""))
dev.off()

allcolors <- c(RColorBrewer::brewer.pal(8, "Set2"))
color_v=allcolors[1:length(unique(combined$Class))]
pdf("../Results/Overview2401/2401.overview.tme.class.pdf", width = 5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "Class", cols = color_v)+ggtitle(""))
dev.off()

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "Biopsy.organ.type", "Class", "UMAP1", "UMAP2")], "../Results/RawData/figureS2b_c.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############overview epithelial cells##############
combined <- readRDS("../3.Annotation/MinorCluster/01.epi/Merged.pca.major.minorcluster.score.rds")
saminfo <- fread("../datasets/sample.qc.doublets.100623.txt", header = T, data.table = F)
saminfo$MiratiID <- paste("RCC", saminfo$`Mirati ID`, sep = "")
cellinfo <- data.frame(ID = rownames(combined@meta.data),
                       HRR = combined@meta.data$GSEID)
cellinfo <- merge(cellinfo, saminfo, by.x = "HRR", by.y = "Sample ID")
rownames(cellinfo) <- cellinfo$ID
combined <- AddMetaData(combined, cellinfo)
combined$MiratiID <- paste("RCC", combined$`Mirati ID`, sep = "")

combined$MinorCluster.main <- strsplit2(combined$MinorCluster, "_")[, 1]
combined$MinorCluster.minor <- strsplit2(strsplit2(combined$MinorCluster, "_")[, 2], "-")[, 1]
combined <- subset(combined, subset = Doublets == "Singlets")

combined$PRCR <- ifelse(combined$`Objective response (PR or CR)` == 1, "PR/CR", "PD/SD")
combined$PRCR <- ifelse(combined$MiratiID == "RCC21", "PR/CR", combined$PRCR)

combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))



allcolors <- c(RColorBrewer::brewer.pal(8, "Set2"))
color_v=allcolors[1:length(unique(combined$Biopsy.organ.type))]
pdf("../Results/Overview2401/2401.overview.epi.biopsy.site.pdf", width = 5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "Biopsy.organ.type", cols = color_v)+ggtitle(""))
dev.off()

allcolors <- c(RColorBrewer::brewer.pal(8, "Set2"))
color_v=allcolors[1:length(unique(combined$Class))]
pdf("../Results/Overview2401/2401.overview.epi.class.pdf", width = 5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "Class", cols = color_v)+ggtitle(""))
dev.off()

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MiratiID",  "Biopsy.organ.type", "seurat_clusters", "Class", "UMAP1", "UMAP2")], "../Results/RawData/figureS2a_b_c_d.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############epithelial infercnv##############
combined <- readRDS("../3.Annotation/MinorCluster/01.epi/Merged.pca.major.minorcluster.score.rds")
saminfo <- fread("../datasets/sample.qc.doublets.100623.txt", header = T, data.table = F)
saminfo$MiratiID <- paste("RCC", saminfo$`Mirati ID`, sep = "")
cellinfo <- data.frame(ID = rownames(combined@meta.data),
                       HRR = combined@meta.data$GSEID)
cellinfo <- merge(cellinfo, saminfo, by.x = "HRR", by.y = "Sample ID")
rownames(cellinfo) <- cellinfo$ID
combined <- AddMetaData(combined, cellinfo)
combined$MiratiID <- paste("RCC", combined$`Mirati ID`, sep = "")

combined$MinorCluster.main <- strsplit2(combined$MinorCluster, "_")[, 1]
combined$MinorCluster.minor <- strsplit2(strsplit2(combined$MinorCluster, "_")[, 2], "-")[, 1]
combined <- subset(combined, subset = Doublets == "Singlets")

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$MinorCluster.minor))]
pdf("../Results/Overview2401/2401.overview.epi.malig.pdf", width = 4.5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "MinorCluster.minor", cols = color_v)+ggtitle(""))
dev.off()

#########For infercnv
combined<-readRDS("../3.Annotation/MinorCluster/01.epi/Merged.pca.major.minorcluster.score.rds")

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$seurat_clusters))]
pdf("../Results/Overview2401/2401.overview.epi.malig.seuratclusetr.pdf", width = 4.5, height = 3)
print(DimPlot(combined, reduction = "umap", label = TRUE, group.by = "seurat_clusters", cols = color_v)+ggtitle(""))
dev.off()

infercnv_obj = readRDS("../3.Annotation/MinorCluster/01.epi/19_HMM_pred.Bayes_NetHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$Stroma
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- test_loc$Epi

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("Stromal",length(normal_loc)),rep("Epithelial",length(test_loc)))
)
anno.df <- merge(anno.df, combined@meta.data, by.x = "CB", by.y = "ID", all.x = T)
head(anno.df)
rownames(anno.df) <- anno.df$CB
anno.df <- anno.df[order(anno.df$seurat_clusters), ]
anno.df$MiratiID <- anno.df$`Mirati ID`

anno.df$class <- anno.df$class.x
anno.df <- anno.df[!(anno.df$class == "Epithelial" & is.na(anno.df$seurat_clusters) == T), ]

gn <- rownames(expr)
geneFile <- read.table("../3.Annotation/MinorCluster/01.epi/gencode.v25.annotation4infercnv.unique.bed",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]

#sort data
anno.df <- anno.df[order(anno.df$class, anno.df$MiratiID, decreasing=F), ]
expr <- expr[, rownames(anno.df)]

#annotation heatmap and color
allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
library(ComplexHeatmap)
library(circlize)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 0.8)))
color_v1=allcolors[1:length(na.omit(unique(anno.df$MiratiID)))]
color_v2=allcolors[1:length(na.omit(unique(anno.df$seurat_clusters)))]
color_v3=allcolors[1:length(unique(anno.df$class))]
color_v4=allcolors[1:length(na.omit(unique(anno.df$Biopsy.organ.type)))]

names(color_v1)=as.character(na.omit(unique(anno.df$MiratiID)))
names(color_v2)=as.character(na.omit(unique(anno.df$seurat_clusters)))
names(color_v3)=as.character(unique(anno.df$class))
names(color_v4)=as.character(na.omit(unique(anno.df$Biopsy.organ.type)))

left_anno <- rowAnnotation(df = anno.df[, c("MiratiID", "seurat_clusters", "class", "Biopsy.organ.type")],
                           col=list(MiratiID=color_v1,
                                    seurat_clusters=color_v2,
                                    class = color_v3, 
                                    Biopsy.organ.type = color_v4),
                           na_col = "white")


#heatmap
ranselect <- sort(sample(1:ncol(expr), 1000))
pdf("../Results/Overview2401/2401.overview.epi.infercnv.res.cluster.heatmap.pdf",width = 10,height = 6)
ht = Heatmap(t(expr)[rownames(anno.df),][ranselect, ], #绘图数据的CB顺序和注释CB顺序保持一致
             col = colorRamp2(c(0,3,6), c("#377EB8","#F0F0F0","#E41A1C")), #10x(0.8, 1, 1.2)   smartseq(0.4, 1, 1.6)
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #sort chr
             column_gap = unit(2, "mm"),
             row_split = anno.df[ranselect, ]$class,
             top_annotation = top_anno,left_annotation = left_anno[ranselect, ], #add annotation
             row_title = NULL,column_title = NULL,
             use_raster = T, raster_device = "png")
draw(ht, heatmap_legend_side = "right")
dev.off()

##############tme cell fraction diff for matched samples##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC24", "RCC27"))
combined <- subset(combined, subset = Collection.point %in% c("Baseline", "C2"))
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)
combined$ModCluster <- factor(combined$ModCluster, c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast"))

cellratio <- prop.table(table(combined$ModCluster, as.character(combined$GSEID), as.character(combined$Collection.point)), margin = 2)
cellratio <- as.data.frame(cellratio)
cellratio <- cellratio[cellratio$Freq != 0, ]
cellratio <- na.omit(cellratio)

cellratio$Var3 <- factor(cellratio$Var3, c("Baseline", "C2"))
cellratio$Var1 <- factor(cellratio$Var1, c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast"))

write.table(cellratio, "../Results/RawData/figureS2g.txt",
            row.names = T, col.names = NA, sep = "\t", quote = F)


p <- ggplot(cellratio, aes(x = Var3, y = Freq, fill = Var3))+
  geom_boxplot()+
  geom_point()+
  #geom_label_repel(data = cellratio, aes(x = Var3, y = Freq, label = Var2))+
  stat_compare_means()+
  theme_classic()+
  ylab("% among all TME cells")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
q <- p + facet_wrap(vars(Var1), scales = "free_y", ncol = 7)

ggsave("../Results/Overview2401/2401.overview.tme.cellfraction.boxplot.pdf", q, width = 10, height = 3.5)
ggsave("../Results/Overview2401/2401.overview.tme.cellfraction.boxplot.big.pdf", q, width = 20, height = 6)

############################Figure S3############################
##############harmony estimate##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")

allcolors <- c(RColorBrewer::brewer.pal(8, "Set2"))
color_v=allcolors[1:length(unique(combined$Biopsy.organ.type))]
pdf("../Results/Overview2401/2401.overview.harmony.epi.biopsy.site.pdf", width = 5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "Biopsy.organ.type", cols = color_v)+ggtitle(""))
dev.off()

allcolors <- c(RColorBrewer::brewer.pal(8, "Set2"))
color_v=allcolors[1:length(unique(combined$MiratiID))]
pdf("../Results/Overview2401/2401.overview.harmony.epi.class.pdf", width = 4, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "MiratiID", cols = color_v)+ggtitle(""))
dev.off()

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MiratiID",  "Biopsy.organ.type", "Class", "UMAP1", "UMAP2")], "../Results/RawData/figureS3a.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############EMT level##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")

library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(dplyr)
library(plyr)
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol)
hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
hall_kegg <- rbind(kegg, hall)
hall_kegg <- hall_kegg[gsub(".*_", "", hall_kegg$gs_name) != "DISEASE", ]

allscores <- split(hall_kegg$gene_symbol, hall_kegg$gs_name)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

selpath <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

selhall_kegg <- hall_kegg[hall_kegg$gs_name %in% selpath, ]

p1 <- FeaturePlot(combined, reduction = "umap",
                  features = selpath, order = TRUE, ncol = 5) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()
ggsave(p1, width = 20, height = 10, filename = "../Results/Overview2401/2401.epi.modulescore.featureplot.pdf")


##############MPs##############
markers <- fread("../4.Analy/02.epi/harmony.epi.cluster.markers.txt", header = T, stringsAsFactors = F, data.table = F)

metapros <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cancer.mps.txt", header = T, stringsAsFactors = F, data.table = F)
metapros$ID <- "ID"
metapros <- melt(metapros, id.vars = "ID")
metapros <- metapros[, -1]


toppath <- c()
mergedpath <- data.frame()
for (celltype in unique(markers$cluster)) {
  pro_up <- enricher(markers[markers$cluster == celltype, ]$gene, TERM2GENE = metapros, pvalueCutoff = 1,minGSSize = 1)
  pro_up@result$Class <- celltype
  toppath <- unique(c(toppath, pro_up@result$ID[1:10]))
  mergedpath <- rbind(mergedpath, pro_up@result)
}

res <- mergedpath[mergedpath$ID %in% toppath, ]

write.table(res, "../Results/RawData/figureS3b.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

res$logp <- -log10(res$p.adjust)

res4sig <- res[res$p.adjust < 0.01, ]
res <- res[res$ID %in% res4sig$ID, ]
res$Class <- factor(res$Class, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))

res <- res[order(res$p.adjust, decreasing = T), ]
res <- res[order(res$Class, decreasing = F), ]
res$ID <- factor(res$ID, unique(res$ID))

p <- ggplot(res, aes(x = Class, y = ID, fill = logp))+
  geom_tile()+
  geom_text(data = res4sig, aes(x = Class, y = ID), label = "*")+
  theme_classic()+
  scale_fill_gradient2(low = "white", high = "red")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Overview2401/2401.epi.cluster.metaprogram.pdf", p, width = 8, height = 6)


##############regulons##############
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
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

cellClusters <- combined@meta.data
close_loom(loom)

selectedResolution <- "Class" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
# plot:
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
pdf("../Results/Overview2401/2401.epi.scenic.tfs.heatmap.pdf")
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=6))
dev.off()

regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC), selectedResolution])
## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

pdf("../Results/Overview2401/2401.epi.scenic.tfs.showrss.pdf")
print(plotly::ggplotly(rssPlot$plot))
dev.off()

common_tfs <- intersect(topRegulators[topRegulators$CellType == "Baseline_PD/SD", ]$Regulon, 
                        topRegulators[topRegulators$CellType == "EOT", ]$Regulon)

selectedResolution <- "CCluster" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pdf("../Results/Overview2401/2401.epi.scenic.tfs.heatmap.ccluster.selected.pdf", width = 6, height = 8)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[common_tfs, ], name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=6))
dev.off()

write.table(regulonActivity_byCellType_Scaled, "../Results/RawData/figureS3c.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############survival##############
surv_info <- fread("../Results/Revised07/2407.epi.signature.C8.survival.tcga.surv.txt", header = T, stringsAsFactors = F, data.table = F)
phenotype <- fread("../../../datasets/03.TCGA_GDC/phenotype/TCGA-KIRC.GDC_phenotype.tsv", header = T, stringsAsFactors = F, data.table = F)
phenotype <- phenotype[, c("submitter_id.samples", "history_of_neoadjuvant_treatment", "tumor_stage.diagnoses")]

surv_info <- merge(surv_info, phenotype, by.x = "sample", by.y = "submitter_id.samples")
tmp <- summary(coxph((Surv(DSS.time, DSS)) ~ Group + history_of_neoadjuvant_treatment + tumor_stage.diagnoses, data = surv_info))
fit <- survfit(Surv(DSS.time, DSS) ~ Group, data = surv_info)
pval <- tmp$coefficients[1, "Pr(>|z|)"]
pfs <- survplot(surv_info, type = "DSS-KIRC", fit = fit, pval = pval)
pdf("../Results/Revised07/2407.epi.signature.C8.survival.tcga.km.multi.pdf", width = 6.5, height = 5)
print(pfs)
dev.off()

write.table(surv_info, "../Results/RawData/figureS3g_i.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


##############figureS3e_f##############

combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
topmarkers <- fread("../4.Analy/02.epi/harmony.epi.cluster.markers.top50.txt", header = T, stringsAsFactors = F, data.table = F)

combined <- AddModuleScore(combined,
                           features = list(topmarkers[topmarkers$cluster == "C8", ]$gene),
                           name="C8_signatures")

tmp.combined <- subset(combined, subset = Collection.point == "Baseline")
tmp.combined$ORR <- tmp.combined$`ORR (%)`

res <- aggregate(.~ MiratiID, tmp.combined@meta.data[, c("C8_signatures1", "ORR", "MiratiID")], median)
p <- ggplot(data = res, mapping = aes(x = C8_signatures1, y = ORR)) +
  geom_point(shape = 21, fill = "#0f993d", color = 'black', size = 3) +
  sm_statCorr()+
  ggtitle("C8_signatures")
ggsave("../Results/Overview2401/2401.epi.signature.orr.C8_signatures.pdf", p, width = 5, height = 4)

write.table(res, "../Results/RawData/figureS3e.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


medscore <- median(combined$C8_signatures1)
combined$C8_group <- ifelse(combined$C8_signatures1 > medscore, "C8.high", "C8.low")
pdf("../Results/Overview2401/2401.epi.signature.C8.group.dimplot.pdf", width = 5.5, height = 4)
DimPlot(combined, reduction = "umap", pt.size = 2.0, label = FALSE, 
        raster=TRUE, group.by = "C8_group")
dev.off()

tmp.combined <- subset(combined, subset = Collection.point == "Baseline")
res <- aggregate(.~ MiratiID, tmp.combined@meta.data[, c("C8_signatures1", "PFS time (months)", "PFS event", "MiratiID")], median)
names(res) <- c("MiratiID", "C8_signatures1", "PFS.time", "PFS.status")
res$PFS.status <- ifelse(res$PFS.status == 1, "P/D", "Censored")
corre <- cor.test(res$C8_signatures1,res$PFS.time,method="pearson")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p1 <- ggplot(res, aes(x = C8_signatures1, y = PFS.time, color = PFS.status))+
  geom_point()+
  ggtitle(plottitle)+
  geom_smooth(method="lm",color="#1a9641") + 
  xlab("C8_signatures")+
  ylab("PFS")+
  theme_classic()
p1
ggsave("../Results/Overview2401/2401.epi.signature.C8.pfs.pdf", p1, width = 5.5, height = 4)


write.table(res, "../Results/RawData/figureS3f.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

############################Figure S4_S5############################
##############marker expression##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")

markers <- fread("../4.Analy/03.t/CD4T.cluster.markers.txt", header = T, stringsAsFactors = F, data.table = F)
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)

skipgenes <- c(markers$gene[grep("\\.", markers$gene)],
               markers$gene[grep("\\-", markers$gene)],
               markers$gene[grep("^RP", markers$gene)],
               markers$gene[grep("^LINC", markers$gene)])

combined$MinorCluster <- factor(combined$MinorCluster, sort(unique(combined$MinorCluster)))
Idents(combined) <- combined$MinorCluster
pdf ("../Results/Overview2401/2401.tcell.cd4t.marker.exp.top10.pdf", height = 4, width = 18)
print(DotPlot(combined, features = setdiff(unique(top10$gene), skipgenes))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "white", mid = "#7C93C3", high = "#750E21", midpoint = 1.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

tmpres <- DotPlot(combined, features = setdiff(unique(top10$gene), skipgenes))
write.table(tmpres$data, "../Results/RawData/figureS4a_b.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")

markers <- fread("../4.Analy/03.t/CD8T.cluster.markers.txt", header = T, stringsAsFactors = F, data.table = F)
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)

skipgenes <- c(markers$gene[grep("\\.", markers$gene)],
               markers$gene[grep("\\-", markers$gene)],
               markers$gene[grep("^RP", markers$gene)],
               markers$gene[grep("^LINC", markers$gene)])

combined$MinorCluster <- factor(combined$MinorCluster, sort(unique(combined$MinorCluster)))
Idents(combined) <- combined$MinorCluster
pdf ("../Results/Overview2401/2401.tcell.cd8t.marker.exp.top10.pdf", height = 4, width = 15)
print(DotPlot(combined, features = setdiff(unique(top10$gene), skipgenes))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "white", mid = "#7C93C3", high = "#750E21", midpoint = 1.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

tmpres <- DotPlot(combined, features = setdiff(unique(top10$gene), skipgenes))
write.table(tmpres$data, "../Results/RawData/figureS5a_b.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############signature##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")

allscores <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd4t.signature.txt", header = T)

cd4yanshuo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/CD4TcellMarkers.Yanshuo.top50.txt", header = T, stringsAsFactors = F, data.table = F)
selmarkers1 <- cd4yanshuo[cd4yanshuo$cluster == "CD4_c1_Treg", ][c(1:30), c(1, 2)]
names(selmarkers1) <- c("Signature", "Genes")
allscores <- rbind(allscores, selmarkers1)

selmarkers1 <- cd4yanshuo[cd4yanshuo$cluster == "CD4_c4_Tstr", ][c(1:50), c(1, 2)]
names(selmarkers1) <- c("Signature", "Genes")
allscores <- rbind(allscores, selmarkers1)

allscores <- split(allscores$Genes, allscores$Signature)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

combined$CD4_c4_Tstr <- ifelse(combined$CD4_c4_Tstr < -0.2, -0.2, combined$CD4_c4_Tstr)
combined$`TCR signaling` <- ifelse(combined$`TCR signaling` < -0.15, -0.15, combined$`TCR signaling`)

p1 <- FeaturePlot(combined, reduction = "umap",
                  features = c("Na\303\257ve", "TCR signaling", "CD4_c4_Tstr", "CD4_c1_Treg", "Cytotoxicity"), 
                  order = TRUE, ncol = 5) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()
ggsave(p1, width = 25, height = 4, filename = "../Results/Overview2401/2401.tcell.cd4t.overview.signature.pdf")

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MinorCluster", "Class", "UMAP1", "UMAP2", "Na\303\257ve", "TCR signaling", "CD4_c4_Tstr", "CD4_c1_Treg", "Cytotoxicity")], "../Results/RawData/figureS4c_d.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)



combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")

allscores <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd8t.signature.txt", header = T)

cd4yanshuo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/CD8TcellMarkers.Yanshuo.top50.txt", header = T, stringsAsFactors = F, data.table = F)
selmarkers1 <- cd4yanshuo[cd4yanshuo$icluster == "CD8_c1_Tex", ][c(1:30), c(1, 2)]
names(selmarkers1) <- c("Signature", "Genes")
allscores <- rbind(allscores, selmarkers1)

selmarkers1 <- cd4yanshuo[cd4yanshuo$icluster == "CD8_c4_Tstr", ][c(1:50), c(1, 2)]
names(selmarkers1) <- c("Signature", "Genes")
allscores <- rbind(allscores, selmarkers1)


allscores <- split(allscores$Genes, allscores$Signature)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

p1 <- FeaturePlot(combined, reduction = "umap",
                  features = c("Na\303\257ve", "TCR Signaling", "CD8_c4_Tstr", "CD8_c1_Tex", "Cytotoxicity"), 
                  order = TRUE, ncol = 5) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()
ggsave(p1, width = 25, height = 4, filename = "../Results/Overview2401/2401.tcell.cd8t.overview.signature.pdf")

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MinorCluster", "Class", "UMAP1", "UMAP2", "Na\303\257ve", "TCR Signaling", "CD8_c4_Tstr", "CD8_c1_Tex", "Cytotoxicity")], "../Results/RawData/figureS5c_d.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

############################Figure S6############################
##############signature marker expression##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")

allscores <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.txt", header = F)
names(allscores) <- c("gs_name", "gene_name")
immsets <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/immune.sets.txt", header = F)
names(immsets) <- c("gs_name", "gene_name")

makgs <- rbind(allscores, immsets)
rownames(makgs) <- makgs$gene_name

Idents(combined) <- combined$MinorCluster
Idents(combined) <- factor(Idents(combined), sort(unique(combined$MinorCluster)))

pdat <- DotPlot(combined, features = unique(makgs$gene_name))
data <- pdat$data

data <- cbind(data, makgs)
data$gs_name <- factor(data$gs_name, unique(data$gs_name))

p <- ggplot(data, aes(x = gene_name, y = id))+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  geom_point(data = data, aes(size=pct.exp, fill = avg.exp.scaled), shape = 21, colour="black", stroke=0.5) +
  scale_fill_gradient2(low = "white", mid = "#7C93C3", high = "#750E21", midpoint = 1.0) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(. ~ gs_name, scales = "free_x", space = "free")
ggsave("../Results/Overview2401/2401.myeloid.overview.markers.signature.pdf", p, height = 5, width = 15)

write.table(data, "../Results/RawData/figureS6a_c.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############signature feature plot##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")
allscores <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.txt", header = F)
allscores <- split(allscores$V2, allscores$V1)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$MinorCluster))]

#c("M1", "M2","Angiogenesis", "Phagocytosis")

p1 <- FeaturePlot(combined, reduction = "umap",
                  features = c("M1", "M2","Angiogenesis", "Phagocytosis"), 
                  order = TRUE, ncol = 4) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()
ggsave(p1, width = 20, height = 4, filename = "../Results/Overview2401/2401.myeloid.overview.signature.pdf")

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MinorCluster", "Class", "UMAP1", "UMAP2", "M1", "M2","Angiogenesis", "Phagocytosis")], "../Results/RawData/figureS6d.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############M1 and M2 signature##############
markers <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/new.myeloid.signature.txt", header = T, stringsAsFactors = F, data.table = F)
write_gmt(markers, "/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.gmt")

combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.rds")
combined$Class <- as.character(combined$Class)
##NR vs R
sub_sce <- subset(combined, subset = Class %in% c("Baseline_PR/CR", "Baseline_PD/SD"))

sub_sce.genes <- wilcoxauc(sub_sce, 'Class')
sub_sce.genes <- sub_sce.genes[sub_sce.genes$group == "Baseline_PD/SD", ]
sub_sce.genes <- sub_sce.genes[order(sub_sce.genes$logFC, decreasing = T), ]

write.table(sub_sce.genes[, c("feature", "logFC")], "../Results/Overview2401/2401.myeloid.overview.m1m2.seled.rnk",
            row.names = F, col.names = F, sep = "\t", quote = F)

#gseapy prerank -r ../Results/Overview2401/2401.myeloid.overview.m1m2.seled.rnk -g /rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.gmt --min-size 5 -o ../Results/Overview2401/gsea_nr_r
#gseapy prerank -r ../Results/Overview2401/2401.myeloid.overview.m1m2.seled.rnk -g /rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/h.all.v2023.2.Hs.symbols.gmt --min-size 5 -o ../Results/Overview2401/gsea_nr_r

##EOT vs Baseline
sub_sce <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC24"))
sub_sce <- subset(sub_sce, subset = Class %in% c("Baseline_PR/CR", "Baseline_PD/SD", "EOT"))
sub_sce$Class <- ifelse(sub_sce$Class %in% c("Baseline_PR/CR", "Baseline_PD/SD"), "Baseline", "EOT")
sub_sce.genes <- wilcoxauc(sub_sce, 'Class')
sub_sce.genes <- sub_sce.genes[sub_sce.genes$group == "EOT", ]
sub_sce.genes <- sub_sce.genes[order(sub_sce.genes$logFC, decreasing = T), ]

skipgenes <- c(sub_sce.genes$feature[grep("\\.", sub_sce.genes$feature)],
               sub_sce.genes$feature[grep("\\-", sub_sce.genes$feature)],
               sub_sce.genes$feature[grep("^RP", sub_sce.genes$feature)],
               sub_sce.genes$feature[grep("^LINC", sub_sce.genes$feature)])

sub_sce.genes <- sub_sce.genes[!(sub_sce.genes$feature %in% skipgenes), ]

write.table(sub_sce.genes[, c("feature", "logFC")], "../Results/Overview2401/2401.myeloid.overview.m1m2.seled.rnk",
            row.names = F, col.names = F, sep = "\t", quote = F)

#gseapy prerank -r ../Results/Overview2401/2401.myeloid.overview.m1m2.seled.rnk -g /rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.gmt --min-size 5 -o ../Results/Overview2401/gsea_eot_base
#gseapy prerank -r ../Results/Overview2401/2401.myeloid.overview.m1m2.seled.rnk -g /rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/h.all.v2023.2.Hs.symbols.gmt --min-size 5 -o ../Results/Overview2401/gsea_eot_base

##############Some signatures##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

markers <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/macrophage.m1.m2.txt", header = T, stringsAsFactors = F, data.table = F)
markers <- markers[markers$Signature == "M2", ]
Idents(combined) <- combined$Class
selgenes <- c("CTSA", "CTSB", "CTSC", "CTSD", "MSR1", "APOE", "C1QA", "C1QB", "CSF1R")
pdf ("../Results/Revised07/2407.myeloid.marker.sel4class.pdf", height = 2, width = 6.5)
print(DotPlot(combined, features = selgenes)+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "#3652AD", mid = "white", high = "#FF004D", midpoint = 0.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

tmpres <- DotPlot(combined, features = selgenes)
write.table(tmpres$data, "../Results/RawData/figureS6e.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Some signatures##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

Idents(combined) <- combined$Class
markers <- FindAllMarkers(combined, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, "../Results/Overview2407/2407.myeloid.class.DEGs.tsv", row.names = F, col.names = T, sep = "\t", quote = F)

markers <- fread("../Results/Overview2407/2407.myeloid.class.DEGs.tsv", header = T, stringsAsFactors = F, data.table = F)

color_v <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072")
pdf("../Results/Overview2407/2407.myeloid.class.DEGs.pdf", width = 5, height = 5)
jjVolcano(diffData = markers,
          tile.col = color_v,
          size  = 3.5,
          fontface = 'italic',
          cluster.order = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT")
)
dev.off()


############################Figure S7############################
##############B cell umap##############
combined <- readRDS("../4.Analy/05.b/B.reclassified.reannot.cells.rds")
allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$MinorCluster))]
pdf("../Results/Overview2401/2401.bcell.overview.pdf", width = 5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "MinorCluster", cols = color_v)+ggtitle("T cells"))
dev.off()

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MinorCluster", "UMAP1", "UMAP2")], "../Results/RawData/figureS7a.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############B cell fraction##############
combined <- readRDS("../4.Analy/05.b/B.reclassified.reannot.cells.rds")
#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Overview2401/2401.bcell."

OR.all.list <- do.tissueDist(cellInfo.tb=meta.tb,
                             meta.cluster = meta.tb$Class,
                             colname.patient = "MiratiID",
                             loc = meta.tb$MinorCluster,
                             out.prefix=sprintf("%s.Class",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../Results/Overview2401/2401.bcell.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Overview2401/2401.bcell.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$variable <- as.character(cellfrac$variable)
cellfrac$variable <- factor(cellfrac$variable, c("Bm_C0", "Bn_C1", "B_C2_ISG15", "Plasma_C3_JCHAIN", "Plasma_C4_IGHA1"))

cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Overview2401/2401.bcell.Class.diff.pdf", p1, width = 5, height = 2.5)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(cellratio$Var2))]
cellratio$Var1 <- factor(cellratio$Var1, c("Bm_C0", "Bn_C1", "B_C2_ISG15", "Plasma_C3_JCHAIN", "Plasma_C4_IGHA1"))

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Overview2401/2401.bcell.Class.diff.bar.pdf", width = 5, height = 2.5)
p
dev.off()




##############B cell markers##############
combined <- readRDS("../4.Analy/05.b/B.reclassified.reannot.cells.rds")

markers <- c("MS4A1", "CD19", "BANK1", "IFI30",
             "IL4R", "FCER2", "TCL1A", "BACH2",
             "IFIT3", "ISG15", "IFITM1", "JCHAIN", "MEB1",
             "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2")

combined$MinorCluster <- factor(combined$MinorCluster, c("Bm_C0", "Bn_C1", "B_C2_ISG15", "Plasma_C3_JCHAIN", "Plasma_C4_IGHA1"))
Idents(combined) <- combined$MinorCluster
pdf ("../Results/Overview2401/2401.bcell.marker.exp.pdf", height = 3, width = 10)
print(DotPlot(combined, features = markers)+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "white", mid = "#7C93C3", high = "#750E21", midpoint = 1.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

tmpres <- DotPlot(combined, features = markers)
write.table(tmpres$data, "../Results/RawData/figureS7c.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


##############Stromal cell umap##############
combined <- readRDS("../4.Analy/06.stroma/Stroma.reclassified.reannot.cells.rds")
allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$MinorCluster))]
pdf("../Results/Overview2401/2401.stroma.overview.pdf", width = 5.5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "MinorCluster", cols = color_v)+ggtitle("T cells"))
dev.off()
combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MinorCluster", "UMAP1", "UMAP2")], "../Results/RawData/figureS7d.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Stromal cell fraction##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
name_list<-c("../4.Analy/06.stroma/Stroma.reclassified.reannot.cells.rds")

combined$MinorCluster <- "Unknown"
for (file_name in name_list){
  sub_sce<-readRDS(file_name)
  ann_table<-data.frame(Cell_ID=rownames(sub_sce@meta.data),
                        MinorCluster=sub_sce$MinorCluster)
  Cluster_list<-unique(ann_table$MinorCluster)
  for (Cluster in Cluster_list){
    Cell_ID_list<-ann_table$Cell_ID[which(ann_table$MinorCluster==Cluster)]
    combined$MinorCluster[which(rownames(combined@meta.data) %in%Cell_ID_list)]<-Cluster
  }
}

Idents(combined) <- combined$MinorCluster
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Overview2401/2401.stroma."

OR.all.list <- do.tissueDist(cellInfo.tb=meta.tb,
                             meta.cluster = meta.tb$Class,
                             colname.patient = "MiratiID",
                             loc = meta.tb$MinorCluster,
                             out.prefix=sprintf("%s.Class",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../Results/Overview2401/2401.stroma.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Overview2401/2401.stroma.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac <- cellfrac[cellfrac$variable != "Unknown", ]
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Overview2401/2401.stroma.Class.diff.pdf", p1, width = 5, height = 4)


combined <- readRDS("../4.Analy/06.stroma/Stroma.reclassified.reannot.cells.rds")
##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(cellratio$Var2))]

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Overview2401/2401.stroma.Class.diff.bar.pdf", width = 5, height = 2.5)
p
dev.off()


##############Stromal cell markers##############
combined <- readRDS("../4.Analy/06.stroma/Stroma.reclassified.reannot.cells.rds")

markers <- c("PECAM1", "VWF", "CDH5",
             "PDGFRA", "FAP", "COL1A1",
             "ACKR1", "SELP", "CLU", "RGS5",
             "NDUFA4L2", "ACTA2", "RGCC", "ESM1",
             "CTHRC1", "TGFB1", "PROX1", "CCL21", "FABP5",
             "PRKG1", "SOX5")

#combined$MinorCluster <- factor(combined$MinorCluster, c("Bm_C0", "Bn_C1", "B_C2_ISG15", "Plasma_C3_JCHAIN", "Plasma_C4_IGHA1"))
Idents(combined) <- combined$MinorCluster
pdf ("../Results/Overview2401/2401.stroma.marker.exp.pdf", height = 3, width = 10)
print(DotPlot(combined, features = markers)+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "white", mid = "#7C93C3", high = "#750E21", midpoint = 1.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

tmpres <- DotPlot(combined, features = markers)
write.table(tmpres$data, "../Results/RawData/figureS7f.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)






########################################Figure S8########################################
##############tme cell statisis with heatmap##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)

combined$ModCluster <- factor(combined$ModCluster, rev(c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast")))

combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))



library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.paired"

source("/rsrch5/home/genomic_med/kyu3/scAnaly/function.R")
OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$ModCluster,
                                 out.prefix=sprintf("%s.Class.lauren",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.paired.Class.lauren.OR.dist.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.paired.Class.lauren.OR.dist.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
#cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.paired.tme.diff.Roe.pdf", p1, width = 3, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$ModCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.paired.tme.diff.fraction.pdf", width = 4, height = 3.5)
p
dev.off()



##############epi cell heatmap##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$CCluster
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.paired"

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$CCluster,
                                 out.prefix=sprintf("%s.Class.lauren",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.paired.epi.CCluster.lauren.OR.dist.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.paired.epi.CCluster.lauren.OR.dist.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.paired.epi.diff.Roe.pdf", p1, width = 3, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$CCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.paired.epi.diff.fraction.pdf", width = 4, height = 3.5)
p
dev.off()

##############selected genes##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
Idents(combined) <- combined$CCluster

combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

selmarkers <- c("B2M", "CD63", "TGFB1", "TGFBI", "SERPINE2", "H19", 
                "TIMP1", "CAV1", "HGF")

Idents(combined) <- combined$Class
pdf ("../Results/Revised07/2407.paired.epi.markers.pdf", height = 2, width = 5)
print(DotPlot(combined, features = unique(selmarkers))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "#3652AD", mid = "white", high = "#FF004D", midpoint = 0.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

tmpres <- DotPlot(combined, features = unique(selmarkers))
write.table(tmpres$data, "../Results/RawData/figureS8c.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############C8 score##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
topmarkers <- fread("../4.Analy/02.epi/harmony.epi.cluster.markers.top50.txt", header = T, stringsAsFactors = F, data.table = F)

combined <- AddModuleScore(combined,
                           features = list(topmarkers[topmarkers$cluster == "C8", ]$gene),
                           name="C8_signatures")

combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

finalres <- combined@meta.data
color_v=c("#8DD3C7", "#FB8072")
my_compa <- list(c("Baseline", "EOT"))
p <- ggplot(finalres, aes(x = Class, y = C8_signatures1, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Revised07/2407.paired.epi.c8signature.pdf", p, width = 3, height = 3)

write.table(finalres, "../Results/RawData/figureS8d.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


##############cell fraction difference##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")
combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.paired.tcell.cd4t."

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$MinorCluster,
                                 out.prefix=sprintf("%s.Class",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.paired.tcell.cd4t.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.paired.tcell.cd4t.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.paired.tcell.cd4t.Class.diff.pdf", p1, width = 3, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.paired.tcell.cd4t.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()



combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")
combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.paired.tcell.cd8t."

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$MinorCluster,
                                 out.prefix=sprintf("%s.Class",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.paired.tcell.cd8t.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.paired.tcell.cd8t.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.paired.tcell.cd8t.Class.diff.pdf", p1, width = 3, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.paired.tcell.cd8t.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()


##############signature radar plot##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")

combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

combined$Class <- as.character(combined$Class)
Idents(combined) <- combined$Class
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers[markers$p_val_adj < 0.01, ]

for (cla in unique(combined$Class)) {
  genelst <- markers[markers$cluster == cla, ]$gene
  refdat <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd4t.signature.txt", header = T, stringsAsFactors = F, data.table = F)
  
  cyto2gene <- refdat[refdat$Signature %in% c("Na\303\257ve", "Cytotoxicity", "TCR signaling"), ]
  
  cd4yanshuo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/CD4TcellMarkers.Yanshuo.top50.txt", header = T, stringsAsFactors = F, data.table = F)
  selmarkers1 <- cd4yanshuo[cd4yanshuo$cluster == "CD4_c1_Treg", ][c(1:30), c(1, 2)]
  selmarkers1$cluster <- "Treg signature"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  selmarkers1 <- cd4yanshuo[cd4yanshuo$cluster == "CD4_c4_Tstr", ][c(1:50), c(1, 2)]
  selmarkers1$cluster <- "Stress response"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  
  resdf <- enricher(genelst, TERM2GENE = cyto2gene, pvalueCutoff = 1,minGSSize = 1)
  resdf <- resdf@result
  
  resdf$ES <- as.numeric(strsplit2(resdf$GeneRatio, "/")[, 1]) / as.numeric(strsplit2(resdf$GeneRatio, "/")[, 2])
  resdf$FDR <- resdf$p.adjust
  
  write.table(resdf, paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
  
}

allclass <- data.frame(ID = c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"),
                       ES = 0,
                       FDR = 1)

##data to plot
cla = "Baseline"
res <- fread(paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res1 <- res[!duplicated(res$ID), ]
rownames(res1) <- res1$ID
res1 <- res1[c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"), ]

cla = "EOT"
res <- fread(paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res4 <- res[!duplicated(res$ID), ]
rownames(res4) <- res4$ID
res4 <- res4[c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"), ]


library(fmsb)
tmpdata <- rbind(res1$ES , res4$ES)
tmpdata = (tmpdata-min(tmpdata))/(max(tmpdata)-min(tmpdata))
data <- rbind(allclass$ES+1, allclass$ES, tmpdata)
colnames(data) <- res1$ID
rownames(data) <- c("1", "2", "Baseline", "EOT")
#data = (data-min(data))/(max(data)-min(data))
data <- as.data.frame(data)
write.table(data, "../Results/RawData/figureS8f.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# plot with default options:
library(RColorBrewer)
allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
colors_border <- c("#8DD3C7", "#FB8072")
library(scales)
colors_in <- alpha(colors_border,0.5)

pdf("../Results/Revised07/2407.paired.tcell.cd4t.signature.radar.pdf", width = 5, height = 5)
radarchart( data , axistype=1 , maxmin=T,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            caxislabels=seq(0,1,0.5), seg=length(seq(0,1,0.5))-1,
            #custom labels
            vlcex=0.8 
)
# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
dev.off()




combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")
combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

combined$Class <- as.character(combined$Class)
Idents(combined) <- combined$Class
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers[markers$p_val_adj < 0.01, ]

for (cla in unique(combined$Class)) {
  genelst <- markers[markers$cluster == cla, ]$gene
  refdat <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd8t.signature.txt", header = T, stringsAsFactors = F, data.table = F)
  
  cd4yanshuo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/CD8TcellMarkers.Yanshuo.top50.txt", header = T, stringsAsFactors = F, data.table = F)
  selmarkers1 <- cd4yanshuo[cd4yanshuo$icluster == "CD8_c1_Tex", ][c(1:30), c(1, 2)]
  selmarkers1$icluster <- "Exhaustion"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  selmarkers1 <- cd4yanshuo[cd4yanshuo$icluster == "CD8_c4_Tstr", ][c(1:50), c(1, 2)]
  selmarkers1$icluster <- "Stress response"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  cyto2gene <- refdat[refdat$Signature %in% c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]
  resdf <- enricher(genelst, TERM2GENE = cyto2gene, pvalueCutoff = 1,minGSSize = 1)
  resdf <- resdf@result
  
  resdf$ES <- as.numeric(strsplit2(resdf$GeneRatio, "/")[, 1]) / as.numeric(strsplit2(resdf$GeneRatio, "/")[, 2])
  resdf$FDR <- resdf$p.adjust
  
  write.table(resdf, paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
  
}

allclass <- data.frame(ID = c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"),
                       ES = 0,
                       FDR = 1)

##data to plot
cla = "Baseline"
res <- fread(paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res1 <- res[!duplicated(res$ID), ]
rownames(res1) <- res1$ID
res1 <- res1[c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]

cla = "EOT"
res <- fread(paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res4 <- res[!duplicated(res$ID), ]
rownames(res4) <- res4$ID
res4 <- res4[c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]


library(fmsb)
tmpdata <- rbind(res1$ES , res4$ES)
tmpdata = (tmpdata-min(tmpdata))/(max(tmpdata)-min(tmpdata))
data <- rbind(allclass$ES+1, allclass$ES, tmpdata)
colnames(data) <- res1$ID
rownames(data) <- c("1", "2", "Baseline", "EOT")
#data = (data-min(data))/(max(data)-min(data))
data <- as.data.frame(data)

write.table(data, "../Results/RawData/figureS8h.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# plot with default options:
library(RColorBrewer)
allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
colors_border <- c("#8DD3C7", "#FB8072")
library(scales)
colors_in <- alpha(colors_border,0.5)

pdf("../Results/Revised07/2407.paired.tcell.cd8t.signature.radar.pdf", width = 5, height = 5)
radarchart( data , axistype=1 , maxmin=T,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            caxislabels=seq(0,1,0.5), seg=length(seq(0,1,0.5))-1,
            #custom labels
            vlcex=0.8 
)
# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
dev.off()


##############cell fraction difference##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")
combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.paired.myeloid."

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$MinorCluster,
                                 out.prefix=sprintf("%s.Class",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.paired.myeloid.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.paired.myeloid.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.paired.myeloid.Class.diff.pdf", p1, width = 4, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")
cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.paired.myeloid.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()







##############signature score##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")
allscores <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.txt", header = F)
allscores <- split(allscores$V2, allscores$V1)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))


allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")
selgenes <- names(allscores)
resdat <- combined@meta.data[, c("Class", selgenes)]
resdat <- melt(resdat, id.vars = "Class")
names(resdat) <- c("Class", "MPs", "Score")

resdat <- resdat[resdat$MPs == "M2", ]

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

my_compa <- list(c("Baseline", "EOT"))
p <- ggplot(resdat, aes(x = Class, y = Score, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(p, width = 3, height = 3, filename = "../Results/Revised07/2407.paired.myeloid.overview.signature.vlnplot.pdf")

write.table(resdat, "../Results/RawData/figureS8j.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############CSF1R expression##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

combined$CSF1R <- FetchData(combined, vars = "CSF1R")$CSF1R
finalres <- combined@meta.data

color_v <- c("#8DD3C7", "#FB8072")
my_compa <- list(c("Baseline", "EOT"))
p <- ggplot(finalres, aes(x = Class, y = CSF1R, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Revised07/2407.paired.myeloid.class.csf1r.pdf", p, width = 3, height = 3)

write.table(resdat, "../Results/RawData/figureS8k.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############IL6pathway##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

allscores <- split(hall$gene_symbol, hall$gs_name)
#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

combined$Class <- as.character(combined$Class)
combined$Class <- ifelse(combined$Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"), "Baseline", combined$Class)
combined <- subset(combined, subset = Class != "C2")
combined <- subset(combined, subset = MiratiID %in% c("RCC16", "RCC18", "RCC20"))
combined$Class <- factor(x = combined$Class, levels = c("Baseline", "EOT"))

finalres <- combined@meta.data
color_v=c("#8DD3C7", "#FB8072")
my_compa <- list(c("Baseline", "EOT"))
p <- ggplot(finalres, aes(x = Class, y = HALLMARK_IL6_JAK_STAT3_SIGNALING, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Revised07/2407.paired.myeloid.il6.jak.pdf", p, width = 3, height = 3)

write.table(resdat, "../Results/RawData/figureS8l.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

########################################Figure S9########################################
##############tme cell statisis with heatmap##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)

combined$ModCluster <- factor(combined$ModCluster, rev(c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast")))

combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.dose4"

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$ModCluster,
                                 out.prefix=sprintf("%s.Class.lauren",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.dose4.Class.lauren.OR.dist.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.dose4.Class.lauren.OR.dist.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.dose4.tme.diff.Roe.pdf", p1, width = 5, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$ModCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(cellratio$Var2))]

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.dose4.tme.diff.fraction.pdf", width = 4, height = 3.5)
p
dev.off()


##############tme cell fraction with orr##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)
combined$ModCluster <- factor(combined$ModCluster, c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast"))
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)
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

write.table(res, "../Results/Revised07/2407.dose4.tme.cellfraction.lm.r2.pva.txt", row.names = F, col.names = T, quote = F, sep = "\t")

res$logp <- -log10(res$Pva)
allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$ModCluster))]

res$logp <- 10*res$logp

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

ggsave("../Results/Revised07/2407.dose4.tme.cellfraction.lm.r2.pva.pdf", p, width = 4, height = 3)


##############epi cell heatmap##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$CCluster
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)
library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.dose4"

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$CCluster,
                                 out.prefix=sprintf("%s.Class.lauren",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.dose4.epi.CCluster.lauren.OR.dist.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.dose4.epi.CCluster.lauren.OR.dist.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.dose4.epi.diff.Roe.pdf", p1, width = 5, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$CCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(cellratio$Var2))]

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.dose4.epi.diff.fraction.pdf", width = 4, height = 3.5)
p
dev.off()



##############epi cell fraction with orr##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$CCluster
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)
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

write.table(res, "../Results/Revised07/2407.dose4.epi.cellfraction.lm.r2.pva.txt", row.names = F, col.names = T, quote = F, sep = "\t")

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

ggsave("../Results/Revised07/2407.dose4.epi.cellfraction.lm.r2.pva.pdf", p, width = 3.5, height = 2.5)



##############selected genes##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
Idents(combined) <- combined$CCluster
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)
combined <- subset(combined, subset = CCluster %in% c("C2", "C8"))

selmarkers <- c("B2M", "CD63", "TGFB1", "TGFBI", "SERPINE2", "H19", 
                "TIMP1", "CAV1", "HGF")

Idents(combined) <- combined$ModCluster
pdf ("../Results/Revised07/2407.dose4.epi.markers.pdf", height = 2, width = 5)
print(DotPlot(combined, features = unique(selmarkers))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "#3652AD", mid = "white", high = "#FF004D", midpoint = 0.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

tmpres <- DotPlot(combined, features = unique(selmarkers))
write.table(tmpres$data, "../Results/RawData/figureS8e.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############selected regulon plot##############
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
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)
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
pdf("../Results/Revised07/2407.dose4.epi.scenic.tfs.heatmap.ccluster.selected.pdf", width = 6, height = 2)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[c("TWIST1(+),", "ZEB1(+),"), ], name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=6))
dev.off()

write.table(regulonActivity_byCellType_Scaled[c("TWIST1(+),", "ZEB1(+),"), ], "../Results/RawData/figureS8f.txt",
            row.names = T, col.names = NA, sep = "\t", quote = F)


##############MP EMT##############
#####meta program
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)

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

ggsave(q, width = 12, height = 4, filename = "../Results/Revised07/2407.dose4.epi.cluster.selmetaprograms.pdf")

write.table(resdat, "../Results/RawData/figureS8g.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############cell fraction difference##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")
#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.dose4.tcell.cd4t."

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$MinorCluster,
                                 out.prefix=sprintf("%s.Class",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.dose4.tcell.cd4t.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.dose4.tcell.cd4t.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.dose4.tcell.cd4t.Class.diff.pdf", p1, width = 5, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(cellratio$Var2))]

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.dose4.tcell.cd4t.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()



combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")
#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.dose4.tcell.cd8t."

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$MinorCluster,
                                 out.prefix=sprintf("%s.Class",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.dose4.tcell.cd8t.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.dose4.tcell.cd8t.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.dose4.tcell.cd8t.Class.diff.pdf", p1, width = 5, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(cellratio$Var2))]

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.dose4.tcell.cd8t.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()


##############signature radar plot##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")
combined$Class <- as.character(combined$Class)

combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)


Idents(combined) <- combined$Class
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers[markers$p_val_adj < 0.01, ]

for (cla in unique(combined$Class)) {
  genelst <- markers[markers$cluster == cla, ]$gene
  refdat <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd4t.signature.txt", header = T, stringsAsFactors = F, data.table = F)
  
  cyto2gene <- refdat[refdat$Signature %in% c("Na\303\257ve", "Cytotoxicity", "TCR signaling"), ]
  
  cd4yanshuo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/CD4TcellMarkers.Yanshuo.top50.txt", header = T, stringsAsFactors = F, data.table = F)
  selmarkers1 <- cd4yanshuo[cd4yanshuo$cluster == "CD4_c1_Treg", ][c(1:30), c(1, 2)]
  selmarkers1$cluster <- "Treg signature"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  selmarkers1 <- cd4yanshuo[cd4yanshuo$cluster == "CD4_c4_Tstr", ][c(1:50), c(1, 2)]
  selmarkers1$cluster <- "Stress response"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  
  resdf <- enricher(genelst, TERM2GENE = cyto2gene, pvalueCutoff = 1,minGSSize = 1)
  resdf <- resdf@result
  
  resdf$ES <- as.numeric(strsplit2(resdf$GeneRatio, "/")[, 1]) / as.numeric(strsplit2(resdf$GeneRatio, "/")[, 2])
  resdf$FDR <- resdf$p.adjust
  
  write.table(resdf, paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
  
}

allclass <- data.frame(ID = c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"),
                       ES = 0,
                       FDR = 1)

##data to plot
cla = "Baseline_PR/CR"
res <- fread(paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res1 <- res[!duplicated(res$ID), ]
rownames(res1) <- res1$ID
res1 <- res1[c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"), ]

cla = "Baseline_PD/SD"
res <- fread(paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res2 <- res[!duplicated(res$ID), ]
rownames(res2) <- res2$ID
res2 <- res2[c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"), ]

cla = "C2"
res <- fread(paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res3 <- res[!duplicated(res$ID), ]
rownames(res3) <- res3$ID
res3 <- res3[c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"), ]

cla = "EOT"
res <- fread(paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res4 <- res[!duplicated(res$ID), ]
rownames(res4) <- res4$ID
res4 <- res4[c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"), ]


library(fmsb)
tmpdata <- rbind(res1$ES , res2$ES, res3$ES, res4$ES)
tmpdata = (tmpdata-min(tmpdata))/(max(tmpdata)-min(tmpdata))
data <- rbind(allclass$ES+1, allclass$ES, tmpdata)
colnames(data) <- res1$ID
rownames(data) <- c("1", "2", "Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT")
#data = (data-min(data))/(max(data)-min(data))
data <- as.data.frame(data)
write.table(data, "../Results/RawData/figureS9i.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# plot with default options:
library(RColorBrewer)
allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
colors_border <- allcolors[1:4]
library(scales)
colors_in <- alpha(colors_border,0.5)

pdf("../Results/Revised07/2407.dose4.tcell.cd4t.signature.radar.pdf", width = 5, height = 5)
radarchart( data , axistype=1 , maxmin=T,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            caxislabels=seq(0,1,0.5), seg=length(seq(0,1,0.5))-1,
            #custom labels
            vlcex=0.8 
)
# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
dev.off()




combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)


combined$Class <- as.character(combined$Class)
Idents(combined) <- combined$Class
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers[markers$p_val_adj < 0.01, ]

for (cla in unique(combined$Class)) {
  genelst <- markers[markers$cluster == cla, ]$gene
  refdat <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd8t.signature.txt", header = T, stringsAsFactors = F, data.table = F)
  
  cd4yanshuo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/CD8TcellMarkers.Yanshuo.top50.txt", header = T, stringsAsFactors = F, data.table = F)
  selmarkers1 <- cd4yanshuo[cd4yanshuo$icluster == "CD8_c1_Tex", ][c(1:30), c(1, 2)]
  selmarkers1$icluster <- "Exhaustion"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  selmarkers1 <- cd4yanshuo[cd4yanshuo$icluster == "CD8_c4_Tstr", ][c(1:50), c(1, 2)]
  selmarkers1$icluster <- "Stress response"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  cyto2gene <- refdat[refdat$Signature %in% c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]
  resdf <- enricher(genelst, TERM2GENE = cyto2gene, pvalueCutoff = 1,minGSSize = 1)
  resdf <- resdf@result
  
  resdf$ES <- as.numeric(strsplit2(resdf$GeneRatio, "/")[, 1]) / as.numeric(strsplit2(resdf$GeneRatio, "/")[, 2])
  resdf$FDR <- resdf$p.adjust
  
  write.table(resdf, paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
  
}

allclass <- data.frame(ID = c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"),
                       ES = 0,
                       FDR = 1)

##data to plot
cla = "Baseline_PR/CR"
res <- fread(paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res1 <- res[!duplicated(res$ID), ]
rownames(res1) <- res1$ID
res1 <- res1[c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]

cla = "Baseline_PD/SD"
res <- fread(paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res2 <- res[!duplicated(res$ID), ]
rownames(res2) <- res2$ID
res2 <- res2[c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]

cla = "C2"
res <- fread(paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res3 <- res[!duplicated(res$ID), ]
rownames(res3) <- res3$ID
res3 <- res3[c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]

cla = "EOT"
res <- fread(paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res4 <- res[!duplicated(res$ID), ]
rownames(res4) <- res4$ID
res4 <- res4[c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]


library(fmsb)
tmpdata <- rbind(res1$ES , res2$ES, res3$ES, res4$ES)
tmpdata = (tmpdata-min(tmpdata))/(max(tmpdata)-min(tmpdata))
data <- rbind(allclass$ES+1, allclass$ES, tmpdata)
colnames(data) <- res1$ID
rownames(data) <- c("1", "2", "Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT")
#data = (data-min(data))/(max(data)-min(data))
data <- as.data.frame(data)

write.table(data, "../Results/RawData/figureS9k.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# plot with default options:
library(RColorBrewer)
allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
colors_border <- allcolors[1:4]
library(scales)
colors_in <- alpha(colors_border,0.5)

pdf("../Results/Revised07/2407.dose4.tcell.cd8t.signature.radar.pdf", width = 5, height = 5)
radarchart( data , axistype=1 , maxmin=T,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            caxislabels=seq(0,1,0.5), seg=length(seq(0,1,0.5))-1,
            #custom labels
            vlcex=0.8 
)
# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
dev.off()



##############cell fraction difference##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")
#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.dose4.myeloid."

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$MinorCluster,
                                 out.prefix=sprintf("%s.Class",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.dose4.myeloid.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.dose4.myeloid.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.dose4.myeloid.Class.diff.pdf", p1, width = 5, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(cellratio$Var2))]

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.dose4.myeloid.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()



##############CSF1R-tumor response##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)

combined$CSF1R <- FetchData(combined, vars = "CSF1R")$CSF1R
finalres <- combined@meta.data

color_v <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072")
my_compa <- list(c("Baseline_PR/CR", "Baseline_PD/SD"), c("Baseline_PR/CR", "EOT"))
p <- ggplot(finalres, aes(x = Class, y = CSF1R, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Revised07/2407.dose4.myeloid.class.csf1r.pdf", p, width = 4, height = 3)

write.table(resdat, "../Results/RawData/figureS9m.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############M2-tumor response##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")
allscores <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.txt", header = F)
allscores <- split(allscores$V2, allscores$V1)
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)


finalres <- combined@meta.data

color_v <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072")
my_compa <- list(c("Baseline_PR/CR", "Baseline_PD/SD"), c("Baseline_PR/CR", "EOT"))
p <- ggplot(finalres, aes(x = Class, y = M2, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Revised07/2407.dose4.myeloid.class.m2.pdf", p, width = 4, height = 3)


##############IL6pathway##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")
allscores <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.txt", header = F)
allscores <- split(allscores$V2, allscores$V1)
combined <- subset(combined, subset = MiratiID %in% c("RCC12", "RCC13", "RCC14", "RCC16"), invert = T)

hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

allscores <- split(hall$gene_symbol, hall$gs_name)
#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)


finalres <- combined@meta.data
color_v <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072")
my_compa <- list(c("Baseline_PR/CR", "Baseline_PD/SD"), c("Baseline_PR/CR", "EOT"))
p <- ggplot(finalres, aes(x = Class, y = HALLMARK_IL6_JAK_STAT3_SIGNALING, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Revised07/2407.dose4.myeloid.il6.jak.pdf", p, width = 4, height = 3)

########################################Figure S10########################################
##############tme cell statisis with heatmap##############
combined <- readRDS("../Results/Overview2311/2312.overview.tme.rds")
combined$ModCluster <- ifelse(combined$ModCluster %in% c("T", "NK"), "T/NK", combined$ModCluster)

combined$ModCluster <- factor(combined$ModCluster, rev(c("T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast")))

combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))



library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.sarco"

source("/rsrch5/home/genomic_med/kyu3/scAnaly/function.R")
OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$ModCluster,
                                 out.prefix=sprintf("%s.Class.lauren",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.sarco.Class.lauren.OR.dist.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.sarco.Class.lauren.OR.dist.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
#cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Others", "Sar_Rha"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.sarco.tme.diff.Roe.pdf", p1, width = 3, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$ModCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Others", "Sar_Rha"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.sarco.tme.diff.fraction.pdf", width = 4, height = 3.5)
p
dev.off()



##############epi cell heatmap##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$CCluster
combined$Class <- "Unknown"
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PR/CR", "Baseline_PR/CR", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "Baseline" & combined$PRCR == "PD/SD", "Baseline_PD/SD", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "C2", "C2", combined$Class)
combined$Class <- ifelse(combined$Collection.point == "EOT", "EOT", combined$Class)
combined$Class <- factor(x = combined$Class, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.sarco"

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$CCluster,
                                 out.prefix=sprintf("%s.Class.lauren",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.sarco.epi.CCluster.lauren.OR.dist.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.sarco.epi.CCluster.lauren.OR.dist.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Others", "Sar_Rha"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.sarco.epi.diff.Roe.pdf", p1, width = 3, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$CCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Others", "Sar_Rha"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.sarco.epi.diff.fraction.pdf", width = 4, height = 3.5)
p
dev.off()

##############selected genes##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
Idents(combined) <- combined$CCluster

combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

selmarkers <- c("B2M", "CD63", "TGFB1", "TGFBI", "SERPINE2", "H19", 
                "TIMP1", "CAV1", "HGF")

Idents(combined) <- combined$Class
pdf ("../Results/Revised07/2407.sarco.epi.markers.pdf", height = 2, width = 5)
print(DotPlot(combined, features = unique(selmarkers))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "#3652AD", mid = "white", high = "#FF004D", midpoint = 0.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()
tmpres <- DotPlot(combined, features = unique(selmarkers))
write.table(tmpres$data, "../Results/RawData/figureS10c.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############C8 score##############
combined <- readRDS("../4.Analy/02.epi/harmony.Merged.pca.for.monocle3.rds")
topmarkers <- fread("../4.Analy/02.epi/harmony.epi.cluster.markers.top50.txt", header = T, stringsAsFactors = F, data.table = F)

combined <- AddModuleScore(combined,
                           features = list(topmarkers[topmarkers$cluster == "C8", ]$gene),
                           name="C8_signatures")

combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

finalres <- combined@meta.data
color_v=c("#8DD3C7", "#FB8072")
my_compa <- list(c("Others", "Sar_Rha"))
p <- ggplot(finalres, aes(x = Class, y = C8_signatures1, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Revised07/2407.sarco.epi.c8signature.pdf", p, width = 3, height = 3)

write.table(finalres, "../Results/RawData/figureS10d.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)



##############cell fraction difference##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")
combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.sarco.tcell.cd4t."

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$MinorCluster,
                                 out.prefix=sprintf("%s.Class",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.sarco.tcell.cd4t.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.sarco.tcell.cd4t.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Others", "Sar_Rha"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.sarco.tcell.cd4t.Class.diff.pdf", p1, width = 3, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Others", "Sar_Rha"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.sarco.tcell.cd4t.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()



combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")
combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.sarco.tcell.cd8t."

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$MinorCluster,
                                 out.prefix=sprintf("%s.Class",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.sarco.tcell.cd8t.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.sarco.tcell.cd8t.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Others", "Sar_Rha"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.sarco.tcell.cd8t.Class.diff.pdf", p1, width = 3, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Others", "Sar_Rha"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.sarco.tcell.cd8t.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()


##############signature radar plot##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")

combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

combined$Class <- as.character(combined$Class)
Idents(combined) <- combined$Class
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers[markers$p_val_adj < 0.01, ]

for (cla in unique(combined$Class)) {
  genelst <- markers[markers$cluster == cla, ]$gene
  refdat <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd4t.signature.txt", header = T, stringsAsFactors = F, data.table = F)
  
  cyto2gene <- refdat[refdat$Signature %in% c("Na\303\257ve", "Cytotoxicity", "TCR signaling"), ]
  
  cd4yanshuo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/CD4TcellMarkers.Yanshuo.top50.txt", header = T, stringsAsFactors = F, data.table = F)
  selmarkers1 <- cd4yanshuo[cd4yanshuo$cluster == "CD4_c1_Treg", ][c(1:30), c(1, 2)]
  selmarkers1$cluster <- "Treg signature"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  selmarkers1 <- cd4yanshuo[cd4yanshuo$cluster == "CD4_c4_Tstr", ][c(1:50), c(1, 2)]
  selmarkers1$cluster <- "Stress response"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  
  resdf <- enricher(genelst, TERM2GENE = cyto2gene, pvalueCutoff = 1,minGSSize = 1)
  resdf <- resdf@result
  
  resdf$ES <- as.numeric(strsplit2(resdf$GeneRatio, "/")[, 1]) / as.numeric(strsplit2(resdf$GeneRatio, "/")[, 2])
  resdf$FDR <- resdf$p.adjust
  
  write.table(resdf, paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
  
}

allclass <- data.frame(ID = c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"),
                       ES = 0,
                       FDR = 1)

##data to plot
cla = "Baseline"
res <- fread(paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res1 <- res[!duplicated(res$ID), ]
rownames(res1) <- res1$ID
res1 <- res1[c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"), ]

cla = "EOT"
res <- fread(paste("../4.Analy/03.t/CD4T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res4 <- res[!duplicated(res$ID), ]
rownames(res4) <- res4$ID
res4 <- res4[c("Na\303\257ve", "Treg signature", "Cytotoxicity", "Stress response", "TCR signaling"), ]


library(fmsb)
tmpdata <- rbind(res1$ES , res4$ES)
tmpdata = (tmpdata-min(tmpdata))/(max(tmpdata)-min(tmpdata))
data <- rbind(allclass$ES+1, allclass$ES, tmpdata)
colnames(data) <- res1$ID
rownames(data) <- c("1", "2", "Baseline", "EOT")
#data = (data-min(data))/(max(data)-min(data))
data <- as.data.frame(data)

write.table(data, "../Results/RawData/figureS10f.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# plot with default options:
library(RColorBrewer)
allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
colors_border <- c("#8DD3C7", "#FB8072")
library(scales)
colors_in <- alpha(colors_border,0.5)

pdf("../Results/Revised07/2407.sarco.tcell.cd4t.signature.radar.pdf", width = 5, height = 5)
radarchart( data , axistype=1 , maxmin=T,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            caxislabels=seq(0,1,0.5), seg=length(seq(0,1,0.5))-1,
            #custom labels
            vlcex=0.8 
)
# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
dev.off()




combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")
combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

combined$Class <- as.character(combined$Class)
Idents(combined) <- combined$Class
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers[markers$p_val_adj < 0.01, ]

for (cla in unique(combined$Class)) {
  genelst <- markers[markers$cluster == cla, ]$gene
  refdat <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd8t.signature.txt", header = T, stringsAsFactors = F, data.table = F)
  
  cd4yanshuo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/CD8TcellMarkers.Yanshuo.top50.txt", header = T, stringsAsFactors = F, data.table = F)
  selmarkers1 <- cd4yanshuo[cd4yanshuo$icluster == "CD8_c1_Tex", ][c(1:30), c(1, 2)]
  selmarkers1$icluster <- "Exhaustion"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  selmarkers1 <- cd4yanshuo[cd4yanshuo$icluster == "CD8_c4_Tstr", ][c(1:50), c(1, 2)]
  selmarkers1$icluster <- "Stress response"
  names(selmarkers1) <- c("Signature", "Genes")
  cyto2gene <- rbind(cyto2gene, selmarkers1)
  
  cyto2gene <- refdat[refdat$Signature %in% c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]
  resdf <- enricher(genelst, TERM2GENE = cyto2gene, pvalueCutoff = 1,minGSSize = 1)
  resdf <- resdf@result
  
  resdf$ES <- as.numeric(strsplit2(resdf$GeneRatio, "/")[, 1]) / as.numeric(strsplit2(resdf$GeneRatio, "/")[, 2])
  resdf$FDR <- resdf$p.adjust
  
  write.table(resdf, paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
  
}

allclass <- data.frame(ID = c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"),
                       ES = 0,
                       FDR = 1)

##data to plot
cla = "Baseline"
res <- fread(paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res1 <- res[!duplicated(res$ID), ]
rownames(res1) <- res1$ID
res1 <- res1[c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]

cla = "EOT"
res <- fread(paste("../4.Analy/03.t/CD8T.", gsub("/", ".", cla), ".cellpolar.genelist.enrich.res.txt", sep = ""), header = T, stringsAsFactors = F, data.table = F)
res <- res[, c("ID", "ES", "FDR")]
res <- rbind(res, allclass)
res4 <- res[!duplicated(res$ID), ]
rownames(res4) <- res4$ID
res4 <- res4[c("Na\303\257ve", "Exhaustion", "Cytotoxicity", "Stress response", "TCR Signaling"), ]


library(fmsb)
tmpdata <- rbind(res1$ES , res4$ES)
tmpdata = (tmpdata-min(tmpdata))/(max(tmpdata)-min(tmpdata))
data <- rbind(allclass$ES+1, allclass$ES, tmpdata)
colnames(data) <- res1$ID
rownames(data) <- c("1", "2", "Baseline", "EOT")
#data = (data-min(data))/(max(data)-min(data))
data <- as.data.frame(data)
write.table(data, "../Results/RawData/figureS10h.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# plot with default options:
library(RColorBrewer)
allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
colors_border <- c("#8DD3C7", "#FB8072")
library(scales)
colors_in <- alpha(colors_border,0.5)

pdf("../Results/Revised07/2407.sarco.tcell.cd8t.signature.radar.pdf", width = 5, height = 5)
radarchart( data , axistype=1 , maxmin=T,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            caxislabels=seq(0,1,0.5), seg=length(seq(0,1,0.5))-1,
            #custom labels
            vlcex=0.8 
)
# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
dev.off()


##############cell fraction difference##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")
combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Revised07/2407.sarco.myeloid."

OR.all.list <- do.tissueDist.tmp(cellInfo.tb=meta.tb,
                                 meta.cluster = meta.tb$Class,
                                 colname.patient = "MiratiID",
                                 loc = meta.tb$MinorCluster,
                                 out.prefix=sprintf("%s.Class",out.prefix),
                                 pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list), "../Results/Revised07/2407.sarco.myeloid.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Revised07/2407.sarco.myeloid.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Others", "Sar_Rha"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Revised07/2407.sarco.myeloid.Class.diff.pdf", p1, width = 4, height = 4)


##cell fraction
cellratio <- prop.table(table(combined$MinorCluster, as.factor(as.character(combined$Class))), margin = 1)
cellratio <- as.data.frame(cellratio)
c <- length(unique(cellratio$Var1))

allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")
cellratio$Var2 <- factor(x = cellratio$Var2, levels = c("Others", "Sar_Rha"))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  scale_fill_manual(values = color_v)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf("../Results/Revised07/2407.sarco.myeloid.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()







##############signature score##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")
allscores <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.txt", header = F)
allscores <- split(allscores$V2, allscores$V1)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))


allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")
selgenes <- names(allscores)
resdat <- combined@meta.data[, c("Class", selgenes)]
resdat <- melt(resdat, id.vars = "Class")
names(resdat) <- c("Class", "MPs", "Score")

resdat <- resdat[resdat$MPs == "M2", ]

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=c("#8DD3C7", "#FB8072")

my_compa <- list(c("Others", "Sar_Rha"))
p <- ggplot(resdat, aes(x = Class, y = Score, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(p, width = 3, height = 3, filename = "../Results/Revised07/2407.sarco.myeloid.overview.signature.vlnplot.pdf")
write.table(resdat, "../Results/RawData/figureS10j.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############CSF1R expression##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

combined$CSF1R <- FetchData(combined, vars = "CSF1R")$CSF1R
finalres <- combined@meta.data

color_v <- c("#8DD3C7", "#FB8072")
my_compa <- list(c("Others", "Sar_Rha"))
p <- ggplot(finalres, aes(x = Class, y = CSF1R, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Revised07/2407.sarco.myeloid.class.csf1r.pdf", p, width = 3, height = 3)

write.table(finalres, "../Results/RawData/figureS10k.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############IL6pathway##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

allscores <- split(hall$gene_symbol, hall$gs_name)
#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

combined <- subset(combined, subset = Class %in% c("Baseline_PD/SD", "Baseline_PR/CR"))
combined$Class <- ifelse(combined$MiratiID %in% c("RCC12", "RCC18", "RCC27"), "Sar_Rha", "Others")
combined$Class <- factor(x = combined$Class, levels = c("Others", "Sar_Rha"))

finalres <- combined@meta.data
color_v=c("#8DD3C7", "#FB8072")
my_compa <- list(c("Others", "Sar_Rha"))
p <- ggplot(finalres, aes(x = Class, y = HALLMARK_IL6_JAK_STAT3_SIGNALING, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Revised07/2407.sarco.myeloid.il6.jak.pdf", p, width = 3, height = 3)

write.table(finalres, "../Results/RawData/figureS10l.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)





