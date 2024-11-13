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


############################Figure 6############################
##############Figure 6a##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")
allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$MinorCluster))]
pdf("../Results/Overview2401/2401.myeloid.overview.pdf", width = 5.5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "MinorCluster", cols = color_v)+ggtitle("T cells"))
dev.off()

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MinorCluster", "UMAP1", "UMAP2")], "../Results/RawData/figure6a.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Figure 6b##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")
#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Overview2401/2401.myeloid."

OR.all.list <- do.tissueDist(cellInfo.tb=meta.tb,
                             meta.cluster = meta.tb$Class,
                             colname.patient = "MiratiID",
                             loc = meta.tb$MinorCluster,
                             out.prefix=sprintf("%s.Class",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../Results/Overview2401/2401.myeloid.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Overview2401/2401.myeloid.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Overview2401/2401.myeloid.Class.diff.pdf", p1, width = 5, height = 4)


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

pdf("../Results/Overview2401/2401.myeloid.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()




##############Figure 6c##############
combined <- readRDS("../4.Analy/04.myeloid/Myeloid.reclassified.reAnnot.cells.rds")
tempObj <- combined
Idents(tempObj) <- tempObj$MinorCluster
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
color_v=allcolors[1:length(unique(meta$MinorCluster))]
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
pdf("../Results/Overview2401/2401.myeloid.umap.desity.pdf", width = 12, height = 3)
print(g)
dev.off()
g <- ggplot(data = meta, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = MinorCluster), size = 0.2) +
  scale_color_manual(values = color_v) +
  facet_wrap(~Class,ncol = 4) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.position = "none")
pdf("../Results/Overview2401/2401.myeloid.umap.group.pdf", width = 14, height = 3)
print(g)
dev.off()



##############Figure 6d##############
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
write.table(combined@meta.data[, c("ID", "MinorCluster", "M1", "M2","Angiogenesis", "Phagocytosis")], "../Results/RawData/figure6d.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


selgenes <- names(allscores)
resdat <- combined@meta.data[, c("MinorCluster", selgenes)]
resdat <- melt(resdat, id.vars = "MinorCluster")
names(resdat) <- c("MinorCluster", "MPs", "Score")

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$MinorCluster))]

p <- ggplot(resdat, aes(x = MinorCluster, y = Score, fill = MinorCluster))+
  geom_boxplot(outlier.size = 0)+
  scale_fill_manual(values = color_v) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
q <- p + facet_wrap(vars(MPs), scales = "free_y", ncol = 2)

ggsave(q, width = 8, height = 4, filename = "../Results/Overview2401/2401.myeloid.overview.signature.vlnplot.pdf")



##############Figure 6e##############
cellchat <- readRDS("../Results/Overview2401/2401.cellchat.res.rds")
# 提取特定细胞类型间通讯关系
df.net2 <- subsetCommunication(cellchat, sources.use=c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'
                                                       , 'C8', 'C9', 'C10') , targets.use= c('CD8T_C1_Tex','CD4T_C2_Treg', 'TAM'))
head(df.net2)
write.table(df.net2, "../Results/RawData/figure6e.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# 提取特定信号通路中的通讯关系
#df.net3 <- subsetCommunication(cellchat, signaling = c("IL16", "ANNEXIN"))
#head(df.net3)

# NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat <- computeCommunProbPathway(cellchat)
cellchat@netP

# 指定细胞类型
#netVisual_bubble(cellchat, sources.use = c(1:5), targets.use = c(3:6), remove.isolate = FALSE )
pdf("../Results/Overview2401/2401.cellchat.bubble.pdf")
netVisual_bubble(cellchat, sources.use = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'
                                           , 'C8', 'C9', 'C10'), targets.use = c('CD8T_C1_Tex','CD4T_C2_Treg', 'TAM'), remove.isolate = FALSE )
dev.off()

pdf("../Results/Overview2401/2401.cellchat.bubble.sel.epi.pdf")
netVisual_bubble(cellchat, sources.use = c('C2', 'C3', 'C7'
                                           , 'C8', 'C5', 'C10'), targets.use = c('CD8T_C1_Tex','CD4T_C2_Treg', 'TAM'),
                 sort.by.source = TRUE, remove.isolate = FALSE )
de
##############Figure 6f##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")
allscores <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/myeloid.signature.txt", header = F)
allscores <- split(allscores$V2, allscores$V1)

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

ggsave("../Results/Overview2407/2407.myeloid.class.m2.pdf", p, width = 4, height = 3.5)

write.table(finalres[, c("ID", "Class", "M2")], "../Results/RawData/figure6f.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Figure 6g##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

combined$CSF1R <- FetchData(combined, vars = "CSF1R")$CSF1R
finalres <- combined@meta.data

allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(sub_sce$MinorCluster))]

color_v <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072")
my_compa <- list(c("Baseline_PR/CR", "Baseline_PD/SD"), c("Baseline_PR/CR", "EOT"))
p <- ggplot(finalres, aes(x = Class, y = CSF1R, fill = Class))+
  geom_violin()+
  geom_boxplot(outlier.colour = NA, width = 0.2)+
  scale_fill_manual(values = color_v)+
  stat_compare_means(method = "t.test", comparisons = my_compa, label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("../Results/Overview2407/2407.myeloid.class.csf1r.pdf", p, width = 4, height = 3.5)

write.table(finalres[, c("ID", "Class", "CSF1R")], "../Results/RawData/figure6g.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)



##############Figure 6h##############
combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

allscores <- split(hall$gene_symbol, hall$gs_name)
#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

expmtr <- combined@meta.data[, c("Class", names(allscores))]
expmtr <- aggregate(.~ Class, expmtr, median)
rownames(expmtr) <- expmtr$Class
expmtr <- t(expmtr[, -1])

rowanno <- data.frame(Color = c("Baseline_PRCR", "Baseline_PDSD", "C2", "EOT"))
rownames(rowanno) <- c("Baseline_PRCR", "Baseline_PDSD", "C2", "EOT")
expmtr <- expmtr[, c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT")]
colnames(expmtr) <- gsub("\\/", "", colnames(expmtr))

get_top5_nonredundant <- function(df) {
  df_long <- df %>% 
    pivot_longer(cols = -Pathway, names_to = "Category", values_to = "Value") %>% 
    group_by(Category) %>% 
    arrange(desc(Value)) %>% 
    slice_head(n = 5)
  
  df_wide <- df_long %>% 
    pivot_wider(names_from = Category, values_from = Value)
  
  df_wide <- df_wide %>% select(Pathway, all_of(names(df)[-which(names(df) == "Pathway")]))
  return(df_wide)
}
data <- as.data.frame(t(scale(t(expmtr))))
data$Pathway <- rownames(data)
write.table(data, "../Results/Overview2407/2407.myeloid.class.hallmark.zscore.txt", row.names = F, col.names = T, quote = F, sep = "\t")

data$Pathway <- rownames(data)
# Get the top 2 non-redundant pathways
top5_pathways <- get_top5_nonredundant(data)
# Print the result
allpathes <- c()
for (col in colnames(data)[-which(colnames(data) == "Pathway")]) {
  cat("Top 5 pathways for", col, ":\n")
  selected_pathways <- top5_pathways %>% filter(!is.na(.data[[col]])) %>% select(Pathway, all_of(col))
  print(selected_pathways)
  cat("\n")
  allpathes <- c(allpathes, selected_pathways$Pathway)
}


finalres <- expmtr[allpathes, ]

color_v <- list(
  Color = c(Baseline_PRCR = "#8DD3C7", Baseline_PDSD = "#FFFFB3", C2 = "#BEBADA",
            EOT = "#FB8072")
)

library(pheatmap)
pdf("../Results/Overview2407/2407.myeloid.class.heatmap.hallmark.pdf", width = 8, height = 4)
pheatmap(as.matrix(finalres),
         scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         annotation_col = rowanno,
         annotation_colors = color_v,
         color =colorRampPalette(c("blue", "white","red"))(50),
         fontsize = 10
)
dev.off()



combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

hall <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol)

allscores <- split(hall$gene_symbol, hall$gs_name)
#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

expmtr <- combined@meta.data[, c("Class", names(allscores))]
expmtr <- aggregate(.~ Class, expmtr, median)
rownames(expmtr) <- expmtr$Class
expmtr <- t(expmtr[, -1])

rowanno <- data.frame(Color = c("Baseline_PRCR", "Baseline_PDSD", "C2", "EOT"))
rownames(rowanno) <- c("Baseline_PRCR", "Baseline_PDSD", "C2", "EOT")
expmtr <- expmtr[, c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT")]
colnames(expmtr) <- gsub("\\/", "", colnames(expmtr))

get_top5_nonredundant <- function(df) {
  df_long <- df %>% 
    pivot_longer(cols = -Pathway, names_to = "Category", values_to = "Value") %>% 
    group_by(Category) %>% 
    arrange(desc(Value)) %>% 
    slice_head(n = 5)
  
  df_wide <- df_long %>% 
    pivot_wider(names_from = Category, values_from = Value)
  
  df_wide <- df_wide %>% select(Pathway, all_of(names(df)[-which(names(df) == "Pathway")]))
  return(df_wide)
}
data <- as.data.frame(t(scale(t(expmtr))))
data$Pathway <- rownames(data)
write.table(data, "../Results/Overview2407/2407.myeloid.class.kegg.zscore.txt", row.names = F, col.names = T, quote = F, sep = "\t")

data$Pathway <- rownames(data)
# Get the top 2 non-redundant pathways
top5_pathways <- get_top5_nonredundant(data)
# Print the result
allpathes <- c()
for (col in colnames(data)[-which(colnames(data) == "Pathway")]) {
  cat("Top 5 pathways for", col, ":\n")
  selected_pathways <- top5_pathways %>% filter(!is.na(.data[[col]])) %>% select(Pathway, all_of(col))
  print(selected_pathways)
  cat("\n")
  allpathes <- c(allpathes, selected_pathways$Pathway)
}


finalres <- expmtr[allpathes, ]

color_v <- list(
  Color = c(Baseline_PRCR = "#8DD3C7", Baseline_PDSD = "#FFFFB3", C2 = "#BEBADA",
            EOT = "#FB8072")
)

library(pheatmap)
pdf("../Results/Overview2407/2407.myeloid.class.heatmap.kegg.pdf", width = 10, height = 4)
pheatmap(as.matrix(finalres),
         scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         annotation_col = rowanno,
         annotation_colors = color_v,
         color =colorRampPalette(c("blue", "white","red"))(50),
         fontsize = 10
)
dev.off()




combined <- readRDS("../4.Analy/04.myeloid/Seled.Myeloid.reclassified.reAnnot.subcells.TAM.rds")

hall <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, gene_symbol)

allscores <- split(hall$gene_symbol, hall$gs_name)
#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

expmtr <- combined@meta.data[, c("Class", names(allscores))]
expmtr <- aggregate(.~ Class, expmtr, median)
rownames(expmtr) <- expmtr$Class
expmtr <- t(expmtr[, -1])

rowanno <- data.frame(Color = c("Baseline_PRCR", "Baseline_PDSD", "C2", "EOT"))
rownames(rowanno) <- c("Baseline_PRCR", "Baseline_PDSD", "C2", "EOT")
expmtr <- expmtr[, c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT")]
colnames(expmtr) <- gsub("\\/", "", colnames(expmtr))

get_top5_nonredundant <- function(df) {
  df_long <- df %>% 
    pivot_longer(cols = -Pathway, names_to = "Category", values_to = "Value") %>% 
    group_by(Category) %>% 
    arrange(desc(Value)) %>% 
    slice_head(n = 5)
  
  df_wide <- df_long %>% 
    pivot_wider(names_from = Category, values_from = Value)
  
  df_wide <- df_wide %>% select(Pathway, all_of(names(df)[-which(names(df) == "Pathway")]))
  return(df_wide)
}
data <- as.data.frame(t(scale(t(expmtr))))
data$Pathway <- rownames(data)
write.table(data, "../Results/Overview2407/2407.myeloid.class.reactome.zscore.txt", row.names = F, col.names = T, quote = F, sep = "\t")

data$Pathway <- rownames(data)
# Get the top 2 non-redundant pathways
top5_pathways <- get_top5_nonredundant(data)
# Print the result
allpathes <- c()
for (col in colnames(data)[-which(colnames(data) == "Pathway")]) {
  cat("Top 5 pathways for", col, ":\n")
  selected_pathways <- top5_pathways %>% filter(!is.na(.data[[col]])) %>% select(Pathway, all_of(col))
  print(selected_pathways)
  cat("\n")
  allpathes <- c(allpathes, selected_pathways$Pathway)
}


finalres <- expmtr[allpathes, ]

color_v <- list(
  Color = c(Baseline_PRCR = "#8DD3C7", Baseline_PDSD = "#FFFFB3", C2 = "#BEBADA",
            EOT = "#FB8072")
)

library(pheatmap)
pdf("../Results/Overview2407/2407.myeloid.class.heatmap.reactome.pdf", width = 12, height = 4)
pheatmap(as.matrix(finalres),
         scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         annotation_col = rowanno,
         annotation_colors = color_v,
         color =colorRampPalette(c("blue", "white","red"))(50),
         fontsize = 10
)
dev.off()






