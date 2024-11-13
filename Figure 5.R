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


############################Figure 5############################
##############Figure 5a_f##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")
allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$MinorCluster))]
pdf("../Results/Overview2401/2401.tcell.cd4t.overview.pdf", width = 5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "MinorCluster", cols = color_v)+ggtitle("T cells"))
dev.off()

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MinorCluster", "UMAP1", "UMAP2")], "../Results/RawData/figure5a.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")
allcolors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))
color_v=allcolors[1:length(unique(combined$MinorCluster))]
pdf("../Results/Overview2401/2401.tcell.cd8t.overview.pdf", width = 5, height = 3)
print(DimPlot(combined, reduction = "umap", label = FALSE, group.by = "MinorCluster", cols = color_v)+ggtitle("T cells"))
dev.off()

combined$UMAP1 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 1])
combined$UMAP2 <- as.numeric(Embeddings(object = combined, reduction = "umap")[, 2])
write.table(combined@meta.data[, c("ID", "MinorCluster", "UMAP1", "UMAP2")], "../Results/RawData/figure5f.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

##############Figure 5b_g##############
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(dplyr)
library(plyr)
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol)
kegg <- kegg[kegg$gs_name == "KEGG_CELL_CYCLE", ]
names(kegg) <- c("Signature", "Genes")

combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")
metapros <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd4t.signature.txt", header = T, stringsAsFactors = F, data.table = F)
metapros <- rbind(metapros, kegg)

allscores <- split(metapros$Genes, metapros$Signature)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

expmtr <- combined@meta.data[, c(35, 73:90)]
expmtr <- aggregate(.~ MinorCluster, expmtr, median)
rownames(expmtr) <- expmtr$MinorCluster

write.table(expmtr, "../Results/RawData/figure5b.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


expmtr <- expmtr[, -1]


sortsignatures <- c("Na\303\257ve", "Activation/Effector function",   ##Differentiation
                    "TCR signaling", "Cytotoxicity", "Cytokine/Cytokine receptor", 
                    "Chemokine/Chemokine receptor", "Stress response", "Adhesion", "IFN response", 
                    "Treg signature", "Costimulatory molecules",    ##Function
                    "OXPHOS", "Glycolysis", "Lipid metabolism",    ##Metabolism
                    "Pro-apoptosis", "Anti-apoptosis",    ##Apoptosis
                    "KEGG_CELL_CYCLE"
)
expmtr <- expmtr[, sortsignatures]

library(pheatmap)
pdf("../Results/Overview2401/2401.tcell.cd4t.signature.pdf", width = 5, height = 4)
pheatmap(as.matrix(t(expmtr)),
         scale = "row",
         cluster_rows = F,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         color =colorRampPalette(c("blue", "white","red"))(50),
         fontsize = 10
)
dev.off()



combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")
metapros <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cd8t.signature.txt", header = T, stringsAsFactors = F, data.table = F)
metapros <- rbind(metapros, kegg)

allscores <- split(metapros$Genes, metapros$Signature)

#"Naive" "TCR Signaling" "Activation:Effector function" "Exhaustion" "Cytotoxicity"
combined <- AddModuleScore(combined,
                           features = allscores,
                           name="TCELLSIGNATURE")
colnames(combined@meta.data)[grep("TCELLSIGNATURE", colnames(combined@meta.data))] <- names(allscores)

expmtr <- combined@meta.data[, c(35, 73:91)]
expmtr <- aggregate(.~ MinorCluster, expmtr, median)
rownames(expmtr) <- expmtr$MinorCluster
write.table(expmtr, "../Results/RawData/figure5g.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)
expmtr <- expmtr[, -1]

sortsignatures <- c("Na\303\257ve", "Activation:Effector function", "Exhaustion",   ##Differentiation
                    "TCR Signaling", "Cytotoxicity", "Cytokine/Cytokine receptor", 
                    "Chemokine/Chemokine receptor", "Senescence", "Anergy", "Stress response", 
                    "IFN Response", "NFKB Signaling", "MAPK Signaling",    ##Function
                    "Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism",    ##Metabolism
                    "Pro-apoptosis", "Anti-apoptosis",    ##Apoptosis
                    "KEGG_CELL_CYCLE"
)
expmtr <- expmtr[, sortsignatures]

library(pheatmap)
pdf("../Results/Overview2401/2401.tcell.cd8t.signature.pdf", width = 5, height = 5)
pheatmap(as.matrix(t(expmtr)),
         scale = "row",
         cluster_rows = F,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         color =colorRampPalette(c("blue", "white","red"))(50),
         fontsize = 10
)
dev.off()


##############Figure 5c_h##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")
#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Overview2401/2401.tcell.cd4t."

OR.all.list <- do.tissueDist(cellInfo.tb=meta.tb,
                             meta.cluster = meta.tb$Class,
                             colname.patient = "MiratiID",
                             loc = meta.tb$MinorCluster,
                             out.prefix=sprintf("%s.Class",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../Results/Overview2401/2401.tcell.cd4t.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Overview2401/2401.tcell.cd4t.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Overview2401/2401.tcell.cd4t.Class.diff.pdf", p1, width = 5, height = 4)

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

pdf("../Results/Overview2401/2401.tcell.cd4t.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()



combined <- readRDS("../4.Analy/03.t/CD8T.reclassified.reannot.cells.rds")
#combined$CCluster <- factor(combined$CCluster, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))
Idents(combined) <- combined$MinorCluster

library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- combined@meta.data
out.prefix <- "../Results/Overview2401/2401.tcell.cd8t."

OR.all.list <- do.tissueDist(cellInfo.tb=meta.tb,
                             meta.cluster = meta.tb$Class,
                             colname.patient = "MiratiID",
                             loc = meta.tb$MinorCluster,
                             out.prefix=sprintf("%s.Class",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../Results/Overview2401/2401.tcell.cd8t.Class.diff.txt", 
       row.names = T, quote = F)

cellfrac <- fread("../Results/Overview2401/2401.tcell.cd8t.Class.diff.txt", header = T, stringsAsFactors = F, data.table = F)
cellfrac <- melt(cellfrac, id.vars = "V1")
cellfrac$value <- ifelse(cellfrac$value > 3, 3, cellfrac$value)
cellfrac$V1 <- factor(x = cellfrac$V1, levels = c("Baseline_PR/CR", "Baseline_PD/SD", "C2", "EOT"))

p1 <- ggplot(cellfrac, aes(x = V1, y = variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7C93C3", mid = "white", high = "#750E21", midpoint = 1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

ggsave("../Results/Overview2401/2401.tcell.cd8t.Class.diff.pdf", p1, width = 5, height = 4)
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

pdf("../Results/Overview2401/2401.tcell.cd8t.Class.diff.bar.pdf", width = 5, height = 3.5)
p
dev.off()



##############Figure 5d##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")
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
pdf("../Results/Overview2401/2401.tcell.cd4t.umap.desity.pdf", width = 12, height = 3)
print(g)
dev.off()
g <- ggplot(data = meta, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = MinorCluster), size = 0.2) +
  scale_color_manual(values = color_v) +
  facet_wrap(~Class,ncol = 4) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.position = "none")
pdf("../Results/Overview2401/2401.tcell.cd4t.umap.group.pdf", width = 14, height = 3)
print(g)
dev.off()


##############Figure 5e_i##############
combined <- readRDS("../4.Analy/03.t/CD4T.reclassified.reannot.cells.rds")
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

write.table(data, "../Results/RawData/figure5e.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# plot with default options:
library(RColorBrewer)
allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
colors_border <- allcolors[1:4]
library(scales)
colors_in <- alpha(colors_border,0.5)

pdf("../Results/Overview2401/2401.tcell.cd4t.signature.radar.pdf", width = 5, height = 5)
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

write.table(data, "../Results/RawData/figure5i.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# plot with default options:
library(RColorBrewer)
allcolors <- c(RColorBrewer::brewer.pal(12, "Set3"))
colors_border <- allcolors[1:4]
library(scales)
colors_in <- alpha(colors_border,0.5)

pdf("../Results/Overview2401/2401.tcell.cd8t.signature.radar.pdf", width = 5, height = 5)
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




