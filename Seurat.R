library(harmony)
library(Seurat)
library(dplyr)
library(ggsci)
library(ggplot2)
library(sctransform)
library(patchwork)
library(RColorBrewer)
library(SeuratObject)
library(infercnv)
library(tidyverse)
library(AnnoProbe)
library(tidydr)


#Data processing and quality control---

setwd('/data/scrna/')
sample_paths <- list(
  "Patient 1" = "/data/scrna/Patient 1",
  "Patient 2" = "/data/scrna/Patient 2",
  "Patient 3" = "/data/scrna/Patient 3",
  "Patient 4" = "/data/scrna/Patient 4",
  "Patient 5" = "/data/scrna/Patient 5"
)

seurat_list <- lapply(names(sample_paths), function(patient) {
  data <- Read10X(data.dir = sample_paths[[patient]])
  CreateSeuratObject(counts = data, project = patient, min.cells = 3, min.features = 200)
})

scrna <- merge(seurat_list[[1]], y = seurat_list[2:5])


scrna[["percent.mt"]] <- PercentageFeatureSet(object = scrna, pattern = "^MT-")
scrna[["percent.hb"]] <- PercentageFeatureSet(object = scrna, pattern = "^HB[^(P)]")
scrna[["percent.rp"]] <- PercentageFeatureSet(object = scrna, pattern = "^RP[SL]")


pdf('01_vlnplot_Pre_QC.pdf',width = 18,height = 6)
VlnPlot(object = scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb"), ncol = 4,pt.size = 0,
        group.by = 'orig.ident')
dev.off()
saveRDS(object = scrna ,file = "scrna_raw.rds")


scrna <- subset(scrna, 
                subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & 
                  percent.mt < 20 & 
                  percent.hb < 10 & 
                  nCount_RNA < quantile(nCount_RNA,0.95) & nCount_RNA > 1000)

pdf('01_vlnplot_After_QC.pdf',width = 18,height = 6)
VlnPlot(object = scrnaa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb"), ncol = 4,pt.size = 0,
        group.by = 'orig.ident')
dev.off()

saveRDS(object = scrna ,file = "scrna_filtered.rds")






#Standardization, dimensionality, reduction and clustering---

setwd('/data/scrna/01_seurat')


scrna<-readRDS('/data/scrna/scrna_filtered.rds')
scrna <- NormalizeData(scrna, normalization.method = "LogNormalize", scale.factor = 10000)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scrna)
scrna <- ScaleData(scrna, features = all.genes)


g2m_genes <- cc.genes$g2m.genes 
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(scrna))
s_genes <- cc.genes$s.genes 
s_genes <- CaseMatch(search=s_genes, match=rownames(scrna)) 
scrna[["RNA"]] <- JoinLayers(scrna[["RNA"]])
scrna <- CellCycleScoring(scrna, g2m.features=g2m_genes, s.features=s_genes)

scrna <- RunPCA(scrna, features = all.genes)

pdf('01_phase.pdf')
DimPlot(scrna,group.by = 'Phase')
dev.off()


pdf("01_scrna_ElbowPlot.pdf",width = 8,height = 6)
ElbowPlot(object = scrna,ndims = 50)
dev.off()



scrna <- RunHarmony(scrna, "orig.ident")
scrna <- RunUMAP(scrna, reduction = "harmony",dims = 1: 25)
scrna <- RunTSNE(scrna, reduction = "harmony",dims = 1: 25)

pdf("01_harmony_umap.pdf",width = 14,height = 12)
DimPlot(object = scrna, reduction = "umap",group.by = 'orig.ident')
dev.off()

pdf("01_harmony_tsne.pdf",width = 14,height = 12)
DimPlot(object = scrna, reduction = "tsne",group.by = 'orig.ident')
dev.off()


scrna<- FindNeighbors(scrna, reduction = "harmony", dims = 1:25)
scrna<- FindClusters(scrna, reduction = "harmony", resolution = 1.2)

pdf("01_umap_dimplot.pdf",width = 10,height = 8)
DimPlot(object = scrna, reduction = "umap",label = T)
dev.off()


scrna.markers <- FindAllMarkers(scrna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- scrna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(scrna.markers,file = "01_markers-seurat_clusters.csv")
write.csv(top10,file = "01_markers-top10.csv")
saveRDS(scrna,"scrna_not.rds")



pdf("01_umap_dimplot.pdf",width = 10,height = 8)
DimPlot(object = scrna, reduction = "umap",label = T)
dev.off()

pdf("01_tsne_dimplot.pdf",width = 10,height = 8)
DimPlot(object = scrna, reduction = "tsne",label = T)
dev.off()


table(scrna$seurat_clusters)
av <- AverageExpression(scrna,group.by = "seurat_clusters",assays = "RNA")
av=av[[1]]
av <- as.data.frame(av)
head(av)
cg=names(tail(sort(apply(av, 1, sd)),1000))
pdf('01_pheatmap_seurat_clusters.pdf',width = 10,height = 10)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()


#Cell population annotation---

cell_markers <- list(
  DC = c("GPR183", "GZMB"),
  B_cells_Plasma = c("IGKC", "JCHAIN"),
  B_cells = c("CD79A", "CD19", "BANK1",'MS4A1'),
  Mast_cells = c("TPSAB1", "TPSB2", "MS4A2"),
  Myeloid_cells = c("CD14", "CD68", "MS4A4A"),
  T_cells_NK_cells = c("CD3D", "CD3E", "CD2", "CD7", "NKG7"),
  Endothelial_cells = c("VWF", "A2M","PECAM1"),
  Epithelial_cells = c("EPCAM", "KRT18", "KRT19", "MUC1", "MUC5AC"),
  Fibroblasts = c("DCN", 'COL1A1','COL1A2'),
  Smooth_muscle_cells = c("ACTG2", "MYH11"),
  Macrophages = c('C1QA','C1QB')
  
)

str(cell_markers)


for (cell_type in names(cell_markers)) {
  genes <- cell_markers[[cell_type]]
  num_genes <- length(genes)
  ncol <- ceiling(num_genes / 2) 
  num_rows <- ceiling(num_genes / ncol)
  width <- ncol * 5 
  height <- num_rows * 5 
  p <- VlnPlot(scrna, features = genes, ncol = ncol, pt.size = 0)
  pdf_filename <- paste0("vlnplot_", cell_type, ".pdf")
  ggsave(filename = pdf_filename, plot = p, width = width, height = height, units = "in")
}


marker <-  unique(unlist(cell_markers))

pdf('01_cluster_marker_DotPlot.pdf', width = 16, height = 8)
DotPlot(scrna, features = marker, group.by = 'seurat_clusters') + 
  RotatedAxis() +
  scale_x_discrete("") + 
  scale_y_discrete("") +
  scale_color_gradientn(colors = c("#2b559c", "#f7f7f7", "#fdae61", "#a31d1f"), values = c(0, 0.4, 0.7, 1))
dev.off()


new.cluster.ids<-c("Fibroblast"	,
                   "NK&T cells"	,
                   "Epithelial cells"	,
                   "Fibroblast"	,
                   "Endothelial cells"	,
                   "Fibroblast"	,
                   "Smooth muscle cells"	,
                   "Epithelial cells"	,
                   "Fibroblast"	,
                   "Macrophages"	,
                   "Epithelial cells"	,
                   "Endothelial cells"	,
                   "Fibroblast"	,
                   "Fibroblast"	,
                   "Fibroblast"	,
                   "NK&T cells"	,
                   "Epithelial cells"	,
                   "Macrophages"	,
                   "B cells"	,
                   "Mast cells"	,
                   "NK&T cells"	,
                   "Epithelial cells"	,
                   "Fibroblast"	,
                   "Endothelial cells"	,
                   "Fibroblast"	,
                   "Epithelial cells"	,
                   "Smooth muscle cells"	)


names(new.cluster.ids)<-levels(scrna)
scrna<-RenameIdents(scrna,new.cluster.ids)
scrna[['cell_type']] <- scrna@active.ident
table(scrna$cell_type)


mycolors<- c("#2171B5", "#9ECAE1","#FF7f00","#FDD0A2",
             "#238443","#ADDD8E","#6A51A3","#BCBDDC",
             "#E41A1C")
pdf('01_Dimplot_umap_celltype2.pdf',width = 6,height = 6)
DimPlot(scrna,reduction = 'umap',label = F,cols = mycolors,group.by = "cell_type") +
  #scale_colour_npg() +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  #NoLegend() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()
        ,legend.text=element_text(size=12)
  )+labs(title = "")+
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()


pdf('01_Dimplot_tsne_celltype.pdf',width = 6,height = 6)
DimPlot(scrna,reduction = 'tsne',label = F,cols = mycolors,group.by = "cell_type") +
  #scale_colour_npg() +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  #NoLegend() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()
        ,legend.text=element_text(size=12)
  )+labs(title = "")+
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()


marker <- unique(c(
  'LUM',"DCN", 
  "CD7","CD2", "CD3D", "CD3E",  "NKG7",
  "EPCAM", "KRT18", "KRT19", "MUC1",
  'A2M',"VWF", "PECAM1",
  "MYH11", "ACTA2", 'TAGLN','RGS5',
  'C1QA','C1QB','C1QC',
  "IGKC", "CD79A","BANK1",'MS4A1',
  "TPSAB1", "TPSB2", "MS4A2"))


pdf('01_celltype_marker_DotPlot.pdf',width = 12,height = 3.7)
DotPlot(scrna, features = marker,group.by = 'cell_type')+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

saveRDS(scrna,"scrna_Annotated.rds")


