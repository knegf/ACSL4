library(sp)
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(viridisLite)
library(viridis)
library(dplyr)
library(openxlsx)
library(tibble)


## GSE123813 ---------------------------------------------------------------------
### Data preparing ---------------------------------------------------------------------
SCC <- read.delim("scc_scRNA_counts.txt")
SCC_metadata <- read.delim("scc_metadata.txt")
barcode <- paste(sapply(strsplit(sapply(strsplit(SCC_metadata$cell.id, split = "_"), 
                                        function(x) x[1]), split = "[.]"), 
                        function(x) paste(x[2], x[3], sep = "_")), 
                 sapply(strsplit(SCC_metadata$cell.id, split = "_"), function(x) x[2]), 
                 sep = "_")
SCC_metadata$barcode <- barcode
SCC_metadata <- SCC_metadata[!(SCC_metadata$barcode %in% names(which(table(barcode) > 1))),]
SCC <- SCC[,SCC_metadata$cell.id]
colnames(SCC) <- SCC_metadata$barcode
rownames(SCC_metadata) <- SCC_metadata$barcode
SCC_metadata <- SCC_metadata[,c(1,7,2:6)]
SCC <- CreateSeuratObject(SCC, meta.data = SCC_metadata, project = "GSE123813_SCC")

SCC[["percent.mt"]] <- PercentageFeatureSet(SCC, pattern = "^MT-")
SCC[["percent.ribo"]] <- PercentageFeatureSet(SCC, pattern = "^RP[SL]")
SCC[["percent.hb"]] <- PercentageFeatureSet(SCC, pattern = "^HB[^(P)]")
message("Mitochondrial genes:")
rownames(SCC)[grep("^MT-",rownames(SCC),ignore.case = T)]
message("Ribosomal genes:")
rownames(SCC)[grep("^RP[SL]",rownames(SCC),ignore.case = T)]
message("Haemoglobin genes:")
rownames(SCC)[grep("^HB[^(P)]",rownames(SCC),ignore.case = T)]

SCC <- CellCycleScoring(SCC, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

### Data integrating ---------------------------------------------------------------------
SCC.list <- SplitObject(SCC, split.by = "patient")
SCC.list <- lapply(SCC.list, function(x){
  SCTransform(x, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), 
              return.only.var.genes = FALSE, seed.use = 2022, verbose = T)
})
features <- SelectIntegrationFeatures(object.list = SCC.list, nfeatures = 3000)
SCC.list <- PrepSCTIntegration(object.list = SCC.list, anchor.features = features)
integrate_anchors <- FindIntegrationAnchors(object.list = SCC.list, normalization.method = "SCT", anchor.features = features)
seu <- IntegrateData(anchorset = integrate_anchors, new.assay.name = "Integrated", normalization.method = "SCT")

### Dimensional reduction ---------------------------------------------------------------------
seu <- RunPCA(seu, features = features, seed.use = 2022, ndims.print = 1:2)
seu <- RunTSNE(seu, reduction = "pca", dims = 1:30, n.neighbors = 30, seed.use = 2022)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, n.neighbors = 30, seed.use = 2022)

### Clustering ---------------------------------------------------------------------
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30, k.param = 30, prune.SNN = 0.05)
clusters_louvain <- FindClusters(seu@graphs$Integrated_snn, algorithm = 2, 
                                 resolution = seq(0, 0.5, 0.05), random.seed = 2022)
seu$seurat_clusters <- clusters_louvain$res.0.3

### Cluster DE genes ---------------------------------------------------------------------
Idents(seu) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(seu, test.use = "wilcox", only.pos = T, random.seed = 2022)
message("Before significant selection:")
table(cluster_markers$cluster)
cluster.markers.onlypos.FDR0.05 <- subset(cluster_markers, cluster_markers$p_val_adj < 0.05)
message("After filtering")
table(cluster.markers.onlypos.FDR0.05$cluster)
cluster.markers.onlypos.FDR.top100 <- cluster.markers.onlypos.FDR0.05 %>% 
  group_by(cluster) %>% arrange(., -avg_log2FC, p_val_adj, .by_group = T) %>% slice_head(., n = 100)
write.xlsx(cluster.markers.onlypos.FDR.top100, 
           "SCC_cluster_markers_pos_FDR0.05_100topbyavglog2FC.xlsx", 
           append = F, rowNames = F, colNames = T, overwrite = T)

### Annotation ---------------------------------------------------------------------
DotPlot(seu, features = list(`CD4+ Tcells` = c("CD4","IL7R","CD69","ZFP36L2"), 
                             `CD8+ Tcells` = c("CD8A","CD8B","GZMK","GZMA"), 
                             `Treg`        = c("FOXP3","IL32","IL2RA")), 
        group.by = "seurat_clusters", col.min = 0, col.max = 0.75, dot.scale = 8) + 
  ggtitle("Feature dotplot for each cell type") + ylab("Clusters") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_blank(), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.45), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 11, vjust = 1), 
        legend.position = "top", legend.justification = "center") 
FeaturePlot(seu, reduction = "umap", features = c("CD4","IL7R","CD8A","GZMK"), ncol = 2, 
            coord.fixed = T, order = T, keep.scale = "all", min.cutoff = 0, max.cutoff = 4, cols = cividis(2))
FeaturePlot(seu, reduction = "tsne", features = c("CD4","IL7R","CD8A","GZMK"), ncol = 2, 
            coord.fixed = T, order = T, keep.scale = "all", min.cutoff = 0, max.cutoff = 4, cols = cividis(2))

### CD4 T cells treatment DE genes ---------------------------------------------------------------------
CD4Tcells <- subset(seu, seurat_clusters %in% c(1,4))
CD4Tcells_markers <- FindMarkers(CD4Tcells, assay = "RNA", slot = "counts", 
                                 group.by = "treatment", ident.1 = "post", ident.2 = "pre", 
                                 logfc.threshold = 0, test.use = "wilcox", min.pct = 0) %>% rownames_to_column("gene")
average_expression_all <- AverageExpression(CD4Tcells, assay = "RNA", slot = "counts", group.by = "treatment")$RNA %>% 
  as.data.frame() %>% `colnames<-`(paste0("avg_expr_all_",colnames(.))) %>% rownames_to_column("gene")
average_expression_part <- sapply(c("post","pre"), function(y){
  apply(GetAssayData(subset(CD4Tcells, treatment == y), assay = "RNA", slot = "counts"), MARGIN = 1, 
        function(x){ sum(x)/sum(x>0) })
}) %>% as.data.frame() %>% `colnames<-`(paste0("avg_expr_part_",colnames(.))) %>% rownames_to_column("gene")
average_expression_part[is.na(average_expression_part)] <- 0
average_expression_total <- rowMeans(GetAssayData(CD4Tcells, assay = "RNA", slot = "counts")) %>% as.data.frame() %>% 
  `colnames<-`("avg_expr_all") %>% rownames_to_column("gene")
average_expression_section <- apply(GetAssayData(CD4Tcells, assay = "RNA", slot = "counts"), MARGIN = 1, 
                                    function(x){ sum(x)/sum(x>0) }) %>% 
  as.data.frame() %>% `colnames<-`("avg_expr_part") %>% rownames_to_column("gene")
average_expression_section[is.na(average_expression_section)] <- 0
CD4Tcells_markers <- Reduce(function(x, y) merge(x, y, by="gene"), 
                            list(CD4Tcells_markers, average_expression_all, average_expression_total, 
                                 average_expression_part, average_expression_section))
write.xlsx(CD4Tcells_markers, 
           "SCC_CD4Tcell_treatment_markers_all_info.xlsx", 
           append = F, rowNames = F, colNames = T, overwrite = T)

### CD8 T cells treatment DE genes ---------------------------------------------------------------------
CD8Tcells <- subset(seu, seurat_clusters %in% c(2,3,8))
CD8Tcells_markers <- FindMarkers(CD8Tcells, assay = "RNA", slot = "counts", 
                                 group.by = "treatment", ident.1 = "post", ident.2 = "pre", 
                                 logfc.threshold = 0, test.use = "wilcox", min.pct = 0) %>% rownames_to_column("gene")
average_expression_all <- AverageExpression(CD8Tcells, assay = "RNA", slot = "counts", group.by = "treatment")$RNA %>% 
  as.data.frame() %>% `colnames<-`(paste0("avg_expr_all_",colnames(.))) %>% rownames_to_column("gene")
average_expression_part <- sapply(c("post","pre"), function(y){
  apply(GetAssayData(subset(CD8Tcells, treatment == y), assay = "RNA", slot = "counts"), MARGIN = 1, 
        function(x){ sum(x)/sum(x>0) })
}) %>% as.data.frame() %>% `colnames<-`(paste0("avg_expr_part_",colnames(.))) %>% rownames_to_column("gene")
average_expression_part[is.na(average_expression_part)] <- 0
average_expression_total <- rowMeans(GetAssayData(CD8Tcells, assay = "RNA", slot = "counts")) %>% as.data.frame() %>% 
  `colnames<-`("avg_expr_all") %>% rownames_to_column("gene")
average_expression_section <- apply(GetAssayData(CD8Tcells, assay = "RNA", slot = "counts"), MARGIN = 1, 
                                    function(x){ sum(x)/sum(x>0) }) %>% 
  as.data.frame() %>% `colnames<-`("avg_expr_part") %>% rownames_to_column("gene")
average_expression_section[is.na(average_expression_section)] <- 0
CD8Tcells_markers <- Reduce(function(x, y) merge(x, y, by="gene"), 
                            list(CD8Tcells_markers, average_expression_all, average_expression_total, 
                                 average_expression_part, average_expression_section))
write.xlsx(CD8Tcells_markers, 
           "SCC_CD8Tcell_treatment_markers_all_info.xlsx", 
           append = F, rowNames = F, colNames = T, overwrite = T)



## GSE154795 ---------------------------------------------------------------------
### Data preparing ---------------------------------------------------------------------
GBM <- readRDS("GSE154795_GBM.AllCell.Integrated.Scaled.ClusterRes.0.1.rds")
GBM <- UpdateSeuratObject(GBM)
TILs <- subset(GBM, idents = 4) # based one the publication, Cluster4 are the lympheloid

### Dimensional reduction ---------------------------------------------------------------------
TILs <- RunPCA(TILs, seed.use = 2022, ndims.print = 1:2)
TILs <- RunTSNE(TILs, reduction = "pca", dims = 1:20, n.neighbors = 20, seed.use = 2022)
TILs <- RunUMAP(TILs, reduction = "pca", dims = 1:20, n.neighbors = 20, seed.use = 2022)

### Clustering ---------------------------------------------------------------------
TILs <- FindNeighbors(TILs, reduction = "pca", dims = 1:20, k.param = 20, prune.SNN = 0.05)
clusters_louvain <- FindClusters(TILs@graphs$integrated_snn, algorithm = 2, 
                                 resolution = seq(0, 0.5, 0.05), random.seed = 2022)
TILs$seurat_clusters <- clusters_louvain$res.0.3
DimPlot(TILs, reduction = "umap", group.by = "seurat_clusters", label=T, label.size=6)+NoLegend()
DimPlot(TILs, reduction = "tsne", group.by = "seurat_clusters", label=T, label.size=6)+NoLegend()

### Cluster DE genes ---------------------------------------------------------------------
Idents(TILs) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(TILs, test.use = "wilcox", only.pos = T, random.seed = 2022)
message("Before significant selection:")
table(cluster_markers$cluster)
cluster.markers.onlypos.FDR0.05 <- subset(cluster_markers, cluster_markers$p_val_adj < 0.05)
message("After filtering")
table(cluster.markers.onlypos.FDR0.05$cluster)
cluster.markers.onlypos.FDR.top100 <- cluster.markers.onlypos.FDR0.05 %>% 
  group_by(cluster) %>% arrange(., -avg_log2FC, p_val_adj, .by_group = T) %>% slice_head(., n = 100)
write.xlsx(cluster.markers.onlypos.FDR.top100, 
           "GBM_TILs_cluster_markers_pos_FDR0.05_100topbyavglog2FC.xlsx", 
           append = F, rowNames = F, colNames = T, overwrite = T)

### Annotation ---------------------------------------------------------------------
DotPlot(TILs, features = list(`CD4+ Tcells` = c("CD4","IL7R","CD69","ZFP36L2"), 
                              `CD8+ Tcells` = c("CD8A","CD8B","GZMK","GZMA"), 
                              `Treg`        = c("FOXP3","IL32","IL2RA")), 
        assay = "RNA", group.by = "seurat_clusters", col.min = 0, col.max = 1.5, dot.scale = 8) + 
  ggtitle("Feature dotplot for each cell type") + ylab("Clusters") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_blank(), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.45), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 11, vjust = 1), 
        legend.position = "top", legend.justification = "center") 
FeaturePlot(seu, reduction = "umap", features = c("CD4","IL7R","CD8A","GZMK"), ncol = 2, 
            coord.fixed = T, order = T, keep.scale = "all", min.cutoff = 0, max.cutoff = 4, cols = cividis(2))
FeaturePlot(seu, reduction = "tsne", features = c("CD4","IL7R","CD8A","GZMK"), ncol = 2, 
            coord.fixed = T, order = T, keep.scale = "all", min.cutoff = 0, max.cutoff = 4, cols = cividis(2))

### CD4 T cells treatment DE genes ---------------------------------------------------------------------
CD4Tcells <- subset(TILs, seurat_clusters %in% c(2,5,8))
CD4Tcells$treatment <- gsub("GBM.PD1","treatment",gsub("GBM.rec","nontreatment",gsub("GBM.new","nontreatment",CD4Tcells$condition)))
CD4Tcells_markers <- FindMarkers(CD4Tcells, assay = "RNA", slot = "counts", 
                                 group.by = "treatment", ident.1 = "treatment", ident.2 = "nontreatment", 
                                 logfc.threshold = 0, test.use = "wilcox", min.pct = 0) %>% rownames_to_column("gene")
average_expression_all <- AverageExpression(CD4Tcells, assay = "RNA", slot = "counts", group.by = "treatment")$RNA %>% 
  as.data.frame() %>% `colnames<-`(paste0("avg_expr_all_",colnames(.))) %>% rownames_to_column("gene")
average_expression_part <- sapply(c("treatment","nontreatment"), function(y){
  apply(GetAssayData(subset(CD4Tcells, treatment == y), assay = "RNA", slot = "counts"), MARGIN = 1, 
        function(x){ sum(x)/sum(x>0) })
}) %>% as.data.frame() %>% `colnames<-`(paste0("avg_expr_part_",colnames(.))) %>% rownames_to_column("gene")
average_expression_part[is.na(average_expression_part)] <- 0
average_expression_total <- rowMeans(GetAssayData(CD4Tcells, assay = "RNA", slot = "counts")) %>% as.data.frame() %>% 
  `colnames<-`("avg_expr_all") %>% rownames_to_column("gene")
average_expression_section <- apply(GetAssayData(CD4Tcells, assay = "RNA", slot = "counts"), MARGIN = 1, 
                                    function(x){ sum(x)/sum(x>0) }) %>% 
  as.data.frame() %>% `colnames<-`("avg_expr_part") %>% rownames_to_column("gene")
average_expression_section[is.na(average_expression_section)] <- 0
CD4Tcells_markers <- Reduce(function(x, y) merge(x, y, by="gene"), 
                            list(CD4Tcells_markers, average_expression_all, average_expression_total, 
                                 average_expression_part, average_expression_section))
write.xlsx(CD4Tcells_markers, 
           "GBM_CD4Tcell_treatment_markers_all_info.xlsx", 
           append = F, rowNames = F, colNames = T, overwrite = T)

### CD8 T cells treatment DE genes ---------------------------------------------------------------------
CD8Tcells <- subset(TILs, seurat_clusters %in% c(1,6))
CD8Tcells$treatment <- gsub("GBM.PD1","treatment",gsub("GBM.rec","nontreatment",gsub("GBM.new","nontreatment",CD8Tcells$condition)))
CD8Tcells_markers <- FindMarkers(CD8Tcells, assay = "RNA", slot = "counts", 
                                 group.by = "treatment", ident.1 = "treatment", ident.2 = "nontreatment", 
                                 logfc.threshold = 0, test.use = "wilcox", min.pct = 0) %>% rownames_to_column("gene")
average_expression_all <- AverageExpression(CD8Tcells, assay = "RNA", slot = "counts", group.by = "treatment")$RNA %>% 
  as.data.frame() %>% `colnames<-`(paste0("avg_expr_all_",colnames(.))) %>% rownames_to_column("gene")
average_expression_part <- sapply(c("treatment","nontreatment"), function(y){
  apply(GetAssayData(subset(CD8Tcells, treatment == y), assay = "RNA", slot = "counts"), MARGIN = 1, 
        function(x){ sum(x)/sum(x>0) })
}) %>% as.data.frame() %>% `colnames<-`(paste0("avg_expr_part_",colnames(.))) %>% rownames_to_column("gene")
average_expression_part[is.na(average_expression_part)] <- 0
average_expression_total <- rowMeans(GetAssayData(CD8Tcells, assay = "RNA", slot = "counts")) %>% as.data.frame() %>% 
  `colnames<-`("avg_expr_all") %>% rownames_to_column("gene")
average_expression_section <- apply(GetAssayData(CD8Tcells, assay = "RNA", slot = "counts"), MARGIN = 1, 
                                    function(x){ sum(x)/sum(x>0) }) %>% 
  as.data.frame() %>% `colnames<-`("avg_expr_part") %>% rownames_to_column("gene")
average_expression_section[is.na(average_expression_section)] <- 0
CD8Tcells_markers <- Reduce(function(x, y) merge(x, y, by="gene"), 
                            list(CD8Tcells_markers, average_expression_all, average_expression_total, 
                                 average_expression_part, average_expression_section))
write.xlsx(CD8Tcells_markers, 
           "GBM_CD8Tcell_treatment_markers_all_info.xlsx", 
           append = F, rowNames = F, colNames = T, overwrite = T)







