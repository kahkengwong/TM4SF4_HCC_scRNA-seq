library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(viridis)
library(tidyverse)
library(data.table)

# Set file names, initialize a list to store Seurat objects, loop through each file, and merge all Seurat objects
file_names <- c("GSM3980137_Adult-Gall-Bladder1.txt", "GSM3980136_Adult-Gall-Bladder2.txt")
seurat_objects <- list()
for (i in seq_along(file_names)) {
  data <- fread(file_names[i])
  data <- as_tibble(data)
  data <- column_to_rownames(data, var = "GENE")
  data <- as.data.frame(data)
  rownames(data) <- make.names(rownames(data), unique = TRUE)
  data[] <- lapply(data, as.numeric)
  regular_matrix <- as.matrix(data)
  sparse_matrix <- Matrix(regular_matrix, sparse = TRUE)
  seurat_obj <- CreateSeuratObject(counts = sparse_matrix, project = paste0("GSM3980137_Adult-Gall-Bladder", i), min.cells = 3, min.features = 200)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nCount_RNA > 0)
  seurat_objects[[i]] <- seurat_obj
}

merged_seurat_obj <- merge(seurat_objects[[1]], y = seurat_objects[-1], project = "Merged")

# Mitochondrial gene content filtering
n_cells_before_mt <- ncol(merged_seurat_obj)
merged_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_obj, pattern = "^mt-")
max_mt_percent <- quantile(merged_seurat_obj@meta.data$percent.mt, 0.99)
merged_seurat_obj <- subset(merged_seurat_obj, subset = percent.mt < max_mt_percent)
n_cells_after_mt <- ncol(merged_seurat_obj)
n_excluded_mt <- n_cells_before_mt - n_cells_after_mt

# Normalize data and identify highly variable genes
merged_seurat_obj <- NormalizeData(merged_seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale data and regress out cell cycle effects
s_genes_default <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")
g2m_genes_default <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")
merged_seurat_obj <- CellCycleScoring(merged_seurat_obj, s.features = s_genes_default, g2m.features = g2m_genes_default, set.ident = TRUE)
merged_seurat_obj <- ScaleData(merged_seurat_obj, vars.to.regress = c("S.Score", "G2M.Score"))

# Doublets removal
sce_obj <- as.SingleCellExperiment(merged_seurat_obj)
samples <- merged_seurat_obj@meta.data$orig.ident
doublet_scores <- scDblFinder(sce_obj, samples = samples, k = 30, nfeatures = 2000)
batch_thresholds <- tapply(doublet_scores$scDblFinder.score, samples, function(x) quantile(x, probs = 0.95))
doublet_cells <- colnames(sce_obj[, which(doublet_scores$scDblFinder.score > batch_thresholds[samples])])
doublet_cells <- colnames(sce_obj)[mapply(function(x, y) x > batch_thresholds[y], doublet_scores$scDblFinder.score, samples)]
merged_seurat_obj$doublet <- colnames(merged_seurat_obj) %in% doublet_cells
merged_seurat_obj <- subset(merged_seurat_obj, subset = doublet == FALSE)
merged_seurat_obj$filtered <- "filtered"

# Batch effects correction
merged_seurat_obj_before_integration <- FindVariableFeatures(merged_seurat_obj, selection.method = "vst", nfeatures = 2000)
merged_seurat_obj_before_integration <- RunPCA(merged_seurat_obj, verbose = FALSE)
merged_seurat_obj_before_integration <- RunUMAP(merged_seurat_obj_before_integration, dims = 1:15, verbose = FALSE)
sample_list <- SplitObject(merged_seurat_obj_before_integration, split.by = "orig.ident")
for (i in 1:length(sample_list)) {
  sample_list[[i]] <- SCTransform(sample_list[[i]], verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = sample_list, dims = 1:7, verbose = FALSE) 
merged_seurat_obj_integrated <- IntegrateData(anchorset = anchors, dims = 1:7, verbose = FALSE)
merged_seurat_obj_integrated <- ScaleData(merged_seurat_obj_integrated, verbose = FALSE)
merged_seurat_obj_integrated <- RunPCA(merged_seurat_obj_integrated, verbose = FALSE)
merged_seurat_obj_integrated <- FindNeighbors(merged_seurat_obj_integrated, dims = 1:15)
merged_seurat_obj_integrated <- FindClusters(merged_seurat_obj_integrated, resolution = 0.5)
merged_seurat_obj_integrated <- RunUMAP(merged_seurat_obj_integrated, dims = 1:15, verbose = FALSE)

# Dimensionality reduction, find clusters, and UMAP plots
ElbowPlot(merged_seurat_obj_integrated)
merged_seurat_obj_integrated <- FindNeighbors(merged_seurat_obj_integrated, dims = 1:20)
merged_seurat_obj_integrated <- FindClusters(merged_seurat_obj_integrated, resolution = 1)
cluster_markers <- FindAllMarkers(merged_seurat_ob
umap_data <- as.data.frame(Embeddings(merged_seurat_obj_integrated, "umap"))
umap_data$cluster_id <- Idents(merged_seurat_obj_integrated)
umap_data_mean <- aggregate(. ~ cluster_id, data = umap_data, FUN = mean)
plasma_func <- colorRampPalette(viridis::viridis(100, direction = -1, option = "plasma"))
portion <- 0.8
n_colors <- round(length(unique(umap_data$cluster_id)) / portion)
plasma_colors <- plasma_func(n_colors)

umap_plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = as.factor(cluster_id))) +
    geom_point(size = 0.3, alpha = 0.95) +
    scale_color_manual(values = plasma_colors) + 
    geom_text(data = umap_data_mean, aes(label = cluster_id, x = UMAP_1, y = UMAP_2),
              color = "black", size = 3, fontface = "bold", check_overlap = TRUE) + 
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        panel.background = element_rect(fill = "white")
    ) +
    labs(title = "UMAP plot colored by cluster", x = "UMAP_1", y = "UMAP_2", color = "Cluster") +
    guides(color = guide_legend(override.aes = list(size = 3)))

print(umap_plot)

# Plot UMAP with TM4SF4 expression levels
gene_colors_alpha <- c(scales::alpha("lightgray", 0.85), scales::alpha("lightpink", 0.85), scales::alpha("#FF6666", 0.85), 
                       scales::alpha("#BC2727", 0.85), scales::alpha("#660000", 0.85))

FeaturePlot(merged_seurat_obj_integrated, features = "TM4SF4", min.cutoff = 'q10', max.cutoff = 'q90',
            pt.size = 0.2, cols = gene_colors_alpha) +
  theme(
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      panel.background = element_rect(fill = "white")
  ) +
  labs(title = "UMAP plot colored by TM4SF4 expression", x = "UMAP_1", y = "UMAP_2")