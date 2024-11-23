library(dplyr)
library(Seurat)
library(patchwork)
library(SingleR)
library(celldex)
library(plotly)
library(RColorBrewer)
library(htmlwidgets)

genes <- c("Pgam5")

sc_RNAseq <- function(srt){
  srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(srt), 10)
  
  all.genes <- rownames(srt)
  srt <- ScaleData(srt, features = all.genes)
  
  srt <- RunPCA(srt, features = VariableFeatures(object = srt))
  srt <- FindNeighbors(srt, dims = 1:10)
  srt <- FindClusters(srt, resolution = 0.5)
  
  srt <- RunUMAP(srt, dims = 1:10)
  umap_unlabel <- DimPlot(srt, reduction = "umap")
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  srt.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
  cluster0.markers <- FindMarkers(srt, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
  VlnPlot(srt, features = genes) + 
    theme(axis.text.y = element_text(size = 10))
  
  srt.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  
  heat <- DoHeatmap(srt, features = top10$gene) + NoLegend()
  
  return(list(top10=top10, srt=srt, heat=heat, umap_unlabel=umap_unlabel))
}

create_plots <- function(seurat_obj){
  
  exprs <- seurat_obj@assays$RNA@counts
  singler.results <- SingleR(test = exprs, ref = hpca.se, assay.type.test=1,
                             labels = hpca.se$label.main)
  
  seurat_obj <- AddMetaData(seurat_obj, metadata = singler.results$labels, col.name = "SingleR")

  gene_data <- sapply(genes, function(gene) seurat_obj[["RNA"]]@data[gene, ])
  umap.data <- cbind(FetchData(seurat_obj, vars = "SingleR"), gene_data)
  umap <- UMAPPlot(seurat_obj, group.by = "SingleR")
  umap <- ggplotly(umap)
  
  feature_plots <- lapply(genes, function(gene) FeaturePlot(seurat_obj, features = gene))
  vln <- VlnPlot(object = seurat_obj, features = genes, group.by = "SingleR")
  
  return(list(umap=umap, feature_plots=feature_plots, srt=seurat_obj, vln=vln))
}

#hpca.se <- HumanPrimaryCellAtlasData()
hpca.se <- MouseRNAseqData()
path <- "/Users/matthewpace/scRNAseq_data/mouse/GSE211584/raw/"
directories <- c("ACLR_28d_1", "ACLR_28d_2", "ACLR_7d_1", "ACLR_7d_2", 
                 "Sham_28d", "Sham_7d")
datasets <- list()
results <- list()
plots <- list()

# Load datasets
for (directory in directories) {
    datasets[[directory]] <- Read10X(data.dir = paste0(path, directory, "/"))
    srt_list[[directory]] <- CreateSeuratObject(counts = srt.data, project = "srt3k",
                                                min.cells = 3, min.features = 200)
      srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-")

  srt <- subset(srt, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
}

# Merge datasets
healthy_list <- srt_list[1:3]
OA_list <- srt_list[4:6]
healthy_data <- merge(healthy_list[[1]], y = healthy_list[-1],
                      add.cell.ids = c("Healthy1", "Healthy2", "Healthy3"))
OA_data <- merge(OA_list[[1]], y = OA_list[-1], add.cell.ids = c("OA1", "OA2", "OA3"))

# Process RNAseq
for (data_name in names(datasets)) {
  results[[paste0(data_name, "_results")]] <- sc_RNAseq(datasets[[data_name]])
}

# Create plots
for (result_name in names(results)) {
  plots[[paste0(result_name)]] <- create_plots(results[[result_name]]$srt)
}

# Save results
for (result_name in names(results)) {
  print(result_name)
  directory <- paste0("/Users/matthewpace/Desktop/", result_name)
  dir.create(directory)
  
  png_file <- paste0(directory, "/", result_name, "_heat.png")
  png(png_file)
  plot(results[[result_name]]$heat)
  dev.off()
  
  png_file <- paste0(directory, "/", result_name, "_umap.png")
  png(png_file)
  plot(results[[result_name]]$umap_unlabel)
  dev.off()
  
  csv_file <- paste0(directory, "/top10.csv")
  write.csv(results[[result_name]]$top10, csv_file, row.names = FALSE)
  
  htmlwidgets::saveWidget(plots[[result_name]]$umap, paste0(directory, "/", result_name, "_umap.html"))
  
  for (i in 1:length(genes)) {
    gene <- genes[i]
    png_file <- paste0(directory, "/", result_name, "_", gene, ".png")
    png(png_file)
    plot(plots[[result_name]]$feature_plots[[i]])
    dev.off()
  }
  
  png_file <- paste0(directory, "/", result_name, "_violin.png")
  png(png_file, width = 1400, height = 950)
  plot(plots[[result_name]]$vln)
  dev.off()
}

