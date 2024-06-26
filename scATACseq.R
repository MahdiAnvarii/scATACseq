install.packages("Signac")
install.packages("EnsDb.Hsapiens.v75")
install.packages("biovizBase")
BiocManager::install("biovizBase")
library(Signac)
library(EnsDb.Hsapiens.v75)
library(Seurat)
library(tidyverse)
library(biovizBase)

# Reading data
fragments_file <- read.delim('../Data/atac_v1_pbmc_10k_fragments.tsv.gz', header = F, nrows = 10)
fragments_file

counts <- Read10X_h5('../Data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5')
counts[1:10,1:10]

chromatin_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":","-"),
  fragments = "../Data/atac_v1_pbmc_10k_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

str(chromatin_assay)

meta.data <- read.csv('../Data/atac_v1_pbmc_10k_singlecell.csv', header = T, row.names = 1)
view(meta.data)

Seurat_obj <- CreateSeuratObject(
  counts = chromatin_assay,
  meta.data = meta.data,
  assay = 'ATAC'
)

str(Seurat_obj)

# Adding gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
annotations
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
Annotation(Seurat_obj) <- annotations
Seurat_obj@assays$ATAC@annotation

# QC
# 1.Nucleosome signal and nucleosme binding pattern
# 2.TSS enrichment score
# 3.Total number of fragments in peaks
# 4.Fraction of fragments in peaks
# 5.Ratio reads in genomic blacklist regions
Seurat_obj <- NucleosomeSignal(Seurat_obj)
Seurat_obj <- TSSEnrichment(object = Seurat_obj, fast = F)
Seurat_obj$blacklist_ratio <- Seurat_obj$blacklist_region_fragments / Seurat_obj$peak_region_fragments
Seurat_obj$pct_reads_in_peaks <- Seurat_obj$peak_region_fragments / Seurat_obj$passed_filters * 100

view(Seurat_obj@meta.data)
colnames(Seurat_obj@meta.data)

DS1 <- DensityScatter(Seurat_obj, x='nCount_ATAC', y='TSS.enrichment' , quantiles = T, log_x = T)
ggsave("DS_nCount_TSS.pdf", plot = DS1, device = "pdf", height = 7, width = 10)

DS2 <- DensityScatter(Seurat_obj, x='nucleosome_signal', y='TSS.enrichment' , quantiles = T, log_x = T)
ggsave("DS_signal_TSS.pdf", plot = DS2, device = "pdf", height = 7, width = 10)

vln <- VlnPlot(Seurat_obj,
        features = c('nCount_ATAC', 'nFeature_ATAC' ,'TSS.enrichment','nucleosome_signal' , 'blacklist_ratio','pct_reads_in_peaks'),
        pt.size = 0.1, ncol = 6)
ggsave("vln_plot.pdf", plot = vln, device = "pdf", height = 7, width = 14)

# Filtering
dim(Seurat_obj)
Seurat_obj <- subset(x = Seurat_obj, 
                     subset = nCount_ATAC > 3000
                     & nCount_ATAC < 30000
                     & pct_reads_in_peaks > 15 
                     & nucleosome_signal < 4
                     & TSS.enrichment > 3
                     & blacklist_ratio < 0.05)
dim(Seurat_obj)

### Preprocess workflow
# Normalization
Seurat_obj <- RunTFIDF(Seurat_obj)

# Collect highly variable features
Seurat_obj <- FindTopFeatures(Seurat_obj, min.cutoff = 'q0')

# Linear dimensionality reduction
Seurat_obj <- RunSVD(Seurat_obj)
depth_cor <- DepthCor(Seurat_obj)
ggsave("DepthCor_plot.pdf", plot = depth_cor, device = "pdf", height = 7, width = 12)

# Non-Linear
Seurat_obj <- RunUMAP(Seurat_obj , dims = 2:30 , reduction = 'lsi')
Seurat_obj <- FindNeighbors(Seurat_obj , dims = 2:30 , reduction = 'lsi')
Seurat_obj <- FindClusters(Seurat_obj , algorithm = 3)
dimplot1 <- DimPlot(Seurat_obj, label = T) + NoLegend()
ggsave("clusters_dim_plot.pdf", plot = dimplot1, device = "pdf", width = 16 , height = 10)
