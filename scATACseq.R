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
