#!/usr/bin/env Rscript
library(Seurat)
library(scRank)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Inputs
seuratObj <- args[1]
column <- args[2]
species <- args[3]
target <- args[4]

if (seuratObj == 'AML_object.rda') {
    load(seuratObj)
    seuratObj <- seuratObj[c(VariableFeatures(seuratObj)[1:200], target),]
} else {
    seuratObj <- readRDS(seuratObj)
}

# Keep at max 2000 cells per cell type
downsampled_cells <- seuratObj@meta.data %>% tibble::rowid_to_column("cell_id") %>%
  group_by(!!sym(column)) %>%
  slice_sample(n = 250) %>%
  pull(cell_id)

ncells <- length(downsampled_cells)

seuratObj <- seuratObj[, downsampled_cells]

project <- seuratObj@project.name

# Run scRank
obj <- CreateScRank(input = seuratObj,
                    species = species, 
                    cell_type = column,
                    target = target)

obj@para$gene4use <- unique(c(obj@para$gene4use, "BCL2", "CDK8", "CDK19", "FLT3", "IDH1", "IDH2", "SHH", "TRDMT1"))

n_genes <- length(obj@para$gene4use)

obj <- Constr_net(obj, min_cells = 100, select_ratio = 1)

# Save the object
saveRDS(obj, file = paste0(project, "_", species, "_", target, "_", n_genes, "_", ncells, ".rds"))
