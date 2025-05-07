#!/usr/bin/env Rscript
library(Seurat)
library(scRank)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Inputs
seuratObj <- args[1]
targets <- args[2]
column <- args[3]
species <- args[4]
n_cells <- as.integer(args[5])


targets <- readLines(targets)
target <- targets[1]

if (seuratObj == 'AML_object.rda') {
    load(seuratObj)
    seuratObj <- seuratObj[c(VariableFeatures(seuratObj)[1:200], target),]
} else {
    seuratObj <- readRDS(seuratObj)
}

# Downsample cells by celltype
downsampled_cells <- seuratObj@meta.data %>% tibble::rowid_to_column("id_cell") %>%
  group_by(!!sym(column)) %>%
  slice_sample(n = n_cells) %>%
  pull(id_cell)

ncells <- length(downsampled_cells)
seurat_downsample <- seuratObj[, downsampled_cells]

split_obj <- SplitObject(seurat_downsample, split.by = column)

# Create scRank object
sc_obj <- lapply(split_obj, function(seuobj){
  seuobj <- FindVariableFeatures(seuobj, nfeatures = 2000)
  obj <- CreateScRank(input = seuobj,
                          species = species, 
                          cell_type = column,
                          target = target)
  obj@para$gene4use <- unique(c(obj@para$gene4use, targets))
  return(obj)
})


clean_name <- function(name) {
  gsub("[^A-Za-z0-9_\\-]", "_", name)  # Replace any non-safe character with "_"
}

# Save each object with a cleaned file name
invisible(lapply(names(sc_obj), function(name) {
  file_name <- paste0(clean_name(name), ".RDS")
  saveRDS(sc_obj[[name]], file = file_name)
}))

