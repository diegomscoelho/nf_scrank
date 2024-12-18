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

print(seuratObj)

if (seuratObj == 'AML_object.rda') {
    load(seuratObj)
} else {
    seuratObj <- readRDS(seuratObj)
}


# Run scRank
obj <- CreateScRank(input = seuratObj,
                    species = species, 
                    cell_type = column,
                    target = target)
obj <- Constr_net(obj)
obj <- rank_celltype(obj)
obj <- init_mod(obj)

# Save the object
saveRDS(obj, file = paste0(species, "_", target, ".rds"))
