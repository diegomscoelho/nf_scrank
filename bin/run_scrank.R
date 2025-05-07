#!/usr/bin/env Rscript
library(Seurat)
library(scRank)
library(dplyr)
library(GENIE3)

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
n_cores <- args[2]
n_cores <- as.integer(n_cores)

cell_type <- sub(".RDS", "", seuratObj)

sc_obj <- readRDS(seuratObj)

start_time <- Sys.time()

weight <- GENIE3(as.matrix(sc_obj@data$seuratObj[sc_obj@para$gene4use]@assays$RNA$counts), nCores = n_cores)
weight <- weight[colnames(weight),]

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

n_cells <- dim(sc_obj@data$seuratObj)[2]

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_GENIE3_", n_cells, ".rds"))
