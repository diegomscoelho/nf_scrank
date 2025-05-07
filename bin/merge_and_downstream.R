#!/usr/bin/env Rscript
library(Seurat)
library(scRank)
library(dplyr)
library(GENIE3)


# rds_files <- c("/home/jamorim/data/scrank_results/Astrocytes_weight_GENIE3_50.rds", "/home/jamorim/data/scrank_results/Microglia_PVM_weight_GENIE3_500.rds")
# targets <- c("EGFR", "ACE")
# ad <- "/home/jamorim/scripts/nf_scrank/SEAD_res0.5_dementia_subtype.rds"
# column <- "broad_cell_type"
# species <- "human"
 
args <- commandArgs(trailingOnly = TRUE)

ad <- args[1]
targets <- args[2]
species <- args[3]
column <- args[4]
rds_files <- args[5:length(args)]


cell_types <- sub("_weight.*", "", basename(rds_files))
target <- targets[1]

sc_objs <- lapply(rds_files, readRDS)
seuratObj <- readRDS(ad)

obj <- CreateScRank(input = seuratObj,
                    species = species, 
                    cell_type = column,
                    target = target)

obj@net <- sc_objs
names(obj@net) <- cell_types

obj@para$ct.keep = names(obj@net)

saveRDS(obj, "merged_obj.RDS")

for (target_sc in targets){
  obj@para$target <- target_sc
  obj <- rank_celltype(obj)
  rank_df <- rbind(cell_types, obj@cell_type_rank$perb_score)
  write.table(rank_df, paste0("perbscore_", target_sc, ".txt"), quote = F, row.names = F, col.names = F)
}



