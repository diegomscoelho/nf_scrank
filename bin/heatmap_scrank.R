library(vroom)
library(dplyr)
require(ComplexHeatmap)
require(circlize)

perb_score <- read.table("../scripts/nf_scrank/results/rank_scores/perbscore_all_targets.txt", header = T)

perb_score <- perb_score %>%
  mutate(perb_score = ifelse(perb_score <= 0, NA, perb_score)) %>%
  mutate(perb_score = log10(perb_score))

# Optional: remove NAs introduced by log10
perb_score <- na.omit(perb_score)

# Rescale across the entire dataset (not grouped)
perb_score <- perb_score %>%
  mutate(perb_score = scales::rescale(perb_score, to = c(0, 1)))

heatmap_matrix <- reshape2::acast(perb_score, cell_type ~ target, value.var = "perb_score")

Heatmap(
  heatmap_matrix,
  name = "log10\nscaled\nperturbation\nscore",               # Name for the color legend
  row_title = "Cell Types",           # Title for rows
  column_title = "Alzheimer first trial",           # Title for columns
  cluster_rows = F,
  cluster_columns = F,
  show_heatmap_legend = T
)
