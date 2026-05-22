
library(tidyverse)
library(UpSetR)
library(pheatmap)
library(RColorBrewer)
library(grid)

tissue_files <- list(
  "Breast"      = "Breast.csv",
  "Lung"        = "BronchusLung.csv",
  "Colon"       = "Colon.csv",
  "Liver"       = "Liver.csv",
  "Stomach"     = "Stomach.csv",
  "Prostate"    = "Prostate.csv",
  "Thyroid"     = "Thyroid.csv",
  "Kidney"      = "Kidney.csv"
)

upset_out   <- "upset_essential_reactions.png"
heatmap_out <- "pathway_heatmap.png"

tissue_data <- lapply(names(tissue_files), function(tissue) {
  path <- tissue_files[[tissue]]
  df   <- read.csv(path)
  tibble(tissue = tissue, Reaction = df$SingleLethal)
})

all_data <- bind_rows(tissue_data)

all_data <- all_data %>%
  left_join(map_df, by = "Reaction")


all_reactions <- unique(all_data$Reaction)
tissue_names  <- names(tissue_files)

binary_mat <- map_dfc(tissue_names, function(t) {
  as.integer(all_reactions %in% filter(all_data, tissue == t)$Reaction)
}) %>%
  set_names(tissue_names) %>%
  as.data.frame()

rownames(binary_mat) <- all_reactions


# =============================================================================
# 4. UPSET PLOT
# =============================================================================


tissue_colours <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(tissue_names)),
  tissue_names
)
png(upset_out, width = 3200, height = 1800, res = 300)

upset(
  binary_mat,
  sets              = tissue_names,
  sets.bar.color    = tissue_colours,
  order.by          = "freq",
  decreasing        = TRUE,
  mb.ratio          = c(0.6, 0.4),
  number.angles     = 0,
  point.size        = 3,
  line.size         = 1,
  text.scale        = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.1),
  mainbar.y.label   = "Number of Singlelethal Reactions",
  sets.x.label      = "Total Essential\nReactions per Tissue",
  keep.order        = FALSE,
  show.numbers      = "yes"
)
dev.off()


# =============================================================================
# PATHWAY HEATMAP
# =============================================================================

# Count essential reactions per subsystem per tissue
pathway_counts <- all_data %>%
  filter(!is.na(Subsystem)) %>%
  count(tissue, Subsystem, name = "n_essential") %>%
  pivot_wider(names_from = tissue, values_from = n_essential, values_fill = 0) %>%
  column_to_rownames("Subsystem") %>%
  as.matrix()

pathway_counts <- pathway_counts[-17, ] #remove unassigned pathway row
min_threshold <- 2
pathway_counts <- pathway_counts[apply(pathway_counts, 1, max) >= min_threshold, ]


heatmap_colours <- colorRampPalette(c("white", "#C6DBEF", "#2171B5", "#08306B"))(100)

png(heatmap_out, width = 3200, height = 2400, res = 300)

pheatmap(
  pathway_counts,
  color             = heatmap_colours,
  border_color      = "grey85",
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  clustering_method = "ward.D2",
  scale             = "none",          
  fontsize_row      = 9,
  fontsize_col      = 10,
  angle_col         = 45,
  main              = "Essential Reactions\nacross Cancer Tissue Models",
  cellwidth         = 15,
  cellheight        = 10
)

dev.off()

# Per-tissue counts
per_tissue <- all_data %>%
  group_by(tissue) %>%
  summarise(n_essential = n_distinct(Reaction), .groups = "drop") %>%
  arrange(desc(n_essential))

# Shared across ALL 8 tissues
core_essential <- all_reactions[rowSums(binary_mat) == length(tissue_names)]

# Unique to exactly one tissue
unique_essential <- map(tissue_names, function(t) {
  all_reactions[binary_mat[[t]] == 1 & rowSums(binary_mat) == 1]
}) %>% set_names(tissue_names)

cat("\nTissue-unique essential reactions:\n")
iwalk(unique_essential, ~cat(" ", .y, ":", length(.x), "\n"))


