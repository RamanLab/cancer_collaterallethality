library(readxl)
library(dplyr)
library(tidyr)

## =========================
## 1. READ FILES
## =========================

read_excel("recon_genes_reactions_subsytems.xlsx")
recon_genes_reactions_subsytems = recon_genes_reactions_subsytems[,1:3]
map_df = recon_genes_reactions_subsytems
## =========================
## 2. EXPAND GENE COLUMN
## =========================

map_expanded <- map %>%
  filter(!is.na(Genes)) %>%
  separate_rows(Genes, sep = ",\\s*") %>%   # split "A, B, C"
  mutate(Genes = trimws(Genes)) %>%
  distinct()


bg_genes <- unique(map_expanded$Genes)
pathways <- unique(map_expanded$Subsystem)

# make a gmt file for GSEA from the recon model gene pathway mapping 
gmt_df <- map_expanded %>%
  group_by(Subsystem) %>%
  summarise(
    genes = paste(unique(Genes), collapse = "\t"),
    .groups = "drop"
  )

gmt_lines <- paste(
  gmt_df$Subsystem,
  "NA",
  gmt_df$genes,
  sep = "\t"
)

writeLines(gmt_lines, "subsystem_gene_sets.gmt")

gmt_rxn <- map_expanded %>%
  filter(!is.na(Reaction), Reaction != "") %>%   # clean
  group_by(Subsystem) %>%
  summarise(
    reactions = paste(unique(Reaction), collapse = "\t"),
    .groups = "drop"
  )

gmt_rxn_lines <- paste(
  gmt_rxn$Subsystem,
  "NA",
  gmt_rxn$reactions,
  sep = "\t"
)

writeLines(gmt_rxn_lines, "subsystem_reaction_sets.gmt")

#DEF of reactions after flux sampling 
library(limma)
library(tibble)

files <- list.files(
  "/home/maziya/cancer_collateral/flux_sampling/flux_sampling_DEF_limma",
  pattern = "\\.csv$",
  full.names = TRUE
)

camera_results_all <- list()

for(file in files){
  
  tissue <- tools::file_path_sans_ext(basename(file))
  
  object <- read.csv(file, row.names = 1, check.names = FALSE)
  object <- object[rowSums(object) > 0, ]
  object <- log2(object + 1)
  object <- object[apply(object, 1, function(x) all(is.finite(x))), ]
  groups <- ifelse(grepl("[A-Z]N\\d+$", colnames(object)), "Healthy",
                   ifelse(grepl("[A-Z]T\\d+$", colnames(object)), "Tumor", NA))
  
  groups <- factor(groups, levels = c("Healthy", "Tumor"))
  design <- model.matrix(~0 + groups)
  colnames(design) <- c("Healthy", "Tumor")
  cont_matrix <- makeContrasts(TumorvsHealthy = Tumor - Healthy, levels = design)
  
  # --- limma differential flux ---
  fit <- lmFit(object, design)
  fit <- contrasts.fit(fit, cont_matrix)
  fit <- eBayes(fit)
  t_stats <- fit$t[, "TumorvsHealthy"] 
  top_rxns <- topTable(fit, number = Inf, adjust = "BH")
  top_rxns <- rownames_to_column(top_rxns, "Reaction")
  
  write.csv(top_rxns,paste0(tissue, "_differential_flux.csv"),row.names = FALSE)
  
  # --- camera pathway enrichment ---
  index <- ids2indices(rxn_sets, id = rownames(object))
  index <- index[lengths(index) >= 10]
  
  cam_result <- cameraPR(
    statistic = t_stats,
    index     = index
  )%>%
    rownames_to_column("pathway") %>%
    mutate(tissue = tissue)
  
  camera_results_all[[tissue]] = cam_result
  
  write.csv(cam_result,paste0(tissue, "_camera_enrichment.csv"),row.names = FALSE)
  cat("Completed:", tissue, "\n")
}

# combine all tissues into one summary table
camera_summary <- bind_rows(camera_results_all)
write.csv(camera_summary,"all_tissues_camera_enrichment.csv",row.names = FALSE)
library(patchwork)

plot_list <- list()

file_diff <- list.files(
  "/home/maziya/cancer_collateral/flux_sampling/flux_sampling_DEF_limma/differential_flux",
  pattern = "\\.csv$",
  full.names = TRUE
)
legend_plot <- ggplot() +
  geom_point(aes(x = 1, y = 1, color = "Downregulated"), size = 3) +
  geom_point(aes(x = 1, y = 1, color = "Upregulated"), size = 3) +
  scale_color_manual(
    values = c("Upregulated" = "red",
               "Downregulated" = "blue"),
    name = ""
  ) +
  theme_void() +
  theme(
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.background = element_rect(color = "black"),
    legend.background = element_rect(fill = "white")
  )
for (file in file_diff) {
  tissue <- gsub("_differential_flux\\.csv", "", basename(file))
  res <- read.csv(file,check.names = FALSE)
  
  res <- res %>%
    mutate(
      logP = -log10(adj.P.Val),
      status = case_when(
        adj.P.Val < 0.05 & logFC > 0.5  ~ "Up",
        adj.P.Val < 0.05 & logFC < -0.5 ~ "Down",
        TRUE ~ "NS"
      )
    )
  top_up <- res %>%
    filter(status == "Up") %>%
    arrange(desc(logFC)) %>%  
    slice_head(n = 5)
  
  top_down <- res %>%
    filter(status == "Down") %>%
    arrange(logFC) %>%         
    slice_head(n = 5)
  top_up <- top_up %>%
    mutate(logFC_r = round(logFC, 3),
           logP_r  = round(logP, 3)) %>%
    distinct(logFC_r, logP_r, .keep_all = TRUE)
  
  top_down <- top_down %>%
    mutate(logFC_r = round(logFC, 3),
           logP_r  = round(logP, 3)) %>%
    distinct(logFC_r, logP_r, .keep_all = TRUE)
  
  top_labels <- bind_rows(top_up, top_down)
  
  p <- ggplot(res, aes(logFC, logP)) +
    
    geom_point(color = "grey70", alpha = 0.6, size = 1.5) +
    
    geom_point(data = subset(res, status == "Up"), color = "red", size = 2) +
    geom_point(data = subset(res, status == "Down"), color = "blue", size = 2) +
    
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    geom_text_repel(
      data = top_up,
      aes(label = Reaction),
      nudge_x = 5,                  
      direction = "y",
      hjust = 0,
      segment.curvature = -0.5,
      segment.angle = 20,
      segment.ncp = 3,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "black",
      fill = "white",
      size = 3,
      max.overlaps = Inf
    ) +
    
    geom_text_repel(
      data = top_down,
      aes(label = Reaction),
      nudge_x = -5,                
      direction = "y",
      hjust = 1,
      segment.curvature = 0.5,
      segment.angle = 20,
      segment.ncp = 3,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "black",
      fill = "white",
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      title = tissue,  
      x = expression(log[2] ~ "Fold Change"),
      y = expression(-log[10] * (FDR))
    ) +
    theme(
      text = element_text(family = "mono"),
      plot.title = element_text(hjust = 0.5,face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white", colour = NA))
  p <- p + theme(legend.position = "none")
  
  plot_list[[tissue]] <- p
  
}


plots_grid <- do.call(grid.arrange, c(plot_list, ncol = 4))

panel_plot <- cowplot::plot_grid(
  plots_grid,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.06),
  align = "v"
) +
  theme(
    plot.background = element_rect(fill = "white", colour = "white")
  )


ggsave("all_tissues_volcano_panel_Fig4.png",
       panel_plot, bg = "white",
       width = 16, height = 8, dpi = 600)


#load the all_tissues_camera enrichment csv file 
all_tissues_camera_enrichment = all_tissues_camera_enrichment %>%
  mutate(
    log10FDR  = -log10(FDR),
    Enrichment = factor(Direction, levels = c("Up", "Down"))
  ) %>%
  filter(FDR < 0.05) 


plot_dot = all_tissues_camera_enrichment %>%  ggplot(aes(x = tissue, y = pathway, size = log10FDR, colour = Enrichment)) +
  geom_point(alpha = 0.85) +
  scale_colour_manual(values = c("Up" = "#D85A30", "Down" = "#1D9E75")) +
  scale_size_continuous(
    name   = "-log10(FDR)",
    range  = c(3, 12),
    breaks = c(5, 10, 15, 20, 25)
  ) +
  theme_light(base_size = 14) +
  theme(
    text = element_text(family = "mono"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    x      = "Tissue",
    y      = "Pathway",
    title  = "Enriched pathways in tumors"
  )

ggsave("camera_pathwayenrichment_Fig5.png", plot = plot_dot,width = 10, height = 10)