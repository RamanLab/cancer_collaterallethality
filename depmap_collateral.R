library(depmap)
library(dplyr)
library(purrr)
library(readr)
library(tidyverse)

#retrieve datasets from depmap
chronos_data <- depmap_crispr()
expression_data <- depmap_TPM()
mutation_data <- depmap_mutationCalls()
copy_number_data <- depmap_copyNumber()
#earlier copy_number <=-0.5 meant for deletions now with log_copy_number 
#available doing this modification
# copy_number_data <- copy_number_data%>%
  # mutate(relative_copy_number = 2^log_copy_number - 1)
metadata <- depmap_metadata()

#our tissues and cancers of interest
tissues <- c("breast", "kidney", "liver", "colon", "lung", "thyroid", "stomach", "prostate")
cancer_types <- c(
  "Breast Cancer",
  "Colon/Colorectal Cancer",
  "Gastric Cancer",
  "Kidney Cancer",
  "Liver Cancer",
  "Lung Cancer",
  "Prostate Cancer",
  "Thyroid Cancer"
)

depmap_ids <- metadata %>%
  filter(
    primary_disease %in% cancer_types,
    tolower(sample_collection_site) %in% tissues
  ) %>%
  pull(depmap_id)


#filter all depmap datasets for only our tissues of interest
chronos_data <- chronos_data %>% filter(depmap_id %in% depmap_ids)
expression_data <- expression_data %>% filter(depmap_id %in% depmap_ids)
mutation_data <- mutation_data %>% filter(depmap_id %in% depmap_ids)
copy_number_data <- copy_number_data %>% filter(depmap_id %in% depmap_ids)
metadata <- metadata %>% filter(depmap_id %in% depmap_ids)

#check collateral lethal genes against depmap
evaluate_pair <- function(input_tissue, geneA, geneB) {
  depmap_ids_tissue <- metadata %>%
    filter(tolower(sample_collection_site) == tolower(input_tissue)) %>%
    pull(depmap_id)
  
  #Identify cell lines with geneB mutated
  geneB_inactivated <- mutation_data %>%
    filter(depmap_id %in% depmap_ids_tissue,
           gene_name == geneB,
           var_annotation != "silent") %>%
    pull(depmap_id) %>%
    c(copy_number_data %>%
        filter(depmap_id %in% depmap_ids_tissue,
               gene_name == geneB,
               log_copy_number <= 1) %>%
        pull(depmap_id)) %>%
    unique()
  
  if (length(geneB_inactivated) == 0) {
    return(tibble(
      Input_Tissue = input_tissue,
      Gene_A = geneA,
      Gene_B = geneB,
      N_inactivated_CLs = 0,
      N_CRISPR_essential = 0,
      N_essential_expressed = 0,
      Label = "notvalid"
    ))
  }
  
#Check if geneA is essential and expressed in geneB-inactivated lines
geneA_data <- chronos_data %>%
    filter(depmap_id %in% geneB_inactivated,
           gene_name == geneA) %>%
    left_join(expression_data %>%
                filter(depmap_id %in% geneB_inactivated,
                       gene_name == geneA),
              by = c("depmap_id", "gene_name")) %>%
    mutate(
      is_essential = dependency < -0.5,
      is_expressed = rna_expression > 4,  #rna_expression is log2(TPM+1) > 2.5 ~ TPM > 5.7
      meets_criteria = is_essential & is_expressed
    )
  
n_essential <- sum(geneA_data$is_essential, na.rm = TRUE)
n_essential_expressed <- sum(geneA_data$meets_criteria, na.rm = TRUE)
  
label <- ifelse(n_essential_expressed / length(geneB_inactivated) >= 0.3, "valid", "notvalid")
  
  tibble(
    Input_Tissue = input_tissue,
    Gene_A = geneA,
    Gene_B = geneB,
    N_inactivated_CLs = length(geneB_inactivated),
    N_CRISPR_essential = n_essential,
    N_essential_expressed = n_essential_expressed,
    Label = label
  )
}


#geneA is the collateral lethal/essential gene in cancer
#geneB is the gene which is inactivated or mutated in cancer
#but active in the normal cells

gene_group_input <- read.csv("SupplementaryTable1.csv",sep=',')
  
#do for all pairwise combinations
gene_pairs <- gene_group_input %>%
    mutate(
      geneA = str_split(geneA, ",\\s*"),
      geneB = str_split(geneB, ",\\s*")
    ) %>%
    unnest(geneA) %>%
    unnest(geneB) %>%
    filter(geneA != geneB)%>%
    rename(input_tissue = Tissue) %>%
    filter(tolower(input_tissue) %in% tissues)%>%
    distinct(input_tissue,geneA,geneB)
  
#Evaluate all gene pairs
results <- gene_pairs %>%
  pmap_dfr(evaluate_pair) %>%
  arrange(Input_Tissue, Gene_A, Gene_B)
write_csv(results, "SupplementaryTable2.csv")
