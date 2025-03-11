library(readr)
library(fgsea)
library(purrr)
library(dplyr)
library(tidyverse)
library(ggplot2)

# load in the matrisome gene set
matrisome_all <- read.csv("matrisome_list.csv")

# separate into different ECM categories
collagens <- matrisome_all[matrisome_all$Matrisome.Category == "Collagens",]
collagen_genes <- collagens$Zebrafish.Gene.Symbol

ecm_regulators <- matrisome_all[matrisome_all$Matrisome.Category == "ECM Regulators",]
ecm_regulator_genes <- ecm_regulators$Zebrafish.Gene.Symbol

secreted_factors <- matrisome_all[matrisome_all$Matrisome.Category == "Secreted Factors",]
secreted_factors_genes <- secreted_factors$Zebrafish.Gene.Symbol

affiliated_proteins <- matrisome_all[matrisome_all$Matrisome.Category == "ECM-affiliated Proteins",]
affiliated_proteins_genes <- affiliated_proteins$Zebrafish.Gene.Symbol

proteoglycans <- matrisome_all[matrisome_all$Matrisome.Category == "Proteoglycans",]
proteoglycans_genes <- proteoglycans$Zebrafish.Gene.Symbol

ecm_glycoproteins <- matrisome_all[matrisome_all$Matrisome.Category == "ECM Glycoproteins",]
ecm_glycoproteins_genes <- ecm_glycoproteins$Zebrafish.Gene.Symbol

# Create gene sets (gene subsets)
gene_sets <- list(
  Collagens = collagen_genes,
  ECM_Regulators = ecm_regulator_genes,
  Secreted_Factors = secreted_factors_genes,
  ECM_Affiliated_Proteins = affiliated_proteins_genes,
  Proteoglycans = proteoglycans_genes,
  ECM_Glycoproteins = ecm_glycoproteins_genes
)

# --------------- helper functions !!!!!!! ---------------#
# Function to prepare ranked list by filtering p-value and sorting by log2FC
prepare_ranked_list <- function(df) {
  rankings <- sign(df$avg_log2FC)*(-log10(df$p_val_adj))
  names(rankings) <- rownames(df)
  rankings <- sort(rankings, decreasing = TRUE)
  return(rankings)
}

# Function to run GSEA for each timepoint
run_gsea <- function(ranked_list, gene_sets) {
  gsea_results <- fgsea(pathways = gene_sets,
                        stats = ranked_list+rnorm(length(ranked_list), sd=0.001),
                        scoreType = 'std',
                        minSize = 10,
                        eps = 0,
                        maxSize = 500,
                        nproc = 1)
  return(gsea_results)
}

# Extract NES for each timepoint
get_nes <- function(gsea_results) {
  nes <- gsea_results$NES 
  nes_df <- data.frame(Pathway = gsea_results$pathway, NES = nes)
  return(nes_df)
}

# ----------------- run on hepatocytes ------------------ #
hep0dpa_ranked <- prepare_ranked_list(hep.0)
hep1dpa_ranked <- prepare_ranked_list(hep.1)
hep2dpa_ranked <- prepare_ranked_list(hep.2)
hep3dpa_ranked <- prepare_ranked_list(hep.3)
hep7dpa_ranked <- prepare_ranked_list(hep.7)

# Run GSEA for each timepoint
gsea_hep0dpa <- run_gsea(hep0dpa_ranked, gene_sets)
gsea_hep1dpa <- run_gsea(hep1dpa_ranked, gene_sets)
gsea_hep2dpa <- run_gsea(hep2dpa_ranked, gene_sets)
gsea_hep3dpa <- run_gsea(hep3dpa_ranked, gene_sets)
gsea_hep7dpa <- run_gsea(hep7dpa_ranked, gene_sets)

# Extract NES for each timepoint
nes_hep0dpa <- get_nes(gsea_hep0dpa)
nes_hep1dpa <- get_nes(gsea_hep1dpa)
nes_hep2dpa <- get_nes(gsea_hep2dpa)
nes_hep3dpa <- get_nes(gsea_hep3dpa)
nes_hep7dpa <- get_nes(gsea_hep7dpa)

# Merge NES data from all timepoints
merged_nes_hep <- reduce(list(nes_hep0dpa, nes_hep1dpa, nes_hep2dpa, nes_hep3dpa, nes_hep7dpa), full_join, by = "Pathway")

# Rename columns to indicate timepoints
colnames(merged_nes_hep) <- c("Pathway", "hep0dpa", "hep1dpa", "hep2dpa", "hep3dpa", "hep7dpa")

# Filter for pathways that are present in all timepoints 
common_pathways_hep <- merged_nes_hep[complete.cases(merged_nes_hep[, -1]), ]

# Plot NES over time for common pathways
long_nes_hep <- common_pathways_hep %>%
  pivot_longer(cols = starts_with("hep"), names_to = "Timepoint", values_to = "NES")

ggplot(long_nes_hep, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line() +
  geom_point() +
  labs(title = "NES Over Regeneration In Hepatocytes for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()

# ----------------- run on HSCs ------------------ #
HSC0dpa_ranked <- prepare_ranked_list(HSC.0)
HSC1dpa_ranked <- prepare_ranked_list(HSC.1)
HSC2dpa_ranked <- prepare_ranked_list(HSC.2)
HSC3dpa_ranked <- prepare_ranked_list(HSC.3)
HSC7dpa_ranked <- prepare_ranked_list(HSC.7)

# Run GSEA for each timepoint
gsea_HSC0dpa <- run_gsea(HSC0dpa_ranked, gene_sets)
gsea_HSC1dpa <- run_gsea(HSC1dpa_ranked, gene_sets)
gsea_HSC2dpa <- run_gsea(HSC2dpa_ranked, gene_sets)
gsea_HSC3dpa <- run_gsea(HSC3dpa_ranked, gene_sets)
gsea_HSC7dpa <- run_gsea(HSC7dpa_ranked, gene_sets)

# Extract NES for each timepoint
nes_HSC0dpa <- get_nes(gsea_HSC0dpa)
nes_HSC1dpa <- get_nes(gsea_HSC1dpa)
nes_HSC2dpa <- get_nes(gsea_HSC2dpa)
nes_HSC3dpa <- get_nes(gsea_HSC3dpa)
nes_HSC7dpa <- get_nes(gsea_HSC7dpa)

# Merge NES data from all timepoints
merged_nes_HSC <- reduce(list(nes_HSC0dpa, nes_HSC1dpa, nes_HSC2dpa, nes_HSC3dpa, nes_HSC7dpa), full_join, by = "Pathway")

# Rename columns to indicate timepoints
colnames(merged_nes_HSC) <- c("Pathway", "HSC0dpa", "HSC1dpa", "HSC2dpa", "HSC3dpa", "HSC7dpa")

# Filter for pathways that are present in all timepoints 
common_pathways_HSC <- merged_nes_HSC[complete.cases(merged_nes_HSC[, -1]), ]

# Plot NES over time for common pathways
long_nes_HSC <- common_pathways_HSC %>%
  pivot_longer(cols = starts_with("HSC"), names_to = "Timepoint", values_to = "NES")

ggplot(long_nes_HSC, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line() +
  geom_point() +
  labs(title = "NES Over Regeneration In HSCs for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()