library(readr)
library(fgsea)
library(purrr)
library(dplyr)
library(tidyverse)
library(ggplot2)

# load in the matrisome gene set
matrisome_all <- read.csv("/Users/sm2949/Desktop/Dr_Matrisome_Masterlist_Nauroy et al_2017.xlsx - Dr_Matrisome_Masterlist.csv")

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

# create gene sets (gene subsets)
gene_sets <- list(
  Collagens = collagen_genes,
  ECM_Regulators = ecm_regulator_genes,
  Secreted_Factors = secreted_factors_genes,
  ECM_Affiliated_Proteins = affiliated_proteins_genes,
  Proteoglycans = proteoglycans_genes,
  ECM_Glycoproteins = ecm_glycoproteins_genes
)

# combining all terms together
gene_sets <- list(Matrisome = matrisome_all$Zebrafish.Gene.Symbol)

# --------------- helper functions !!!!!!! ---------------#
# function to prepare ranked list by filtering p-value and sorting by log2FC
prepare_ranked_list <- function(df) {
  rankings <- sign(df$avg_log2FC)*(-log10(df$p_val_adj))
  names(rankings) <- rownames(df)
  rankings <- sort(rankings, decreasing = TRUE)
  return(rankings)
}

# function to run GSEA for each timepoint
run_gsea <- function(ranked_list, gene_sets) {
  gsea_results <- fgsea(pathways = gene_sets,
                        stats = ranked_list+rnorm(length(ranked_list), sd=0.001),
                        scoreType = 'std',
                        minSize = 5,
                        eps = 0,
                        maxSize = 5000,
                        nproc = 1)
  return(gsea_results)
}

# extract NES for each timepoint
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
hep4dpa_ranked <- prepare_ranked_list(hep.4)
hep7dpa_ranked <- prepare_ranked_list(hep.7)

# run GSEA for each timepoint
gsea_hep0dpa <- run_gsea(hep0dpa_ranked, gene_sets)
gsea_hep1dpa <- run_gsea(hep1dpa_ranked, gene_sets)
gsea_hep2dpa <- run_gsea(hep2dpa_ranked, gene_sets)
gsea_hep3dpa <- run_gsea(hep3dpa_ranked, gene_sets)
gsea_hep4dpa <- run_gsea(hep4dpa_ranked, gene_sets)
gsea_hep7dpa <- run_gsea(hep7dpa_ranked, gene_sets)

# extract NES for each timepoint
nes_hep0dpa <- get_nes(gsea_hep0dpa)
nes_hep1dpa <- get_nes(gsea_hep1dpa)
nes_hep2dpa <- get_nes(gsea_hep2dpa)
nes_hep3dpa <- get_nes(gsea_hep3dpa)
nes_hep4dpa <- get_nes(gsea_hep4dpa)
nes_hep7dpa <- get_nes(gsea_hep7dpa)

# merge NES data from all timepoints
merged_nes_hep <- reduce(list(nes_hep0dpa, nes_hep1dpa, nes_hep2dpa, nes_hep3dpa, nes_hep4dpa, nes_hep7dpa), full_join, by = "Pathway")

# rename columns to indicate timepoints
colnames(merged_nes_hep) <- c("Pathway", "hep0dpa", "hep1dpa", "hep2dpa", "hep3dpa", "hep4dpa", "hep7dpa")

# filter for pathways that are present in all timepoints 
common_pathways_hep <- merged_nes_hep[complete.cases(merged_nes_hep[, -1]), ]

# plot NES over time for common pathways
long_nes_hep <- common_pathways_hep %>%
  pivot_longer(cols = starts_with("hep"), names_to = "Timepoint", values_to = "NES")

# plot
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
HSC4dpa_ranked <- prepare_ranked_list(HSC.4)
HSC7dpa_ranked <- prepare_ranked_list(HSC.7)

# run GSEA for each timepoint
gsea_HSC0dpa <- run_gsea(HSC0dpa_ranked, gene_sets)
gsea_HSC1dpa <- run_gsea(HSC1dpa_ranked, gene_sets)
gsea_HSC2dpa <- run_gsea(HSC2dpa_ranked, gene_sets)
gsea_HSC3dpa <- run_gsea(HSC3dpa_ranked, gene_sets)
gsea_HSC4dpa <- run_gsea(HSC4dpa_ranked, gene_sets)
gsea_HSC7dpa <- run_gsea(HSC7dpa_ranked, gene_sets)

# extract NES for each timepoint
nes_HSC0dpa <- get_nes(gsea_HSC0dpa)
nes_HSC1dpa <- get_nes(gsea_HSC1dpa)
nes_HSC2dpa <- get_nes(gsea_HSC2dpa)
nes_HSC3dpa <- get_nes(gsea_HSC3dpa)
nes_HSC4dpa <- get_nes(gsea_HSC4dpa)
nes_HSC7dpa <- get_nes(gsea_HSC7dpa)

# merge NES data from all timepoints
merged_nes_HSC <- reduce(list(nes_HSC0dpa, nes_HSC1dpa, nes_HSC2dpa, nes_HSC3dpa, nes_HSC4dpa, nes_HSC7dpa), full_join, by = "Pathway")

# rename columns to indicate timepoints
colnames(merged_nes_HSC) <- c("Pathway", "HSC0dpa", "HSC1dpa", "HSC2dpa", "HSC3dpa", "HSC4dpa", "HSC7dpa")

# filter for pathways that are present in all timepoints 
common_pathways_HSC <- merged_nes_HSC[complete.cases(merged_nes_HSC[, -1]), ]

# plot NES over time for common pathways
long_nes_HSC <- common_pathways_HSC %>%
  pivot_longer(cols = starts_with("HSC"), names_to = "Timepoint", values_to = "NES")

# plot
ggplot(long_nes_HSC, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line() +
  geom_point() +
  labs(title = "NES Over Regeneration In HSCs for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()

# ----------------- run on endothelial cells ------------------ #
end0dpa_ranked <- prepare_ranked_list(end.0)
end1dpa_ranked <- prepare_ranked_list(end.1)
end2dpa_ranked <- prepare_ranked_list(end.2)
end3dpa_ranked <- prepare_ranked_list(end.3)
end4dpa_ranked <- prepare_ranked_list(end.4)
end7dpa_ranked <- prepare_ranked_list(end.7)

# run GSEA for each timepoint
gsea_end0dpa <- run_gsea(end0dpa_ranked, gene_sets)
gsea_end1dpa <- run_gsea(end1dpa_ranked, gene_sets)
gsea_end2dpa <- run_gsea(end2dpa_ranked, gene_sets)
gsea_end3dpa <- run_gsea(end3dpa_ranked, gene_sets)
gsea_end4dpa <- run_gsea(end4dpa_ranked, gene_sets)
gsea_end7dpa <- run_gsea(end7dpa_ranked, gene_sets)

# extract NES for each timepoint
nes_end0dpa <- get_nes(gsea_end0dpa)
nes_end1dpa <- get_nes(gsea_end1dpa)
nes_end2dpa <- get_nes(gsea_end2dpa)
nes_end3dpa <- get_nes(gsea_end3dpa)
nes_end4dpa <- get_nes(gsea_end4dpa)
nes_end7dpa <- get_nes(gsea_end7dpa)

# merge NES data from all timepoints
merged_nes_end <- reduce(list(nes_end0dpa, nes_end1dpa, nes_end2dpa, nes_end3dpa, nes_end4dpa,nes_end7dpa), full_join, by = "Pathway")

# rename columns to indicate timepoints
colnames(merged_nes_end) <- c("Pathway", "end0dpa", "end1dpa", "end2dpa", "end3dpa", "end4dpa","end7dpa")

# filter for pathways that are present in all timepoints 
common_pathways_end <- merged_nes_end[complete.cases(merged_nes_end[, -1]), ]

# plot NES over time for common pathways
long_nes_end <- common_pathways_end %>%
  pivot_longer(cols = starts_with("end"), names_to = "Timepoint", values_to = "NES")

# plot
ggplot(long_nes_end, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line() +
  geom_point() +
  labs(title = "NES Over Regeneration In Endothelial Cells for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()

# ------------------------- creating a plot of ALL cell types-------------------------- #
long_nes_HSC$CellType <- "HSC"
long_nes_end$CellType <- "Endothelial"
long_nes_hep$CellType <- "Hepatocyte"

# combine into one df
long_nes_all <- bind_rows(long_nes_HSC, long_nes_end, long_nes_hep)

# remove the cell name from the timepoint column
long_nes_all$Timepoint <- gsub("^[a-zA-Z]+", "", long_nes_all$Timepoint)

ggplot(long_nes_all, aes(x = Timepoint, y = NES, color = CellType, group = CellType)) +
  geom_line(size = 1) +  
  geom_point(size = 1.2) +   
  scale_color_manual(values = c("Hepatocyte" = "#011a51", 
                                "Endothelial" = "#fab129", 
                                "HSC" = "#f06c00")) + 
  labs(title = "NES Over Regeneration for Matrisome Gene Sets Across Cell Types", 
       x = "Timepoint (dpa)", 
       y = "Normalized Enrichment Score (NES)", 
       color = "Cell Type") +
  theme_minimal(base_size = 14) +  # Increase base font size for readability
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), # Rotate x-axis labels
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.title = element_text(size = 12, face = "bold"), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Centered title
    legend.position = "right", 
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a box around the plot
  )

# -------------------- create the heatmap --------------------- #
# read in the file containing leading gene edge info for each cell type
leading_edges <- read_csv('/Users/sm2949/Desktop/2019SingleCell_leadingEdgeGenes.csv')

# extract the genes for each timepoint
leading_0dpa_genes <- leading_edges$'0dpa'
leading_1dpa_genes <- leading_edges$'1dpa'
leading_2dpa_genes <- leading_edges$'2dpa'
leading_3dpa_genes <- leading_edges$'3dpa'
leading_4dpa_genes <- leading_edges$'4dpa'
leading_7dpa_genes <- leading_edges$'7dpa'

