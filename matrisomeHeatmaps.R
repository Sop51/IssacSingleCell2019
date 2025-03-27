library(readr)
library(pheatmap)

# load in the matrisome gene set
matrisome_all <- read.csv('/Users/sm2949/Desktop/Dr_Matrisome_Masterlist_Nauroy et al_2017.xlsx - Dr_Matrisome_Masterlist.csv')

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

all_matrisome_genes <- list(Matrisome = matrisome_all$Zebrafish.Gene.Symbol)

# helper functions !!!!!!!
filter_matrisome_genes <- function(df) {
  matrisome_genes <- all_matrisome_genes$Matrisome
  df_filtered <- df[df$avg_log2FC > 1 & df$p_val_adj < 0.05, ]  # Filter for significance
  df_filtered <- df_filtered[rownames(df_filtered) %in% matrisome_genes , ]  # Keep only matrisome genes
  return(df_filtered)
}

generate_heatmap_matrix <- function(timepoint_dfs) {
  # Extract all unique gene names across timepoints
  all_genes <- unique(unlist(lapply(timepoint_dfs, rownames)))
  
  # Create a matrix to store log2FC values
  heatmap_matrix <- data.frame(Gene = all_genes)
  
  # Fill in log2FC values for each timepoint
  for (timepoint in names(timepoint_dfs)) {
    df <- timepoint_dfs[[timepoint]]
    
    heatmap_matrix[[timepoint]] <- df$avg_log2FC[match(heatmap_matrix$Gene, rownames(df))]
  }
  
  # Replace NAs with 0
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Convert to matrix format
  rownames(heatmap_matrix) <- heatmap_matrix$Gene
  heatmap_matrix <- as.matrix(heatmap_matrix[, -1])
  
  return(heatmap_matrix)
}

# ----------------- run on hepatocytes ------------------ #
hep0dpa_ranked <- filter_matrisome_genes(hep.0)
hep1dpa_ranked <- filter_matrisome_genes(hep.1)
hep2dpa_ranked <- filter_matrisome_genes(hep.2)
hep3dpa_ranked <- filter_matrisome_genes(hep.3)
hep4dpa_ranked <- filter_matrisome_genes(hep.4)
hep7dpa_ranked <- filter_matrisome_genes(hep.7)

# list of ranked timepoint dataframes
timepoint_hep_dfs <- list(
  "0dpa" = hep0dpa_ranked,
  "1dpa" = hep1dpa_ranked,
  "2dpa" = hep2dpa_ranked,
  "3dpa" = hep3dpa_ranked,
  "4dpa" = hep4dpa_ranked,
  "7dpa" = hep7dpa_ranked
)

hep_heatmap_matrix <- generate_heatmap_matrix(timepoint_hep_dfs)

# generate heatmap
pheatmap(
  hep_heatmap_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),  
  cluster_rows = TRUE, 
  cluster_cols = FALSE, 
  scale = "row",  
  fontsize_row = 7,  
  fontsize_col = 10,
  main = "Differentially Expressed Matrisome Genes in Hepatocytes Across Regeneration"
)
