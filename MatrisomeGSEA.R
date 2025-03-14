library(readr)
library(fgsea)
library(purrr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(edgeR)
library(ComplexHeatmap)

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

# NOTE: THIS SECTION IS FOR HEPATOCYTES 
# extract the genes for each timepoint
leading_0dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "hepatocytes", '0dpa'])
leading_1dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "hepatocytes", '1dpa'])
leading_2dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "hepatocytes", '2dpa'])
leading_3dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "hepatocytes", '3dpa'])
leading_4dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "hepatocytes", '4dpa'])
leading_7dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "hepatocytes", '7dpa'])

# combine into one list 
all_leading_genes_hep <- c(leading_0dpa_genes, 
                    leading_1dpa_genes, 
                    leading_2dpa_genes, 
                    leading_3dpa_genes, 
                    leading_4dpa_genes, 
                    leading_7dpa_genes)

# create a copy of the seurat object
zf_heat <- zf

# aggregate expression per timepoint per cell type
cts <- AggregateExpression(zf, 
                           group.by = c("cell.type.9.short", "timepoint"),
                           assays = 'RNA',
                           slot = 'counts',
                           return.seurat = FALSE)

# convert to a data frame
cts <- as.data.frame(cts$RNA)

# transpose the data frame
cts.t <- t(cts)

# get values where to split
splitrows <- gsub('_.*', '', rownames(cts.t))

# split the data frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitrows))

# fix the columns and transpose
cts.split.modified <- lapply(cts.split, function(x){
  # for each element remove the cell type name
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})

# normalize the data and generate for hepatocytes
# get the counts matrix
counts_hep <- cts.split.modified$Hep

# generate sample level metadata
colData <- data.frame(timepoint = colnames(counts_hep))

# set the condition to a factor vairable
groups <- as.factor(colData$timepoint)

# create the DGElist object
y <- DGEList(counts=counts_hep,group=groups)
# normalize by library size
y <- normLibSizes(y)
# calculate the normalization factors
y <- calcNormFactors(y)
normalized_hep <- cpm(y)
normalized_hep <- as.data.frame(normalized_hep)

# subset the count matrix for only the wanted genes
leading_edge_hep <- normalized_hep[rownames(normalized_hep) %in% all_leading_genes_hep, ]
leading_edge_hep <- as.data.frame(leading_edge_hep)
leading_edge_hep <- leading_edge_hep[, !colnames(leading_edge_hep) %in% c("mock")]
leading_edge_hep <- as.matrix(leading_edge_hep)

# Create a metadata frame to hold timepoints and KEGG categories for the genes
gene_metadata <- leading_edges %>% 
  filter(`cell type` == "hepatocytes") %>%
  select('0dpa', '1dpa', '2dpa', '3dpa', '4dpa', '7dpa', 'category') %>%
  gather(key = "timepoint", value = "gene", -category) %>%
  filter(gene %in% rownames(leading_edge_hep)) %>%
  distinct(gene, category, .keep_all = TRUE)

# Create a new data frame that matches the order of genes in `leading_edge_hep`
ordered_gene_metadata <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_hep)) %>%
  arrange(match(gene, rownames(leading_edge_hep))) %>%
  select(gene, category)

# Set the gene column as the row names
ordered_gene_metadata <- as.data.frame(ordered_gene_metadata)
rownames(ordered_gene_metadata) <- ordered_gene_metadata$gene
ordered_gene_metadata <- ordered_gene_metadata %>% select(-gene)

annotation_colors <- list(
  category = c("collagens" = "#f4f1de", 
               "ecm affiliated proteins" = "#eab69f",
               "ecm glycoproteins" = "#e07a5f", 
               "ecm regulators" = "#3d405b",
               "proteoglycans" = "#81b29a", 
               "secreted factors" = "#f2cc8f")
)

ha <- HeatmapAnnotation(category = ordered_gene_metadata$category,
                        which = 'row',
                        col = annotation_colors)

timepoint_split <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_hep)) %>%
  arrange(match(gene, rownames(leading_edge_hep))) %>%
  pull(timepoint)

# Z-score normalize the matrix by row (genes)
leading_edge_hep_zscored <- t(scale(t(leading_edge_hep)))

# Generate heatmap with annotations
Heatmap(leading_edge_hep_zscored, right_annotation = ha, 
        cluster_columns = FALSE, cluster_rows = TRUE,
        rect_gp = gpar(col = "white", lwd = 2),
        column_title = "Leading Edge Matrisome Genes Per Pathway In Hepatocytes",
        heatmap_legend_param = list(
          title = gt_render("<span style='color:black'>**Expression**</span>"), 
          at = c(-2, 0, 2), 
          labels = gt_render(c("Low Expression", "No Change", "High Expression"))
        ))


# NOTE: THIS SECTION IS FOR ENDOTHELIAL CELLS
# extract the genes for each timepoint
leading_0dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "endothelial", '0dpa'])
leading_1dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "endothelial", '1dpa'])
leading_2dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "endothelial", '2dpa'])
leading_3dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "endothelial", '3dpa'])
leading_4dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "endothelial", '4dpa'])
leading_7dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "endothelial", '7dpa'])

# combine into one list 
all_leading_genes_end <- c(leading_0dpa_genes, 
                           leading_1dpa_genes, 
                           leading_2dpa_genes, 
                           leading_3dpa_genes, 
                           leading_4dpa_genes, 
                           leading_7dpa_genes)

# normalize the data and generate for endothelial cells
# get the counts matrix
counts_end <- cts.split.modified$End

# generate sample level metadata
colData <- data.frame(timepoint = colnames(counts_end))

# set the condition to a factor vairable
groups <- as.factor(colData$timepoint)

# create the DGElist object
y <- DGEList(counts=counts_end,group=groups)
# normalize by library size
y <- normLibSizes(y)
# calculate the normalization factors
y <- calcNormFactors(y)
normalized_end <- cpm(y)
normalized_end <- as.data.frame(normalized_end)

# subset the count matrix for only the wanted genes
leading_edge_end <- normalized_end[rownames(normalized_end) %in% all_leading_genes_end, ]
leading_edge_end <- as.data.frame(leading_edge_end)
leading_edge_end <- leading_edge_end[, !colnames(leading_edge_end) %in% c("mock")]
leading_edge_end <- as.matrix(leading_edge_end)
# remove genes where all expression values are 0 across all columns
leading_edge_end <- leading_edge_end[rowSums(leading_edge_end != 0) > 0, ]

# Create a metadata frame to hold timepoints and KEGG categories for the genes
gene_metadata <- leading_edges %>% 
  filter(`cell type` == "endothelial") %>%
  select('0dpa', '1dpa', '2dpa', '3dpa', '4dpa', '7dpa', 'category') %>%
  gather(key = "timepoint", value = "gene", -category) %>%
  filter(gene %in% rownames(leading_edge_end)) %>%
  distinct(gene, category, .keep_all = TRUE)

# Create a new data frame that matches the order of genes
ordered_gene_metadata <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_end)) %>%
  arrange(match(gene, rownames(leading_edge_end))) %>%
  select(gene, category)

# Set the gene column as the row names
ordered_gene_metadata <- as.data.frame(ordered_gene_metadata)
rownames(ordered_gene_metadata) <- ordered_gene_metadata$gene
ordered_gene_metadata <- ordered_gene_metadata %>% select(-gene)

annotation_colors <- list(
  category = c("collagens" = "#f4f1de", 
               "ecm affiliated proteins" = "#eab69f",
               "ecm glycoproteins" = "#e07a5f", 
               "ecm regulators" = "#3d405b",
               "proteoglycans" = "#81b29a", 
               "secreted factors" = "#f2cc8f")
)

ha <- HeatmapAnnotation(category = ordered_gene_metadata$category,
                        which = 'row',
                        col = annotation_colors)

timepoint_split <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_end)) %>%
  arrange(match(gene, rownames(leading_edge_end))) %>%
  pull(timepoint)

# Z-score normalize the matrix by row (genes)
leading_edge_end_zscored <- t(scale(t(leading_edge_end)))

# Generate heatmap with annotations
Heatmap(leading_edge_end_zscored, right_annotation = ha, 
        cluster_columns = FALSE, cluster_rows = TRUE,
        rect_gp = gpar(col = "white", lwd = 2),
        column_title = "Leading Edge Matrisome Genes Per Pathway In Endothelial Cells",
        heatmap_legend_param = list(
          title = gt_render("<span style='color:black'>**Expression**</span>"), 
          at = c(-2, 0, 2), 
          labels = gt_render(c("Low Expression", "No Change", "High Expression"))
        ))

# NOTE: THIS SECTION IS FOR HSCs
# extract the genes for each timepoint
leading_0dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "HSC", '0dpa'])
leading_1dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "HSC", '1dpa'])
leading_2dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "HSC", '2dpa'])
leading_3dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "HSC", '3dpa'])
leading_4dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "HSC", '4dpa'])
leading_7dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "HSC", '7dpa'])

# combine into one list 
all_leading_genes_hsc <- c(leading_0dpa_genes, 
                           leading_1dpa_genes, 
                           leading_2dpa_genes, 
                           leading_3dpa_genes, 
                           leading_4dpa_genes, 
                           leading_7dpa_genes)

# normalize the data and generate for endothelial cells
# get the counts matrix
counts_hsc <- cts.split.modified$HSC

# generate sample level metadata
colData <- data.frame(timepoint = colnames(counts_hsc))

# set the condition to a factor vairable
groups <- as.factor(colData$timepoint)

# create the DGElist object
y <- DGEList(counts=counts_hsc,group=groups)
# normalize by library size
y <- normLibSizes(y)
# calculate the normalization factors
y <- calcNormFactors(y)
normalized_hsc <- cpm(y)
normalized_hsc <- as.data.frame(normalized_hsc)

# subset the count matrix for only the wanted genes
leading_edge_hsc <- normalized_hsc[rownames(normalized_hsc) %in% all_leading_genes_hsc, ]
leading_edge_hsc <- as.data.frame(leading_edge_hsc)
leading_edge_hsc <- leading_edge_hsc[, !colnames(leading_edge_hsc) %in% c("mock")]
leading_edge_hsc <- as.matrix(leading_edge_hsc)
# remove genes where all expression values are 0 across all columns
leading_edge_hsc <- leading_edge_hsc[rowSums(leading_edge_hsc != 0) > 0, ]

# Create a metadata frame to hold timepoints and KEGG categories for the genes
gene_metadata <- leading_edges %>% 
  filter(`cell type` == "HSC") %>%
  select('0dpa', '1dpa', '2dpa', '3dpa', '4dpa', '7dpa', 'category') %>%
  gather(key = "timepoint", value = "gene", -category) %>%
  filter(gene %in% rownames(leading_edge_hsc)) %>%
  distinct(gene, category, .keep_all = TRUE)

# Create a new data frame that matches the order of genes
ordered_gene_metadata <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_hsc)) %>%
  arrange(match(gene, rownames(leading_edge_hsc))) %>%
  select(gene, category)

# Set the gene column as the row names
ordered_gene_metadata <- as.data.frame(ordered_gene_metadata)
rownames(ordered_gene_metadata) <- ordered_gene_metadata$gene
ordered_gene_metadata <- ordered_gene_metadata %>% select(-gene)

annotation_colors <- list(
  category = c("collagens" = "#f4f1de", 
               "ecm affiliated proteins" = "#eab69f",
               "ecm glycoproteins" = "#e07a5f", 
               "ecm regulators" = "#3d405b",
               "proteoglycans" = "#81b29a", 
               "secreted factors" = "#f2cc8f")
)

ha <- HeatmapAnnotation(category = ordered_gene_metadata$category,
                        which = 'row',
                        col = annotation_colors)

timepoint_split <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_hsc)) %>%
  arrange(match(gene, rownames(leading_edge_hsc))) %>%
  pull(timepoint)

# Z-score normalize the matrix by row (genes)
leading_edge_hsc_zscored <- t(scale(t(leading_edge_hsc)))

# Generate heatmap with annotations
Heatmap(leading_edge_hsc_zscored, right_annotation = ha, 
        cluster_columns = FALSE, cluster_rows = TRUE,
        rect_gp = gpar(col = "white", lwd = 2),
        column_title = "Leading Edge Matrisome Genes Per Pathway In HSCs",
        heatmap_legend_param = list(
          title = gt_render("<span style='color:black'>**Expression**</span>"), 
          at = c(-2, 0, 2), 
          labels = gt_render(c("Low Expression", "No Change", "High Expression"))
        ))

