library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyverse)
library(SeuratDisk)

zf <- LoadH5Seurat("/Users/sm2949/Desktop/SeuratProject.h5Seurat")

# make a copy of the seurat object
zf_plotting <- zf

# count number of cells per cell type
cell_counts <- zf@meta.data %>%
  count(cell.type.9.long) %>%
  arrange(desc(n))

# bar plot
ggplot(cell_counts, aes(x = reorder(cell.type.9.long, n), y = n, fill = cell.type.9.long)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 4.5) +  # add cell counts above bars
  scale_fill_viridis_d(option = "D", direction = -1) +  # colorblind-friendly discrete palette
  labs(
    title = "Number of Cells per Cell Type",
    x = "Cell Type",
    y = "Cell Count",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"  # hide legend
  ) +
  ylim(0, max(cell_counts$n) * 1.1)  # add space above bars for labels

# create a violin plot of anxa2a expression -----
cells <- subset(zf, idents = "Endothelial Cell")
# Extract data from Seurat object
data <- FetchData(cells, vars = c("anxa2a", "timepoint"))

# convert timepoint to factor for correct ordering
data$timepoint <- factor(data$timepoint, levels = c("mock", "0dpa", "1dpa", "2dpa", "3dpa", "4dpa", "7dpa"))

# Plot violin plot
ggplot(data, aes(x = timepoint, y = anxa2a, fill = timepoint)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Add boxplot for clarity
  scale_fill_manual(values = c("mock" = "gray", "0dpa" = "#ccd5ae", "1dpa" = "#e9c46a", "2dpa" = "#2a9d8f", "3dpa" = "#e63946", "7dpa" = "#eab69f")) +
  labs(title = "anxa2a Expression Across Timepoints in Endothelial Cells", x = "Timepoint", y = "Expression") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 

