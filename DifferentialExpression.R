library(Seurat)
library(SeuratDisk)

zf <- LoadH5Seurat("/Users/sm2949/Desktop/SeuratProject.h5Seurat")

# make a copy of the seurat object
zf_de <- zf

# create a new column holding cell type and timepoint information
zf_de$celltype.status <- paste(zf_de$cell.type.9.short, zf_de$timepoint, sep = "_")

# set the identity to this new column
Idents(zf_de) <- "celltype.status"

# -------------- DE on hepatocytes --------------- #
hep.0 <- FindMarkers(zf_de, ident.1 = "Hep_0dpa", ident.2 = "Hep_mock", verbose = FALSE)
hep.1 <- FindMarkers(zf_de, ident.1 = "Hep_1dpa", ident.2 = "Hep_mock", verbose = FALSE)
hep.2 <- FindMarkers(zf_de, ident.1 = "Hep_2dpa", ident.2 = "Hep_mock", verbose = FALSE)
hep.3 <- FindMarkers(zf_de, ident.1 = "Hep_3dpa", ident.2 = "Hep_mock", verbose = FALSE)
hep.4 <- FindMarkers(zf_de, ident.1 = "Hep_4dpa", ident.2 = "Hep_mock", verbose = FALSE)
hep.7 <- FindMarkers(zf_de, ident.1 = "Hep_7dpa", ident.2 = "Hep_mock", verbose = FALSE)
# -------------- DE on hsc --------------- #
HSC.0 <- FindMarkers(zf_de, ident.1 = "HSC_0dpa", ident.2 = "HSC_mock", verbose = FALSE)
HSC.1 <- FindMarkers(zf_de, ident.1 = "HSC_1dpa", ident.2 = "HSC_mock", verbose = FALSE)
HSC.2 <- FindMarkers(zf_de, ident.1 = "HSC_2dpa", ident.2 = "HSC_mock", verbose = FALSE)
HSC.3 <- FindMarkers(zf_de, ident.1 = "HSC_3dpa", ident.2 = "HSC_mock", verbose = FALSE)
HSC.4 <- FindMarkers(zf_de, ident.1 = "HSC_4dpa", ident.2 = "HSC_mock", verbose = FALSE)
HSC.7 <- FindMarkers(zf_de, ident.1 = "HSC_7dpa", ident.2 = "HSC_mock", verbose = FALSE)
# -------------- DE on endothelial cells --------------- #
end.0 <- FindMarkers(zf_de, ident.1 = "End_0dpa", ident.2 = "End_mock", verbose = FALSE)
end.1 <- FindMarkers(zf_de, ident.1 = "End_1dpa", ident.2 = "End_mock", verbose = FALSE)
end.2 <- FindMarkers(zf_de, ident.1 = "End_2dpa", ident.2 = "End_mock", verbose = FALSE)
end.3 <- FindMarkers(zf_de, ident.1 = "End_3dpa", ident.2 = "End_mock", verbose = FALSE)
end.4 <- FindMarkers(zf_de, ident.1 = "End_4dpa", ident.2 = "End_mock", verbose = FALSE)
end.7 <- FindMarkers(zf_de, ident.1 = "End_7dpa", ident.2 = "End_mock", verbose = FALSE)

