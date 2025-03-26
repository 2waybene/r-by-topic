##==============================================================================
##  File:   scRNAseqAddCellLabel.R
##  Author: Jianying Li
##  History: 03/26/2025
##  Comment: credit goes to google search, all about Idents
##  
##===============================================================================

suppressPackageStartupMessages(library(Seurat, quietly =TRUE))
suppressPackageStartupMessages(library(dplyr, quietly =TRUE))
suppressPackageStartupMessages(library(scran, quietly =TRUE))
suppressPackageStartupMessages(library(BiocSingular, quietly =TRUE))
suppressPackageStartupMessages(library(plotly, quietly = TRUE))
suppressPackageStartupMessages(library(patchwork,quietly = TRUE))
suppressPackageStartupMessages(library(tidyverse,quietly = TRUE))

##=================================================
##  processed data
##=================================================
## This is a big file and I can only use googleDrive
## EC2 has 1 GB limitation!! 


Merged.obj <- readRDS("X:/project2025/FelipeHT22Project/preprocessedData/Merged.obj_dimDuc_WTonly.RDS")
DimPlot(Merged.obj, reduction = 'umap', label = TRUE , repel = TRUE)

seurat_obj <- Merged.obj
Idents(seurat_obj) 
##=================================
##  sample code stars
##=================================

# Assuming you have a Seurat object named 'pbmc'

# Get current identities
current_idents <- Idents(object = pbmc)

# Set identities to a metadata column named 'orig.ident'
Idents(object = pbmc) <- "orig.ident"

# Set identities for cells 1:10 to "CD4 T cells"
Idents(object = pbmc, cells = 1:10) <- "CD4 T cells"

# Rename the identity class "CD4 T cells" to "T Helper cells"
pbmc <- RenameIdents(object = pbmc, `CD4 T cells` = "T Helper cells")

##=================================
##  sample code ends
##=================================


head(seurat_obj[[]])
# Rename cell identities


Idents(seurat_obj) <- seurat_obj@meta.data$seurat_clusters
Idents(seurat_obj) 
seurat_obj <- RenameIdents(seurat_obj,
             '0' = "Microglia" ,'1' = "Astrocytes", '2' = "Microglia" ,
             '3' = "Oligodendrocytes", '4' = "Oligodendrocytes", '6' = "Astrocytes" ,
             '7' = "VECs" , '8' = "Astrocytes" , '9'="NPCs",'10' = "Oligodendrocytes", '11' = "Neuron",
             '12' = "Oligodendrocytes", '13' = "Oligodendrocytes", 
             '14' = "BAM", '15'="OPCs", '16' = "Neuron",
             '17' = "Pericytes",
             '19' = "Oligodendrocytes_OPCs","20" = "GABAergic_neuron", '21' = "VECs"  
)
Idents(seurat_obj) 

head(seurat_obj[[]])
new_metadata <- data.frame(
  cell_type = Idents(object = seurat_obj),
  row.names = colnames(seurat_obj) # Ensure rownames match cell names
)

seurat_obj <- AddMetaData(seurat_obj, new_metadata, col.name = "cell_type")
head(seurat_obj[[]])

##  Now, I have successfully add the cell type into the "metadata"
