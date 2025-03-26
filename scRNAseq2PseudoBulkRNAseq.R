##==============================================================================
##  File:   scRNAseq2PseudoBulkRNAseq.R
##  Author: Jianying Li
##  History: 03/24/2025
##  Comment: 
##  X:\project2025\FelipeHT22Project\rScripts\GSE242666_SeuratAnalysis_WT_only.rmd
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

##=================================================
##  label the cell population
##=================================================
Merged.obj <- RenameIdents(Merged.obj,
                           '0' = "Microglia" ,'1' = "Astrocytes", '2' = "Microglia" ,
                           '3' = "Oligodendrocytes", '4' = "Oligodendrocytes", '6' = "Astrocytes" ,
                           '7' = "VECs" , '8' = "Astrocytes" , '9'="NPCs",'10' = "Oligodendrocytes", '11' = "Neuron",
                           '12' = "Oligodendrocytes", '13' = "Oligodendrocytes", 
                           '14' = "BAM", '15'="OPCs", '16' = "Neuron",
                           '17' = "Pericytes",
                           '19' = "Oligodendrocytes_OPCs","20" = "GABAergic_neuron", '21' = "VECs"  
)

## missing "5" and "18"
DimPlot(Merged.obj, reduction = 'umap', label = TRUE , repel = TRUE)
head(Merged.obj[[]])

##=========================================================
##  examine the cell count by population distribution
##=========================================================

cell.count <- as.data.frame(table(Merged.obj@active.ident))
colnames(cell.count)[1] <- "Population"


countsMX <- GetAssayData(object = Merged.obj, slot = "counts", assay="RNA")


##======================================================================
##  Here, I am creating a ScaleFactor to compensate the "read depth" 
##  unblance among different cell populations
##  right or wrong?? I don't know
##======================================================================
cell.count <- cbind(cell.count, as.vector(list ("ScaleFactor" = sum(cell.count$Freq)/cell.count$Freq)))

sum(cell.count$Freq)
freq = as.numeric(cell.count$Freq)/sum(cell.count$Freq)

cell.count <- cbind(cell.count, as.data.frame(list(RelativeFreq = freq)))

ggplot(cell.count, aes(x = Population, y = RelativeFreq, fill = Population)) + geom_col() + coord_flip()


##====================================================
## Created cell specific data matrix for WT only
##====================================================
## Astrocytes.srat <- subset (Merged.obj, idents = "Astrocytes")
srat <- subset (Merged.obj, idents = "Astrocytes")
counts <- GetAssayData(object = srat, slot = "counts", assay="RNA")
metadata <- srat@meta.data
counts <- counts[ ,which(metadata$orig.ident == "WT")]
metadata <- metadata[which(metadata$orig.ident == "WT"),]
gene.cnt.collapse = as.vector(rowSums(counts))
names(gene.cnt.collapse) <- row.names(counts)





##=================================================
##  create data matrix, there is NO replicate
##=================================================

selected.populations <-  c( "Microglia" ,  "Oligodendrocytes", 
                            "VECs", "NPCs", "Neuron","BAM", "OPCs",
                            "Pericytes", "Oligodendrocytes_OPCs", "GABAergic_neuron" )

##  stupid but is working
##  get one population first "Astrocytes", then loop through other populations
##  and get them all

srat <- subset (Merged.obj, idents = "Astrocytes")
counts <- GetAssayData(object = srat, slot = "counts", assay="RNA")
metadata <- srat@meta.data
counts <- counts[ ,which(metadata$orig.ident == "WT")]
metadata <- metadata[which(metadata$orig.ident == "WT"),]
gene.cnt.collapse = as.vector(rowSums(counts))
names(gene.cnt.collapse) <- row.names(counts)
f <- cell.count$ScaleFactor[cell.count$Population == "Astrocytes"]
dm <- as.data.frame(ceiling(gene.cnt.collapse*f))
colnames(dm) <- "Astrocytes"



for (i in selected.populations)
{
  
  srat <- subset (Merged.obj, idents = i)
  counts <- GetAssayData(object = srat, slot = "counts", assay="RNA")
  metadata <- srat@meta.data
  counts <- counts[ ,which(metadata$orig.ident == "WT")]
  metadata <- metadata[which(metadata$orig.ident == "WT"),]
  gene.cnt.collapse = as.vector(rowSums(counts))
  names(gene.cnt.collapse) <- row.names(counts)
  f <- cell.count$ScaleFactor[cell.count$Population == i]
  dm.temp <- as.data.frame(ceiling(gene.cnt.collapse*f))
  #colnames(dm) <- "Astrocytes"
  #dm.temp <- as.data.frame(gene.cnt.collapse)
  colnames(dm.temp) = i
  dm <- cbind(dm, dm.temp)
}

dim(dm)
head(dm)
is.na(dm)
write.csv(dm, file ="X:/project2025/FelipeHT22Project/preprocessedData/GSE242666_WT_count_by_population.csv")
write.table(dm, file ="X:/project2025/FelipeHT22Project/preprocessedData/GSE242666_WT_count_by_population_4_decon.txt", sep="\t")

