library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/Islets_Harmony_12_20181216.Rda"))
(load(file='../SingleR/data/immgen.rda'))
(load(file ='data/ref_GSE80673_mouse.rnaseq.RData'))

object@scale.data = NULL; GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();

singler = CreateSinglerObject(object@data, annot = NULL, project.name=object@project.name,
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              ref.list = list(immgen, ref),normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL,
                              numCores = SingleR.numCores/2)

# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object@data)
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = object@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_Islets_12F_20181216.RData")