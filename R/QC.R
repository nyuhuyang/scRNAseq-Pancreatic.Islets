########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20181204_sample_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("test",paste0("test",1:2)))
df_samples %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample.id[sample_n]
conditions <- df_samples$conditions[sample_n]

# check missing data
(current <- list.files("data/single_cell_data")[!grepl(".Rda|RData",list.files("data"))])
(missing_data <- sample.id[!(sample.id %in% current)])

if(length(missing_data)>0){
        # Move files from Download to ./data and rename them
        for(missing_dat in missing_data){
                old.pth  <- paste("~/Downloads", missing_dat,"outs",
                                  "filtered_gene_bc_matrices","mm10",
                                  sep = "/")
                list.of.files <- list.files(old.pth)
                new.sample <- gsub("_","-",missing_dat)
                new.folder <- paste("./data", missing_dat,"outs",
                                    "filtered_gene_bc_matrices",
                                    "mm10",sep = "/")
                if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
                # copy the files to the new folder
                file.copy(paste(old.pth, list.of.files, sep = "/"), new.folder)
        }
}

## Load the dataset
object_raw <- list()
object_Seurat <- list()
for(i in 1:length(samples)){
        object_raw[[i]] <- Read10X(data.dir = paste0("./data/single_cell_data/",
                                                     sample.id[i]))
        colnames(object_raw[[i]]) <- paste0(samples[i],"_",colnames(object_raw[[i]]))
        object_Seurat[[i]] <- CreateSeuratObject(object_raw[[i]],
                                               min.cells = 3,
                                               min.genes = 100,
                                               names.delim = "_")
        object_Seurat[[i]]@meta.data$conditions <- conditions[i]
}

#======1.1.2 QC before merge =========================
cell.number <- sapply(object_Seurat, function(x) length(x@cell.names))
QC_list <- lapply(object_Seurat, function(x) as.matrix(x = x@raw.data))
median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))

min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))

QC.list <- cbind(df_samples,cell.number, median.nUMI,median.nGene,min.nUMI,min.nGene,
                 row.names = samples)
write.csv(QC.list,paste0(path,"QC_list.csv"))
QC.list %>% kable() %>% kable_styling()
remove(QC_list,median.nUMI,median.nGene,min.nUMI,min.nGene,QC.list);GC()

#========1.1.3 merge ===================================
object <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), object_Seurat)
remove(object_raw,object_Seurat);GC()

# read and select mitochondial genes
#all.mito.genes <- read.csv("./doc/Mitochondrial.csv",row.names = 1)
#grep("Mitochondria",all.mito.genes$Description) %>% all.mito.genes[.,] %>% 
#        write.csv(file = "./doc/Mitochondrial.csv")
all.mito.genes <- read.csv("./doc/Mitochondrial.csv",row.names = 1) %>% rownames %>%
        tolower %>% Hmisc::capitalize()
(mito.genes <-  rownames(x = object@data) %in% all.mito.genes %>%
        rownames(x = object@data)[.])
percent.mito <- Matrix::colSums(object@raw.data[mito.genes, ])/Matrix::colSums(object@raw.data)
object <- AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito")
object@ident = factor(object@ident,levels = samples)
g1 <- VlnPlot(object = object, features.plot = c("nGene", "nUMI", "percent.mito"),
              nCol = 1,point.size.use = 0.2,,size.x.use = 10, group.by = "ident",
              x.lab.rot = T, do.return = T,return.plotlist =T)

save(g1,file="./data/g1_12_20181215.Rda")
