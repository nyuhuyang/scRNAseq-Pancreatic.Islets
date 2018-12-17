#======1.1.1 Load the bulk RNA-data files and Set up Seurat object =========================
# read sample summary list
df_bulk <- read_excel("doc/20181204_sample_info.xlsx",
                      sheet = "bulk_sample_info")
bulk_counts <- read.delim("data/bulk_data/bulk_counts.txt",
                          row.names =1)
gene_annotation <- read.delim("data/bulk_data/gene_annotation.txt")
rownames(bulk_counts) %<>% match(gene_annotation$ID) %>% 
        gene_annotation[.,"symbol"]

colnames(bulk_counts) %<>% match(paste0("X",df_bulk$SampleID)) %>%
        df_bulk$sample_label[.]

bulk_object <- CreateSeuratObject(bulk_counts,
                                   min.cells = 3,
                                   min.genes = 100,
                                   names.delim = "_")
orig.ident <- colnames(bulk_counts) %>% toupper %>% 
        gsub("-MALE-","-M-",.) %>%
        gsub("-FEMALE-","-F-",.) %>%
        gsub("WT-","GIPR-WT-",.)

bulk_object@meta.data$orig.ident <- orig.ident
bulk_object@meta.data$tests <- "test0"
bulk_object@meta.data$conditions <- gsub("-[0-9]$","",orig.ident)


#======1.4 mito, QC, filteration =========================
all.mito.genes <- read.csv("./doc/Mitochondrial.csv",row.names = 1) %>% rownames %>%
        tolower %>% Hmisc::capitalize()
(mito.genes <-  rownames(x = bulk_object@data) %in% all.mito.genes %>%
                rownames(x = bulk_object@data)[.])
percent.mito <- Matrix::colSums(bulk_object@raw.data[mito.genes, ])/Matrix::colSums(bulk_object@raw.data)
bulk_object <- AddMetaData(bulk_object, metadata = percent.mito, col.name = "percent.mito")

bulk_object <- SetAllIdent(bulk_object, id = "orig.ident")


bulk_object@ident = factor(bulk_object@ident,levels = samples)
g3 <- VlnPlot(bulk_object, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 5,size.x.use = 10, group.by = "ident",
              x.lab.rot = T, do.return = T,return.plotlist =T)
save(g3,file = "./data/g2_12_20181215.Rda")
jpeg(paste0(path,"/S2_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(g3[[1]]+ggtitle("Total gene number in bulk RNA-seq data"))
dev.off()
jpeg(paste0(path,"/S2_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(g3[[2]]+ggtitle("Total counts in bulk RNA-seq data"))
dev.off()
jpeg(paste0(path,"/S2_mito.jpeg"), units="in", width=10, height=7,res=600)
print(g3[[3]]+ggtitle("mito % in bulk RNA-seq data"))
dev.off()

save(bulk_object, file = "./data/Islets_bulk_12_20181215.Rda")
