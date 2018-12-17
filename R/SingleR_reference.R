############################################
# combine mouse.rnaseq and 
############################################
library(SingleR)
library(genefilter)
library(dplyr)
library(magrittr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
####functions===========

data("mouse.rnaseq")
head(mouse.rnaseq$types)
dim(mouse.rnaseq$data)
head(mouse.rnaseq$data[,1:5])
anyNA(mouse.rnaseq$data)
testMMM(mouse.rnaseq$data)
boxplot(mouse.rnaseq$data, main="mouse.rnaseq")#slow!

# remove low quanlity mouse.rnaseq data
par(mfrow=c(2,1))
hist(colMeans(mouse.rnaseq$data),breaks=ncol(mouse.rnaseq$data))
quantile_75 <- apply(mouse.rnaseq$data,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(mouse.rnaseq$data))
rm_samples <- names(quantile_75)[quantile_75<5]
(rm_index <- which(colnames(mouse.rnaseq$data) %in% rm_samples))
mouse.rnaseq_rm <- mouse.rnaseq$data[,-rm_index]
par(mfrow=c(1,1))
boxplot(mouse.rnaseq_rm,main="mouse.rnaseq_rm")#slow!
testMMM(mouse.rnaseq_rm)


sort(unique(mouse.rnaseq$main_types))
types = FineTune(as.character(mouse.rnaseq$types[-rm_index]),
                       main.type = FALSE)
main_types = FineTune(as.character(mouse.rnaseq$main_types[-rm_index]),
                            main.type = TRUE)
sort(unique(main_types))
mouse.rnaseq = CreateSinglerReference(name = 'mouse.rnaseq',
                                          expr = mouse.rnaseq_rm,
                                          types = types, 
                                          main_types = main_types)
save(mouse.rnaseq,file='../SingleR/data/mouse.rnaseq.RData')

# check GSE80673 data==============================
GSE80673 <- read.delim2("data/GSE80673_Huising_mouse_deltaVSbetaVSalpha_STAR_mm10_RPKM.txt")
dim(GSE80673)

# remove NA rows
GSE80673 <- GSE80673[!apply(GSE80673,1, function(x) all(is.na(x))),]

#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
RemoveDup <- function(mat){
        gene_id <- as.matrix(mat[,1])
        mat <- mat[,-1]
        if(!is.matrix(mat)) mat <- sapply(mat,as.numeric)
        rownames(mat) <- 1:nrow(mat)
        mat[is.na(mat)] = 0
        mat <- cbind(mat, "rowSums" = rowSums(mat))
        mat <- mat[order(mat[,"rowSums"],decreasing = T),]
        gene_id <- gene_id[as.numeric(rownames(mat))]
        remove_index <- duplicated(gene_id)
        mat <- mat[!remove_index,]
        rownames(mat) <- gene_id[!remove_index]
        return(mat[,-ncol(mat)])
}
GSE80673 <- RemoveDup(GSE80673)
dim(GSE80673)
head(GSE80673)
testMMM(GSE80673)

# merge MCL and mouse.rnaseq

(load(file='../SingleR/data/mouse.rnaseq.RData'))
dim(mouse.rnaseq$data)
GSE80673_mouse.rnaseq <- merge(log1p(GSE80673),mouse.rnaseq$data,
                         by="row.names",all=FALSE)
rownames(GSE80673_mouse.rnaseq) = GSE80673_mouse.rnaseq$Row.names
GSE80673_mouse.rnaseq <- GSE80673_mouse.rnaseq[-which(colnames(GSE80673_mouse.rnaseq)=="Row.names")]
testMMM(GSE80673_mouse.rnaseq)

colsum <- colSums(GSE80673_mouse.rnaseq)
scale_factor = median(colsum)
GSE80673_mouse.rnaseq = GSE80673_mouse.rnaseq/colsum * scale_factor
testMMM(GSE80673_mouse.rnaseq)


jpeg(paste0(path,"boxplot_GSE80673_mouse.rnaseq.jpeg"), units="in", width=10, height=7,res=600)
boxplot(GSE80673_mouse.rnaseq) #slow
title(main = "boxplot for GSE80673 + mouse.rnaseq")
dev.off()



# GSE80673 cell types======
GSE80673.rnaseq <- readxl::read_excel("doc/20181204_sample_info.xlsx",
                      sheet = "GSE80673")
(main_types <- match(colnames(GSE80673),
                            paste0("X",GSE80673.rnaseq$sample.id)) %>%
        GSE80673.rnaseq$main_types[.])

(types <- match(colnames(GSE80673),
                     paste0("X",GSE80673.rnaseq$sample.id)) %>%
                GSE80673.rnaseq$types[.])

# Create Singler Reference=============
ref = CreateSinglerReference(name = 'GSE80673_mouse.rnaseq',
                             expr = as.matrix(GSE80673_mouse.rnaseq), # the expression matrix
                             types = c(types,mouse.rnaseq$types), 
                             main_types = c(main_types, mouse.rnaseq$main_types))

save(ref,file='data/ref_GSE80673_mouse.rnaseq.RData') # it is best to name the object and the file with the same name.
