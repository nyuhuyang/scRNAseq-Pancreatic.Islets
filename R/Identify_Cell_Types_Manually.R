library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#====== 2.1 Creat Gene Panel From SingleR ==========================================
(load(file='../SingleR/data/mouse.rnaseq.RData'))
(load(file = "data/Islets_Harmony_12_20181216.Rda"))
mouse.rnaseq_main = CreatGenePanelFromSingleR(object = mouse.rnaseq,
                                              main.type = TRUE, species = "Mouse")
write.csv(mouse.rnaseq_main,file = "../SingleR/output/mouse.rnaseq_main.csv")

#====== 2.2 marker gene analysis ==========================================
# combine marker gene table
mouse.rnaseq_main = read.csv("../SingleR/output/mouse.rnaseq_main.csv",
                             row.names =1, header = T, stringsAsFactors = F)
Pancreatic.markers <- readxl::read_excel("doc/Pancreatic.cells.markers.xlsx")
colnames(Pancreatic.markers) <- gsub("\\s","_",colnames(Pancreatic.markers))

mouse.pancreatic <- merge(t(Pancreatic.markers), t(mouse.rnaseq_main),
                          by="row.names",all=TRUE)
rownames(mouse.pancreatic) <- mouse.pancreatic$Row.names
mouse.pancreatic <- mouse.pancreatic[,-which(colnames(mouse.pancreatic) == "Row.names")]
mouse.pancreatic.list <- df2list(t(mouse.pancreatic))


marker.list <- lapply(mouse.pancreatic.list, function(x) {
    MouseGenes(object = object, marker.genes= x, unique = T)
})
marker.df <- list2df(marker.list)
marker.list <- df2list(marker.df[1:9,])
str(marker.list)

FeaturePlot.1 <- function(object = object, x){
    p <- FeaturePlot(object = object, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, do.return =T,
                     cols.use = c("lightgrey","blue"), pt.size = 0.5)
    return(p)
}

dev.off() # to accelerate 
for(i in 1:length(marker.list)){
    p <- FeaturePlot.1(object = object, x = marker.list[[i]])
    p1 <- do.call(plot_grid, p)
    p1 <- p1 + ggtitle(paste(names(marker.list)[i],"markers"))+
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
    jpeg(paste0(path,names(marker.list)[i],".jpeg"),
         units="in", width=10, height=7,res=600)
    print(p1)
    print(paste0(i,":",length(marker.list)))
    dev.off()
}

#====== 2.3 manually label ==========================================
object <- SetAllIdent(object, id = "res.0.6")
table(object@ident)
idents <- as.data.frame(table(object@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Beta_cells",
                     "Beta_cells",
                     "Beta_cells/Delta_cells",
                     "Beta_cells",
                     "Alpha_cells",
                     "Beta_cells",
                     "Alpha_cells",
                     "Beta_cells",
                     "Delta_cells",
                     "Alpha_cells/Delta_cells",
                     "Macrophages/Monocytes",
                     "B_cells",
                     "Endothelial_cells",
                     "Fibroblasts",
                     "Gamma_cells",
                     "Acinar_cells/Delta_cells",
                     "NK_cells/T_cells",
                     "Ductal_cells")
object@ident <- plyr::mapvalues(x = object@ident,
                              from = old.ident.ids,
                              to = new.cluster.ids)
object <- StashIdent(object, save.name = "manual")
#Finally, we can also view the labeling as a table compared to the original identities:

# cell number
kable(table(object@ident, object@meta.data$orig.ident)) %>%
    kable_styling()

# cell percentage
prop.table(x = table(object@ident,
                     object@meta.data$orig.ident),margin = 2) %>%
    kable %>% kable_styling

object@meta.data$orig.ident %>% table() %>% kable() %>% kable_styling()
object@meta.data$manual %>% table() %>% kable() %>% kable_styling()

##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
manual.color = as.vector(singler_colors$manual.color[!is.na(singler_colors$manual.color)])
manual.color[duplicated(manual.color)]
length(manual.color)
object@meta.data[,c("manual")] %>% unique %>% length
object@meta.data[,c("manual")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaColor(object = object, label= "manual", colors = manual.color)
object <- SetAllIdent(object = object, id = "manual")
TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F)

p4 <- TSNEPlot.1(object = object, do.label = F, group.by = "ident", 
                 do.return = TRUE, no.legend = F, 
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 4,force = 3)+
    ggtitle("Manually label clusters by marker genes")+
    theme(text = element_text(size=10),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_manual_legend.jpeg"), units="in", width=10, height=7,
     res=600)
print(p4)
dev.off()

##############################
# subset Seurat
###############################
table(object@meta.data$orig.ident)
table(object@ident)

df_samples <- readxl::read_excel("doc/20181204_sample_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
tests <- paste0("test",1:2)
for(test in tests){
    sample_n = which(df_samples$tests %in% test)
    print(samples <- unique(df_samples$sample[sample_n]))
    
    cell.use <- rownames(object@meta.data)[object@meta.data$old.ident %in% samples]
    subset.object <- SubsetData(object, cells.use = cell.use)
    
    g <- SplitTSNEPlot(subset.object,group.by = "ident",split.by = "conditions",
                       select.plots = c(2,1),
                       no.legend = T,do.label =F,label.size=3,
                       return.plots =T, label.repel = T,force=2)
    jpeg(paste0(path,test,"_conditions_TSNEPlot.jpeg"), units="in", width=10, height=7,
         res=600)
    print(do.call(plot_grid, g))
    dev.off()
}


#======== bar chart cell number ==============
object@meta.data$orig.ident = gsub("GIPR-","",object@meta.data$orig.ident)
Split_object <- SplitSeurat(object, split.by = "tests")
object@meta.data$tests %>% unique %>% sort
Freq <- table(Split_object[[2]]@meta.data$manual, 
              Split_object[[2]]@meta.data$orig.ident)

jpeg(paste0(path,"/bar_chart1.jpeg"), units="in", width=10, height=7,
     res=600)
par(mfrow=c(1, 1), mar=c(5, 4, 4, 0))
barplot(Freq, main="Total numbers of each male samples",
        xlab="samples", ylab="Cell numbers",
        col = manual.color)
        #legend = rownames(Freq),
        #args.legend = list(x = "topleft", bty = "n",inset=c(0, -0.18)))
dev.off()

#======== bar chart cell percentage ==============
Split_object <- SplitSeurat(object, split.by = "tests")
object@meta.data$tests %>% unique %>% sort

# cell percentage
Freq <- table(Split_object[[3]]@meta.data$manual, 
              Split_object[[3]]@meta.data$orig.ident) %>%
    prop.table(margin = 2)

jpeg(paste0(path,"/bar_p_chart2.jpeg"), units="in", width=10, height=7,
     res=600)
par(mfrow=c(1, 1), mar=c(5, 4, 4, 0))
barplot(Freq, main="Total numbers of each female samples",
        xlab="samples", ylab="Cell numbers",
        col = manual.color)
dev.off()
