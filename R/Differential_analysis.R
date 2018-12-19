library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#====== 4.1 Subset ==========================================
(load(file = "data/Islets_Harmony_12_20181216.Rda"))

object <- SetAllIdent(object, id = "tests")
# It looks like fewer beta cells express E350Q variant compared to WT GIPR.
# This doesn’t seem to be the case for the males.  Is this true?
object_female <- SubsetData(object, ident.use = "test2")
markers <- "Gipr"
SplitSingleFeaturePlot(object_female,group.by = "ident",split.by = "conditions",
                       select.plots = c(2,1),
                       no.legend = T,label.size=3,do.print =T,markers = markers,
                       threshold = 0.1)
# cell percentage
prob <- prop.table(x = table(object_female@data[markers,]>0,
                     object_female@meta.data$orig.ident),margin = 2) 
prob %>% kable()  %>% kable_styling
prob <- as.data.frame(prob)

t.test(prob[as.logical(prob$Var1),"Freq"][1:3],
       prob[as.logical(prob$Var1),"Freq"][4:6])

#====Male
object_male <- SubsetData(object, ident.use = "test1")
SplitSingleFeaturePlot(object_male,group.by = "ident",split.by = "conditions",
                       select.plots = c(2,1),
                       no.legend = T,label.size=3,do.print =T,markers = markers,
                       threshold = 0.1)
# cell percentage
prob <- prop.table(x = table(object_male@data[markers,]>0,
                             object_male@meta.data$orig.ident),margin = 2) 
prob %>% kable()  %>% kable_styling
prob <- as.data.frame(prob)

t.test(prob[as.logical(prob$Var1),"Freq"][1:3],
       prob[as.logical(prob$Var1),"Freq"][4:6])

#######################
# geom_density
#######################
markers <- "Gipr"
genders <- c("Male", "Female")
df_samples <- readxl::read_excel("doc/20181204_sample_info.xlsx")
(conditions <- df_samples$Conditions %>% unique)
g <- list()

object <- SubsetData(object, ident.remove = "test0")
for(i in 1:2){
        gender <- conditions[c(i*2-1,i*2)]
        index <- object@meta.data$conditions %in% gender
        data.use <- split(rownames(object@meta.data[index,]),
                          object@meta.data[index,"orig.ident"]) %>% 
                lapply(function(cells_use) {
                        object@data[markers,cells_use] %>% 
                        as.matrix %>% t %>% as.data.frame %>%
                        gather(key = barcode, value = ave.expr) %>%
                        mutate(samples = gsub('\\_.*', '',.[,"barcode"])) %>%
                        mutate(conditions = gsub('-1|-2|-3', '',.[,"samples"]))
                        }) %>% bind_rows
        g[[i]] <- ggplot(data.use, aes(x = ave.expr, color = conditions)) + 
                geom_density(size = 1) +
                scale_y_sqrt() + ylim(0, 1)+
                xlab(paste(markers,"expression (log nUMI)"))+
                ggtitle(genders[i])+
                theme(text = element_text(size=15),
                      #legend.position="none", 
                      legend.position=c(0.3,0.85) ,
                      plot.title = element_text(hjust = 0.5,size = 15,
                                                face = "bold"))
                }
jpeg(paste0(path,"/density_plot.jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid,g))
dev.off()
