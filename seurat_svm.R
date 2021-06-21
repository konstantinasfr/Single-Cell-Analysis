#### INTEGRATED #########################
BiocManager::install('multtest')
install.packages('metap')
library(metap)
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)

###### Importing the data ###################################
sham1 <- read.csv('Dataset/GSM4376680_Sham_1day_plate1_UMIcount.tsv', sep='\t',header=TRUE)
rownames(sham1) <- sham1$GENEID
sham1 <- sham1[grep("ERCC",rownames(sham1),invert=TRUE),-1]
new_col_names <- c()
for (i in 1:384){
  new_col_names<-c(new_col_names,paste0("sham1_", i))
}
colnames(sham1)<- new_col_names

sham2 <- read.csv('Dataset/GSM4376681_Sham_1day_plate2_UMIcount.tsv', sep='\t',header=TRUE)
rownames(sham2) <- sham2$GENEID
sham2 <- sham2[grep("ERCC",rownames(sham2),invert=TRUE),-1]
new_col_names <- c()
for (i in 1:384){
  new_col_names<-c(new_col_names,paste0("sham2_", i))
}
colnames(sham2)<- new_col_names

sham3 <- read.csv('Dataset/GSM4376682_Sham_1day_plate3_UMIcount.tsv', sep='\t',header=TRUE)
rownames(sham3) <- sham3$GENEID
sham3 <- sham3[grep("ERCC",rownames(sham3),invert=TRUE),-1]
new_col_names <- c()
for (i in 1:384){
  new_col_names<-c(new_col_names,paste0("sham3_", i))
}
colnames(sham3)<- new_col_names

sham4 <- read.csv('Dataset/GSM4376683_Sham_1day_plate4_UMIcount.tsv', sep='\t',header=TRUE)
rownames(sham4) <- sham4$GENEID
sham4 <- sham4[grep("ERCC",rownames(sham4),invert=TRUE),-1]
new_col_names <- c()
for (i in 1:384){
  new_col_names<-c(new_col_names,paste0("sham4_", i))
}
colnames(sham4)<- new_col_names

mi1 <- read.csv('Dataset/GSM4614787_MI_3day_plate1_UMIcount.tsv', sep='\t',header=TRUE)
rownames(mi1) <- mi1$GENEID
mi1 <- mi1[grep("ERCC",rownames(mi1),invert=TRUE),-1]
new_col_names <- c()
for (i in 1:384){
  new_col_names<-c(new_col_names,paste0("mi1_", i))
}
colnames(mi1)<- new_col_names

mi2 <- read.csv('Dataset/GSM4614788_MI_3day_plate2_UMIcount.tsv', sep='\t',header=TRUE)
rownames(mi2) <- mi2$GENEID
mi2 <- mi2[grep("ERCC",rownames(mi2),invert=TRUE),-1]
new_col_names <- c()
for (i in 1:384){
  new_col_names<-c(new_col_names,paste0("mi2_", i))
}
colnames(mi2)<- new_col_names

mi3 <- read.csv('Dataset/GSM4614789_MI_3day_plate3_UMIcount.tsv', sep='\t',header=TRUE)
rownames(mi3) <- mi3$GENEID
mi3 <- mi3[grep("ERCC",rownames(mi3),invert=TRUE),-1]
new_col_names <- c()
for (i in 1:384){
  new_col_names<-c(new_col_names,paste0("mi3_", i))
}
colnames(mi3)<- new_col_names

mi4 <- read.csv('Dataset/GSM4614790_MI_3day_plate4_UMIcount.tsv', sep='\t',header=TRUE)
rownames(mi4) <- mi4$GENEID
mi4 <- mi4[grep("ERCC",rownames(mi4),invert=TRUE),-1]
new_col_names <- c()
for (i in 1:384){
  new_col_names<-c(new_col_names,paste0("mi4_", i))
}
colnames(mi4)<- new_col_names

mi5 <- read.csv('Dataset/GSM4614791_MI_3day_plate5_UMIcount.tsv', sep='\t',header=TRUE)
rownames(mi5) <- mi5$GENEID
mi5 <- mi5[grep("ERCC",rownames(mi5),invert=TRUE),-1]
new_col_names <- c()
for (i in 1:384){
  new_col_names<-c(new_col_names,paste0("mi5_", i))
}
colnames(mi5)<- new_col_names

tmp12 <- merge(sham1, sham2, by=0, all=F)
rownames(tmp12) <- tmp12$Row.names; tmp12$Row.names <- NULL
tmp123 <- merge(tmp12, sham3, by=0, all=F)
rownames(tmp123) <- tmp123$Row.names; tmp123$Row.names <- NULL
sham_total <- merge(tmp123, sham4, by=0, all=F)
rownames(sham_total) <- sham_total$Row.names; sham_total$Row.names <- NULL

tmp12 <- merge(mi1, mi2, by=0, all=F)
rownames(tmp12) <- tmp12$Row.names; tmp12$Row.names <- NULL
tmp123 <- merge(tmp12, mi3, by=0, all=F)
rownames(tmp123) <- tmp123$Row.names; tmp123$Row.names <- NULL
tmp1234 <- merge(tmp123, mi4, by=0, all=F)
rownames(tmp1234) <- tmp1234$Row.names; tmp1234$Row.names <- NULL
mi_total <- merge(tmp1234, mi5, by=0, all=F)
rownames(mi_total) <- mi_total$Row.names; mi_total$Row.names <- NULL


###################### sham prepossessing ##################################

sham <- CreateSeuratObject(counts = sham_total, project = "sham", min.cells = 3, min.features = 200)
# filtering process and QC
sham[["percent.mt"]] <- PercentageFeatureSet(sham, pattern = "^mt-")

VlnPlot(sham, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#sham <- subset(sham, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <25)
sham <- subset(sham, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 25 & nCount_RNA >800)

sham = sham@assays[["RNA"]]@counts
nonmito_genes <- grep(pattern = "chrM", x = rownames(x = sham), value = TRUE, invert = TRUE)
sham = sham[nonmito_genes,]

sham <- CreateSeuratObject(counts = sham, project = "sham", min.cells = 3, min.features = 200)


sham <- NormalizeData(sham, normalization.method = "LogNormalize", scale.factor = 10000)
sham <- FindVariableFeatures(sham, selection.method = "vst", nfeatures = 2000)
sham@meta.data[, "protocol"] <- "sham"
###################### mi prepossessing ####################################

mi <- CreateSeuratObject(counts = mi_total, project = "mi", min.cells = 3, min.features = 200)

# filtering process and QC
mi[["percent.mt"]] <- PercentageFeatureSet(mi, pattern = "^mt-")




VlnPlot(mi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

mi <- subset(mi, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt <25 & nCount_RNA >800)

mi = mi@assays[["RNA"]]@counts
nonmito_genes <- grep(pattern = "chrM", x = rownames(x = mi), value = TRUE, invert = TRUE)
mi = mi[nonmito_genes,]

mi <- CreateSeuratObject(counts = mi, project = "mi", min.cells = 3, min.features = 200)

mi <- NormalizeData(mi, normalization.method = "LogNormalize", scale.factor = 10000)
mi <- FindVariableFeatures(mi, selection.method = "vst", nfeatures = 2000)
mi@meta.data[, "protocol"] <- "mi"

#################### Integration #####################
combined = merge(sham, y = mi, add.cell.ids = c("sham", "mi"), project = "protocol")
list <- SplitObject(combined, split.by = "protocol")
features <- SelectIntegrationFeatures(object.list = list,nfeatures = 2000)
reference.list <- list[c("sham", "mi")]
anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:8, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors, dims = 1:8)


attributes(FindNeighbors)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:10)
ElbowPlot(integrated)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:10, k.param = 20)
integrated <- FindClusters(integrated, resolution = 0.08)
table(Idents(integrated))

################### HEREEEE ################################
#save(integrated, file="integrated.rds")

#integrated_seurat <- SetIdent(object = integrated,cells = cardio_sham_seurat,value = 'cardio_sham')

#l<-which(integrated@active.ident==3 & integrated@meta.data[["protocol"]]=="mi")
ee<- which(integrated@active.ident==0)
ee1<- which(integrated@active.ident==1)
ee2<- which(integrated@active.ident==2)
ee3<- which(integrated@active.ident==3)
integrated <- SetIdent(object = integrated, cells = ee,value = 'Fibroblasts')
integrated <- SetIdent(object = integrated, cells = ee1,value = 'Cardiomyocytes')
integrated <- SetIdent(object = integrated, cells = ee2,value = 'Endothelial')
integrated <- SetIdent(object = integrated, cells = ee3,value = 'Immune')
# Visualization

p1 <- DimPlot(integrated, reduction = "umap", group.by = "protocol")
p2 <- DimPlot(integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(integrated, reduction = "umap", split.by = "protocol")

saveRDS(integrated, file = "integrated.rds")

############################################################################
# These are the plots of the cells in which the known cell markers are expressed divided by cell type
# These cells markers helped for the creation of the training data
# Cardiomyocytes
features <- c('Ttn--chr2', 'Myh6--chr14', 'Des--chr1', 'Actc1--chr2', 'Csrp3--chr7', 'Tnnt2--chr1', 'Ryr2--chr13')

RidgePlot(object = integrated, features = features, ncol = 5)
VlnPlot(object = integrated, features = features)
FeaturePlot(object = integrated, features = features)

# Fibroblasts
features <- c('Col1a1--chr11', 'Col3a1--chr1', 'Fbln2--chr6', 'Sparc--chr11', 'Fstl1--chr16', 'Mmp2--chr8', 'Vim--chr2')

RidgePlot(object = integrated, features = features, ncol = 5)
VlnPlot(object = integrated, features = features)
FeaturePlot(object = integrated, features = features)

# Endothelial
features <- c('Pecam1--chr11', 'Fabp4--chr3', 'Tie1--chr4', 'Egfl7--chr2', 'Flt1--chr5', 'Epas1--chr17', 'Emcn--chr3')

RidgePlot(object = integrated, features = features, ncol = 5)
VlnPlot(object = integrated, features = features)
FeaturePlot(object = integrated, features = features)

# Immune
features <- c('Il1b--chr2', 'S100a8--chr3', 'Ilr2', 'S100a9--chr3', 'Mmp9', 'Cebpb--chr2', 'Csf3r--chr4')
features <- c('Il1b--chr2', 'S100a8--chr3', 'S100a9--chr3', 'Cebpb--chr2', 'Csf3r--chr4')

RidgePlot(object = integrated, features = features, ncol = 5)
VlnPlot(object = integrated, features = features)
FeaturePlot(object = integrated, features = features)



# Macrophage
#features <- c('Cd68--chr11', 'Lyz1--chr10', 'Lgals3--chr14', 'Itgam', 'Csfr1', 'Emr1--chr17', 'Cd74--chr18')

#RidgePlot(object = integrated, features = features, ncol = 5)
#VlnPlot(object = integrated, features = features)
#FeaturePlot(object = integrated, features = features)

#DefaultAssay(integrated) <- "RNA"
#cm.markers <- FindConservedMarkers(integrated, ident.1 =7, grouping.var = "protocol", verbose = FALSE)
#head(cm.markers,20)
#FeaturePlot(integrated, features, min.cutoff = "q9")

####################################################################################
######### Creation of training set ###########################
# After inspecting the cells that express the known markers we create the training set by isolating only the cells of interest
# Cardiomyocytes 

cardiom1 <- integrated@assays[["RNA"]]@data['Ttn--chr2',]
cardiom2 <- integrated@assays[["RNA"]]@data['Myh6--chr14',]
cardiom3 <- integrated@assays[["RNA"]]@data['Des--chr1',]
cardiom4 <- integrated@assays[["RNA"]]@data['Actc1--chr2',]
cardiom5 <- integrated@assays[["RNA"]]@data['Csrp3--chr7',]
cardiom6 <- integrated@assays[["RNA"]]@data['Tnnt2--chr1',]

cardiom_cells <- which(cardiom1>3 & cardiom2>3 &cardiom3>3 & cardiom4>3 & cardiom5>3 & cardiom6>3 | (cardiom1>5 & cardiom2>5 &cardiom3>5))


length(cardiom_cells)
cardiom_cells_names <- names(cardiom_cells)
pos_cells = subset(integrated,cells=cardiom_cells)
p3 <- FeaturePlot(pos_cells,'Ttn--chr2')+ggtitle("Cardiomyocytes")

x=c()
x[1] = "Cardiomyocytes"
for (i in 2:length(cardiom_cells)) { 
  x[i]="Cardiomyocytes"
}

Cardiomyocytes <- as.data.frame(x)
row.names(Cardiomyocytes) <- cardiom_cells_names

plot_grid(p1, p3)


# Fibroblasts
fibro1 <- integrated@assays[["RNA"]]@data['Col1a1--chr11',]
fibro2 <- integrated@assays[["RNA"]]@data['Col3a1--chr1',]
fibro3 <- integrated@assays[["RNA"]]@data['Fbln2--chr6',]
fibro4 <- integrated@assays[["RNA"]]@data['Sparc--chr11',]
fibro5 <- integrated@assays[["RNA"]]@data['Fstl1--chr16',]
fibro6 <- integrated@assays[["RNA"]]@data['Vim--chr2',]

fibro_cells <- which(fibro1>2.85 & fibro2>2.85 &fibro3>2.85 & fibro4>2.8 & fibro5>2.8 & fibro6>2.8 | (fibro1>5 & fibro2>5 &fibro3>5))

length(fibro_cells)
fibro_cells_names <- names(fibro_cells)
pos_cells = subset(integrated,cells=fibro_cells)
p3 <- FeaturePlot(pos_cells,'Col1a1--chr11')+ggtitle("Fibroblasts")

x=c()
x[1] = "Fibroblasts"
for (i in 2:length(fibro_cells)) { 
  x[i]="Fibroblasts"
}

Fibroblasts <- as.data.frame(x)
row.names(Fibroblasts) <- fibro_cells_names
Card_Fibr <- rbind(Cardiomyocytes, Fibroblasts)

plot_grid(p1, p3)

# Endothelial
# Flt1--chr5', 'Epas1--chr17', 'Emcn--chr3')
endo1 <- integrated@assays[["RNA"]]@data['Pecam1--chr11',]
endo2 <- integrated@assays[["RNA"]]@data['Fabp4--chr3',]
endo3 <- integrated@assays[["RNA"]]@data['Tie1--chr4',]
endo4 <- integrated@assays[["RNA"]]@data['Egfl7--chr2',]

endo_cells <- which(endo1>1.8 & endo2>1.8 &endo3>1.8 & endo4>1.8)

length(endo_cells)
endo_cells_names <- names(endo_cells)
pos_cells = subset(integrated,cells=endo_cells)
p3 <- FeaturePlot(pos_cells,'Fabp4--chr3')+ggtitle("Endothelial")

x=c()
x[1] = "Endothelial"
for (i in 2:length(endo_cells)) { 
  x[i]="Endothelial"
}

Endothelial <- as.data.frame(x)
row.names(Endothelial) <- endo_cells_names
Card_Fibr_Endo <- rbind(Card_Fibr, Endothelial)

plot_grid(p1, p3)


# Immune
# 'Cebpb--chr2'
immun1 <- integrated@assays[["RNA"]]@data['Il1b--chr2',]
immun2 <- integrated@assays[["RNA"]]@data['S100a8--chr3',]
immun3 <- integrated@assays[["RNA"]]@data['S100a9--chr3',]
immun4 <- integrated@assays[["RNA"]]@data['Csf3r--chr4',]

immun_cells <- which(immun1>2.2 & immun2>2.2 &immun3>2.2 & immun4>2.2)

length(immun_cells)
immun_cells_names <- names(immun_cells)
pos_cells = subset(integrated,cells=immun_cells)
p3 <- FeaturePlot(pos_cells,'Il1b--chr2')+ggtitle("Immune")

ee <- GetAssayData(object = pos_cells, 
             assay = "RNA", slot = "data")

x=c()
x[1] = "Immune"
for (i in 2:length(immun_cells)) { 
  x[i]="Immune"
}

Immune <- as.data.frame(x)
row.names(Immune) <- immun_cells_names
Card_Fibr_Endo_Immun <- rbind(Card_Fibr_Endo, Immune)

plot_grid(p1, p3)

### preparing SVM input ########

expr_tot = GetAssayData(integrated)
expr_pcas = as.data.frame(integrated@reductions[["pca"]]@cell.embeddings)

known_cells <- merge(expr_pcas, Card_Fibr_Endo_Immun, by=0, all=F)
rownames(known_cells) <-known_cells[,1]
known_cells[,1] <- NULL #Removing the first column
known_cells[,9:50] <- NULL #Removing from PC9 to PC50
colnames(known_cells)[9] <- "Type"

unknown_cells <- expr_pcas[!(row.names(expr_pcas) %in% immun_cells_names), ]
unknown_cells <- unknown_cells[!(row.names(unknown_cells) %in% cardiom_cells_names), ]
unknown_cells <- unknown_cells[!(row.names(unknown_cells) %in% endo_cells_names), ]
unknown_cells <- unknown_cells[!(row.names(unknown_cells) %in% fibro_cells_names), ]
unknown_cells[,9:50] <- NULL #Removing from PC9 to PC50


########### train SVM ################
# After creating the traing set we can build the svm

library(e1071)
library(caret)

n <- nrow(known_cells)  # Number of observations
ntrain <- round(n*0.8)  # 80% for training set
set.seed(314)    # Set seed for reproducible results

tindex <- sample(n, ntrain)   # Create a random index
train_known_cells <- known_cells[tindex,]   # Create training set
test_known_cells <- known_cells[-tindex,]   # Create test set

svm1 <- svm(formula = factor(Type)~., data=train_known_cells, cost = 10, gamma = 0.5, type = "C-classification", 
tolerance = 0.01)

summary(svm1)

plot(svm1, train_known_cells, PC_1 ~ PC_3)
plot(svm1, train_known_cells,formula = factor(Type)~.)
plot(svm1, train_known_cells, svSymbol = 1, dataSymbol = 2, symbolPalette = rainbow(4),
     color.palette = terrain.colors)

prediction <- predict(svm1, test_known_cells)
xtab <- table(test_known_cells$Type, prediction)
xtab

# Model Evauation
confusionMatrix(xtab)


x <- predict(svm1, unknown_cells)
prediction_unknown_cells <- as.data.frame(x)
total_labels <- rbind(Card_Fibr_Endo_Immun, prediction_unknown_cells)


########### Visualization SVM using UMPA #####################################

cells_umap <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
cells_umap <- merge(cells_umap, total_labels, by=0, all=F)
rownames(cells_umap) <-cells_umap[,1]
#cells_umap[,1] <- NULL #Removing the first column
colnames(cells_umap)[4] <- "Type"

ggplot(cells_umap) +
  geom_point(aes(x = UMAP_1, y= UMAP_2, color = Type)) 

A <- function(x) substr(x, start = 1, stop = 4)
cells_umap <- data.frame(lapply(cells_umap[1], A), cells_umap[2:4] )
#cells_umap["Row.names"] <- lapply(cells_umap["Row.names"], A)

library(dplyr)

by_vs_am <- cells_umap %>% group_by(Row.names, Type)
by_vs <- by_vs_am %>% summarise(n = n())
#> `summarise()` has grouped output by 'vs'. You can override using the `.groups` argument.
by_vs


write.csv(cells_umap,"cells_annotations_svm.csv", row.names = FALSE)
write.csv(cells_umap,"cells_annotations_svm.xlsx", row.names = TRUE)

saveRDS(cells_umap,file="cells_annotations_svm.Rda")
cells_umap_new <- load("cells_annotations_svm.Rda")


