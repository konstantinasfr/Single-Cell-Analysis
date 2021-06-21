# This code is used only for the creation of the plots that compare the clusterings
# it must be executed after the creation of the clustres so after seurat_svm.R hclus.R and naive_bayes.R


library(sqldf)
install.packages('sqldf')

integrated<- readRDS(file = "integrated.rds")

y=c()
for (i in 1:1735) { 
  y[i]=NA
}


total_cells <- as.data.frame(y)
row.names(total_cells) <- attr(total_labels, "row.names")

# gia seurat
x<-integrated@active.ident
seurat_cells = as.data.frame(x)


#cluster_labels

#total_labels
cardio1 <- filter(seurat_cells, x == "Cardiomyocytes")
cardio1["Cells"] <- attr(cardio1, "row.names")

cardio_sham1<-cardio1[ cardio1$Cells %like% "sham",]
cardio_sham1[,2] <- NULL #Removing the first column
cardio_sham_seurat1 <- attr(cardio_sham1, "row.names")

cardio_mi1<-cardio1[rownames(cardio1) %like% "mi",]
cardio_mi1[,2] <- NULL #Removing the first column
cardio_mi_seurat1 <- attr(cardio_mi1, "row.names")
#cardio1[,2] <- NULL #Removing the first column

fibro1 <- filter(seurat_cells, x == "Fibroblasts" )
fibro1["Cells"] <- attr(fibro1, "row.names")
fibro_seurat1 <- attr(fibro1, "row.names")

endo1 <- filter(seurat_cells, x == "Endothelial")
endo1["Cells"] <- attr(endo1, "row.names")
endo_seurat1 <- attr(endo1, "row.names")

immun1 <- filter(seurat_cells, x == "Immune" )
immun1["Cells"] <- attr(immun1, "row.names")
immun_seurat1 <- attr(immun1, "row.names")




#total_labels_naive

cardio2 <- filter(cluster_labels, x == "Cardiomyocytes")
cardio2["Cells"] <- attr(cardio2, "row.names")

cardio_sham2<-cardio2[ cardio2$Cells %like% "sham",]
cardio_sham2[,2] <- NULL #Removing the first column
cardio_sham_seurat2 <- attr(cardio_sham2, "row.names")

cardio_mi2<-cardio2[rownames(cardio2) %like% "mi",]
cardio_mi2[,2] <- NULL #Removing the first column
cardio_mi_seurat2<- attr(cardio_mi2, "row.names")
#cardio2[,2] <- NULL #Removing the first column

fibro2 <- filter(cluster_labels, x == "Fibroblasts" )
fibro2["Cells"] <- attr(fibro2, "row.names")
fibro_seurat2 <- attr(fibro2, "row.names")

endo2 <- filter(cluster_labels, x == "Endothelial")
endo2["Cells"] <- attr(endo2, "row.names")
endo_seurat2 <- attr(endo2, "row.names")

immun2 <- filter(cluster_labels, x == "Immune" )
immun2["Cells"] <- attr(immun2, "row.names")
immun_seurat2 <- attr(immun2, "row.names")

######################### COMPARE ########################


########### Cardio ##########################
cardio_common <- merge(cardio1, cardio2, by=0, all=F)
cardio_common["x"] <- "common"
rownames(cardio_common) <- cardio_common[,3]
cardio_common[,1:5] <- NULL #Removing from PC9 to PC50

cardio_not_common1 <- anti_join(cardio2, cardio1)
cardio_not_common1["x"] <- "hierarchical"
cardio_not_common1[,2] <- NULL #Removing from PC9 to PC50
cardio_not_common2 <- anti_join(cardio1, cardio2)
cardio_not_common2["x"] <- "seurat"
cardio_not_common2[,2] <- NULL #Removing from PC9 to PC50

cardio_temp <- rbind(cardio_not_common2, cardio_not_common1)
cardio_total <- rbind(cardio_common, cardio_temp)

total_common <- merge(total_cells, cardio_total, by=0, all=T)
rownames(total_common) <- total_common[,1]
total_common[,1:2] <- NULL #Removing from PC9 to PC50

cells_umap <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
cells_umap <- merge(cells_umap, total_common, by=0, all=F)
rownames(cells_umap) <-cells_umap[,1]
#cells_umap[,1] <- NULL #Removing the first column
colnames(cells_umap)[4] <- "Type"

ggplot(cells_umap) +
  geom_point(aes(x = UMAP_1, y= UMAP_2, color = Type)) + ggtitle("Cardiomyocytes")

########### Fibroblasts ##########################
cardio_common <- merge(fibro1, fibro2, by=0, all=F)
cardio_common["x"] <- "common"
rownames(cardio_common) <- cardio_common[,3]
cardio_common[,1:5] <- NULL #Removing from PC9 to PC50

cardio_not_common1 <- anti_join(fibro2, fibro1)
cardio_not_common1["x"] <- "hierarchical"
cardio_not_common1[,2] <- NULL #Removing from PC9 to PC50
cardio_not_common2 <- anti_join(fibro1, fibro2)
cardio_not_common2["x"] <- "seurat"
cardio_not_common2[,2] <- NULL #Removing from PC9 to PC50

cardio_temp <- rbind(cardio_not_common2, cardio_not_common1)
cardio_total <- rbind(cardio_common, cardio_temp)

total_common <- merge(total_cells, cardio_total, by=0, all=T)
rownames(total_common) <- total_common[,1]
total_common[,1:2] <- NULL #Removing from PC9 to PC50

cells_umap <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
cells_umap <- merge(cells_umap, total_common, by=0, all=F)
rownames(cells_umap) <-cells_umap[,1]
#cells_umap[,1] <- NULL #Removing the first column
colnames(cells_umap)[4] <- "Type"

ggplot(cells_umap) +
  geom_point(aes(x = UMAP_1, y= UMAP_2, color = Type)) + ggtitle("Fibroblasts")

########### Endothelial ##########################
cardio_common <- merge(endo1, endo2, by=0, all=F)
cardio_common["x"] <- "common"
rownames(cardio_common) <- cardio_common[,3]
cardio_common[,1:5] <- NULL #Removing from PC9 to PC50

cardio_not_common1 <- anti_join(endo2, endo1)
cardio_not_common1["x"] <- "svm"
cardio_not_common1[,2] <- NULL #Removing from PC9 to PC50
cardio_not_common2 <- anti_join(endo1, endo2)
cardio_not_common2["x"] <- "seurat"
cardio_not_common2[,2] <- NULL #Removing from PC9 to PC50

cardio_temp <- rbind(cardio_not_common2, cardio_not_common1)
cardio_total <- rbind(cardio_common, cardio_temp)

total_common <- merge(total_cells, cardio_total, by=0, all=T)
rownames(total_common) <- total_common[,1]
total_common[,1:2] <- NULL #Removing from PC9 to PC50

cells_umap <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
cells_umap <- merge(cells_umap, total_common, by=0, all=F)
rownames(cells_umap) <-cells_umap[,1]
#cells_umap[,1] <- NULL #Removing the first column
colnames(cells_umap)[4] <- "Type"

ggplot(cells_umap) +
  geom_point(aes(x = UMAP_1, y= UMAP_2, color = Type)) + ggtitle("Endothelial")

########### Immune ##########################
cardio_common <- merge(immun1, immun2, by=0, all=F)
cardio_common["x"] <- "common"
rownames(cardio_common) <- cardio_common[,3]
cardio_common[,1:5] <- NULL #Removing from PC9 to PC50

cardio_not_common1 <- anti_join(immun2, immun1)
cardio_not_common1["x"] <- "svm"
cardio_not_common1[,2] <- NULL #Removing from PC9 to PC50
cardio_not_common2 <- anti_join(immun1, immun2)
cardio_not_common2["x"] <- "seurat"
cardio_not_common2[,2] <- NULL #Removing from PC9 to PC50

cardio_temp <- rbind(cardio_not_common2, cardio_not_common1)
cardio_total <- rbind(cardio_common, cardio_temp)

total_common <- merge(total_cells, cardio_total, by=0, all=T)
rownames(total_common) <- total_common[,1]
total_common[,1:2] <- NULL #Removing from PC9 to PC50

cells_umap <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
cells_umap <- merge(cells_umap, total_common, by=0, all=F)
rownames(cells_umap) <-cells_umap[,1]
#cells_umap[,1] <- NULL #Removing the first column
colnames(cells_umap)[4] <- "Type"

ggplot(cells_umap) +
  geom_point(aes(x = UMAP_1, y= UMAP_2, color = Type)) + ggtitle("Immune")

