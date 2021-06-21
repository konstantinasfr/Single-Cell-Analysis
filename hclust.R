# DEPENDENCIES: integrated 
# must execute seurat_svm.R file first


expr_pcas = as.data.frame(integrated@reductions[["pca"]]@cell.embeddings)
expr_pcas[,9:50] <- NULL #Removing the first column


########## hierrarhical ###########

# Dissimilarity matrix
d <- dist(expr_pcas, method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc <- hclust(d, method = "ward.D" )
# Plot the obtained dendrogram
#plot(hc1, cex = 0.6, hang = -1)
x<-cutree(hc, k = 4)
cluster_labels <- as.data.frame(x)

cluster_labels[cluster_labels["x"]== 2] <- "Immune"
cluster_labels[cluster_labels["x"]== 1] <- "Cardiomyocytes"
cluster_labels[cluster_labels["x"]== 3] <- "Fibroblasts"
cluster_labels[cluster_labels["x"]== 4] <- "Endothelial"

cells_umap <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
cells_umap <- merge(cells_umap, cluster_labels, by=0, all=F)
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




