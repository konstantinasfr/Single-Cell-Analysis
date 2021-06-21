# DEPENDENCIES: integrated
# must execute first seurat_svm.R 

BiocManager::install("writexl")
library(Seurat)
library(ggplot2)
library(cowplot)
library("writexl")

integrated<- readRDS(file = "integrated.rds")

cardio_sham_seurat <- which(integrated@active.ident == 'Cardiomyocytes' & integrated@meta.data[["protocol"]] == "sham" )
cardio_mi_seurat <- which(integrated@active.ident == 'Cardiomyocytes' & integrated@meta.data[["protocol"]] == "mi" )
fibro_seurat <- which(integrated@active.ident == 'Fibroblasts' )
endo_seurat <- which(integrated@active.ident == 'Endothelial' )
immun_seurat <- which(integrated@active.ident == 'Immune' )

integrated_seurat <- SetIdent(object = integrated,cells = cardio_sham_seurat,value = 'cardio_sham')
integrated_seurat <- SetIdent(object = integrated_seurat,cells = cardio_mi_seurat,value = 'cardio_mi')
integrated_seurat <- SetIdent(object = integrated_seurat,cells = fibro_seurat,value = 'fibro')
integrated_seurat <- SetIdent(object = integrated_seurat,cells = endo_seurat,value = 'endo')
integrated_seurat <- SetIdent(object = integrated_seurat,cells = immun_seurat,value = 'immun')

# Visualization
p1 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "protocol")
p2 <- DimPlot(integrated_seurat, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(integrated_seurat, reduction = "umap", split.by = "protocol")

# Markers
markers_cardio_sham_seurat <- FindMarkers(integrated_seurat, ident.1 = "cardio_sham",only.pos = TRUE)
markers_cardio_mi_seurat <- FindMarkers(integrated_seurat, ident.1 = "cardio_mi",only.pos = TRUE)
markers_fibro_seurat <- FindMarkers(integrated_seurat, ident.1 = "fibro",only.pos = TRUE)
markers_endo_seurat <- FindMarkers(integrated_seurat, ident.1 = "endo",only.pos = TRUE)
markers_immun_seurat <- FindMarkers(integrated_seurat, ident.1 = "immun",only.pos = TRUE)

sham_markers_seurat <- FindMarkers(integrated_seurat, ident.1 = "cardio_sham",ident.2 = "cardio_mi",only.pos = TRUE)
mi_markers_seurat <- FindMarkers(integrated_seurat, ident.1 = "cardio_mi",ident.2 = "cardio_sham",only.pos = TRUE )

# Markers Mast
markers_mast_cardio_sham_seurat <- FindMarkers(integrated_seurat, ident.1 = "cardio_sham", test.use = "MAST",only.pos = TRUE)
markers_mast_cardio_mi_seurat <- FindMarkers(integrated_seurat, ident.1 = "cardio_mi", test.use = "MAST",only.pos = TRUE)
markers_mast_fibro_seurat <- FindMarkers(integrated_seurat, ident.1 = "fibro", test.use = "MAST",only.pos = TRUE)
markers_mast_endo_seurat <- FindMarkers(integrated_seurat, ident.1 = "endo", test.use = "MAST",only.pos = TRUE)
markers_mast_immun_seurat <- FindMarkers(integrated_seurat, ident.1 = "immun", test.use = "MAST")

sham_markers_mast_seurat <- FindMarkers(integrated_seurat, ident.1 = "cardio_sham",ident.2 = "cardio_mi", test.use = "MAST",only.pos = TRUE)
mi_markers_mast_seurat <- FindMarkers(integrated_seurat, ident.1 = "cardio_mi",ident.2 = "cardio_sham", test.use = "MAST",only.pos = TRUE)




library("writexl")
# Markers 
m1 <- cbind(" "=rownames(markers_cardio_sham_seurat), markers_cardio_sham_seurat)
#m1 <- m1[order(m1$p_val,decreasing = TRUE),]
write_xlsx(m1,"markers_seurat/wilcoxon/markers_cardio_sham_seurat.xlsx")

m2 <- cbind(" "=rownames(markers_cardio_mi_seurat), markers_cardio_mi_seurat)
#m2 <- m2[order(m2$p_val,decreasing = TRUE),]
write_xlsx(m2,"markers_seurat/wilcoxon/markers_cardio_mi_seurat.xlsx")

m3 <- cbind(" "=rownames(markers_fibro_seurat), markers_fibro_seurat)
#m3 <- m3[order(m3$p_val,decreasing = TRUE),]
write_xlsx(m3,"markers_seurat/wilcoxon/markers_fibro_seurat.xlsx")

m4 <- cbind(" "=rownames(markers_endo_seurat), markers_endo_seurat)
#m4 <- m4[order(m4$p_val,decreasing = TRUE),]
write_xlsx(m4,"markers_seurat/wilcoxon/markers_endo_seurat.xlsx")

m5 <- cbind(" "=rownames(markers_immun_seurat), markers_immun_seurat)
#m5 <- m5[order(m5$p_val,decreasing = TRUE),]
write_xlsx(m5,"markers_seurat/wilcoxon/markers_immun_seurat.xlsx")

m6 <- cbind(" "=rownames(sham_markers_seurat), sham_markers_seurat)
#m6 <- m6[order(m6$p_val,decreasing = TRUE),]
write_xlsx(m6,"markers_seurat/wilcoxon/sham_vs_mi/sham_markers_seurat.xlsx")

m7 <- cbind(" "=rownames(mi_markers_seurat), mi_markers_seurat)
#m7 <- m7[order(m7$p_val,decreasing = TRUE),]
write_xlsx(m7,"markers_seurat/wilcoxon/sham_vs_mi/mi_markers_seurat.xlsx")

# Markers Mast ###############################
m1 <- cbind(" "=rownames(markers_mast_cardio_sham_seurat), markers_mast_cardio_sham_seurat)
#m1 <- m1[order(m1$p_val,decreasing = TRUE),]
write_xlsx(m1,"markers_seurat/mast/markers_mast_cardio_sham_seurat.xlsx")

m2 <- cbind(" "=rownames(markers_mast_cardio_mi_seurat), markers_mast_cardio_mi_seurat)
#m2 <- m2[order(m2$p_val,decreasing = TRUE),]
write_xlsx(m2,"markers_seurat/mast/markers_mast_cardio_mi_seurat.xlsx")

m3 <- cbind(" "=rownames(markers_mast_fibro_seurat), markers_mast_fibro_seurat)
#m3 <- m3[order(m3$p_val,decreasing = TRUE),]
write_xlsx(m3,"markers_seurat/mast/markers_mast_fibro_seurat.xlsx")

m4 <- cbind(" "=rownames(markers_mast_endo_seurat), markers_mast_endo_seurat)
#m4 <- m4[order(m4$p_val,decreasing = TRUE),]
write_xlsx(m4,"markers_seurat/mast/markers_mast_endo_seurat.xlsx")

m5 <- cbind(" "=rownames(markers_mast_immun_seurat), markers_mast_immun_seurat)
#m5 <- m5[order(m5$p_val,decreasing = TRUE),]
write_xlsx(m5,"markers_seurat/mast/markers_mast_immun_seurat.xlsx")

m6 <- cbind(" "=rownames(sham_markers_mast_seurat), sham_markers_mast_seurat)
#m6 <- m6[order(m6$p_val,decreasing = TRUE),]
write_xlsx(m6,"markers_seurat/mast/sham_vs_mi/sham_markers_mast_seurat.xlsx")

m7 <- cbind(" "=rownames(mi_markers_mast_seurat), mi_markers_mast_seurat)
#m7 <- m7[order(m7$p_val,decreasing = TRUE),]
write_xlsx(m7,"markers_seurat/mast/sham_vs_mi/mi_markers_mast_seurat.xlsx")
