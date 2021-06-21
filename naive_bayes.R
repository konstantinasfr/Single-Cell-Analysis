# DEPENDENCIES: known_cells, unknown_cells, integrated 
# must have been created from seurat_svm.R file

# Installing Packages
install.packages("e1071")
install.packages("caTools")
install.packages("caret")

# Loading package
library(e1071)
library(caTools)
library(caret)

n <- nrow(known_cells)  # Number of observations
ntrain <- round(n*0.8)  # 80% for training set
set.seed(314)    # Set seed for reproducible results

tindex <- sample(n, ntrain)   # Create a random index
train_known_cells <- known_cells[tindex,]   # Create training set
test_known_cells <- known_cells[-tindex,]   # Create test set

# Fitting Naive Bayes Model 
# to training dataset
classifier_cl <- naiveBayes(formula = factor(Type)~., data=train_known_cells)
classifier_cl

# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_known_cells)

# Confusion Matrix
cm <- table(test_known_cells$Type, y_pred)
cm

# Model Evauation
confusionMatrix(cm)

row_names_we_lost <- attr(unknown_cells,"row.names")
x <- predict(classifier_cl, unknown_cells)
prediction_unknown_cells_naive <- as.data.frame(x)
row.names(prediction_unknown_cells_naive) <- row_names_we_lost
total_labels_naive <- rbind(Card_Fibr_Endo_Immun, prediction_unknown_cells_naive)

########### Visualization SVM using UMPA #####################################

cells_umap <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
cells_umap <- merge(cells_umap, total_labels_naive, by=0, all=F)
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





