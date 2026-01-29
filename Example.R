# Quick Start Example: Minimal Script to Run All Functions
# This is a simplified script showing the essential workflow

# Load libraries and source functions
library(Seurat)
library(randomForest)

source("Functions.R")

# 1. Data Loading and Preprocessing

#Isoform count matrix

# Import Seurat object
# First, load your Seurat object (e.g., from an RDS file)
Isoform_seurat_obj <- readRDS("Isoform_count_matrix.rds")
# OR create a Seurat object from your data using Seurat's native functions
# Then pass it to import_count_matrix:
Isoform_seurat_obj <- import_count_matrix(Isoform_seurat_obj)

# Pre-process (NormalizeData + FindVariableFeatures)
Isoform_seurat_obj <- preprocess_data(Isoform_seurat_obj, nfeatures = 2000)

# Assign cell types from TSV file
Isoform_seurat_obj <- assign_cell_types(Isoform_seurat_obj,
                                        annotation_file = "Annotation.txt",
                                        barcode_col = "sample",
                                        cell_type_col = "stage")

# Check unique cell types
cell_type_summary <- check_unique_cell_types(Isoform_seurat_obj, plot = TRUE)

# Calculate relative feature usage (RUI)
RIU_result <- calculate_relative_feature_usage(Isoform_seurat_obj, top_n = 100)
RIU= RIU_result$unique_top_features
cat("Unique isoforms across all cell types:", RIU_result$n_unique_features, "\n")

# Scale data and perform dimensional reduction
Isoform_seurat_obj <- scale_and_reduce_dimensions(Isoform_seurat_obj)
# Visualize elbow plot to determine dimensionality:
ElbowPlot(Isoform_seurat_obj)

#Gene count matrix

# Import Seurat object
# First, load your Seurat object (e.g., from an RDS file)
Gene_seurat_obj <- readRDS("Gene_count_matrix.rds")
# OR create a Seurat object from your data using Seurat's native functions
# Then pass it to import_count_matrix:
Gene_seurat_obj <- import_count_matrix(Gene_seurat_obj)

# Pre-process (NormalizeData + FindVariableFeatures)
Gene_seurat_obj <- preprocess_data(Gene_seurat_obj, nfeatures = 2000)

# Assign cell types from TSV file
Gene_seurat_obj <- assign_cell_types(Gene_seurat_obj,
                                     annotation_file = "Annotation.txt",
                                     barcode_col = "sample",
                                     cell_type_col = "stage")

# Check unique cell types
cell_type_summary <- check_unique_cell_types(Gene_seurat_obj, plot = TRUE)

# Calculate relative feature usage (RUI)
RGU_result <- calculate_relative_feature_usage(Gene_seurat_obj, top_n = 100)
RGU= RGU_result$unique_top_features
cat("Unique genes across all cell types:", RGU_result$n_unique_features, "\n")

# Scale data and perform dimensional reduction
Gene_seurat_obj <- scale_and_reduce_dimensions(Gene_seurat_obj)
# Visualize elbow plot to determine dimensionality:
ElbowPlot(Gene_seurat_obj)

# 2. Split cells into training and testing sets 

#Isoform count matrix
split_isoform_matrix <- train_test_split_cells(Isoform_seurat_obj,train_proportion = 0.7, seed = 241)
isoform_train_data <- split_isoform_matrix$train_data
isoform_test_data <- split_isoform_matrix$test_data
train_barcodes <- colnames(isoform_train_data)
test_barcodes  <- colnames(isoform_test_data)

# Extract expression data and labels
isoform_expr_data <- extract_expression_data(isoform_train_data, isoform_test_data , nfeatures = 500)
isoform_train_expr <- isoform_expr_data$train_expr
isoform_test_expr <- isoform_expr_data$test_expr
isoform_train_labels <- isoform_expr_data$train_labels
isoform_test_labels <- isoform_expr_data$test_labels

#Gene count matrix
#Select identical cells for training and testing according to the isoform Seurat object split
gene_train_data <- subset(Gene_seurat_obj, cells = train_barcodes)
gene_test_data  <- subset(Gene_seurat_obj, cells = test_barcodes)

# Extract expression data and labels
gene_expr_data <- extract_expression_data(gene_train_data, gene_test_data , nfeatures = 500)
gene_train_expr <- gene_expr_data$train_expr
gene_test_expr <- gene_expr_data$test_expr
gene_train_labels <- gene_expr_data$train_labels
gene_test_labels <- gene_expr_data$test_labels

# 3. Train global Random Forest model

# Global isoform Rf model
global_isoform_rf_result <- train_random_forest(isoform_train_expr, isoform_train_labels, ntree = 100, seed = 241)
global_isoform_rf_model <- global_isoform_rf_result$model

# Predict cell types
predictions <- predict_cell_types(global_isoform_rf_model, isoform_test_expr, isoform_test_data, 
                                  column_name = "Global_isoform_predicted_cell_type")
isoform_test_data <- predictions$test_data

# Evaluate predictions
evaluation <- evaluate_predictions(test_data = isoform_test_data,
                                   actual_col = "cell_type",
                                   predicted_col = "Global_isoform_predicted_cell_type")
cat("Accuracy global isoform RF:", evaluation$accuracy, "\n")
cat("Accuracy global isoform RF manual calculation:", evaluation$accuracy_manual, "\n")
Global_isform_RF_accuracy <- evaluation$accuracy

# Global gene Rf model
global_gene_rf_result <- train_random_forest(gene_train_expr, gene_train_labels, ntree = 100, seed = 241)
global_gene_rf_model <- global_gene_rf_result$model

# Predict cell types
predictions <- predict_cell_types(global_gene_rf_model, gene_test_expr, gene_test_data, 
                                  column_name = "Global_gene_predicted_cell_type")
gene_test_data <- predictions$test_data

# Evaluate predictions
evaluation <- evaluate_predictions(test_data = gene_test_data,
                                   actual_col = "cell_type",
                                   predicted_col = "Global_gene_predicted_cell_type")
cat("Accuracy global gene RF:", evaluation$accuracy, "\n")
cat("Accuracy global gene RF manual calculation:", evaluation$accuracy_manual, "\n")
Global_gene_RF_accuracy <- evaluation$accuracy

# 4. Compare global level Isoform and Gene Predictions

# Merge isoform and gene predictions
merged_meta <- merge_isoform_gene_predictions(isoform_test_data, gene_test_data,
                                              isoform_pred_col = "Global_isoform_predicted_cell_type",
                                              gene_pred_col = "Global_gene_predicted_cell_type")

# Identify overlapping cells (where both models agree)
overlap_result <- identify_overlapping_cells(isoform_test_data, gene_test_data,
                                             isoform_pred_col = "Global_isoform_predicted_cell_type",
                                             gene_pred_col = "Global_gene_predicted_cell_type",
                                             predicted_global_col= "predicted_global")
global_overlap <- overlap_result$overlap_data
global_non_overlap_isoforms <- overlap_result$non_overlap_isoforms
global_non_overlap_genes <- overlap_result$non_overlap_genes

# Evaluate overlapping predictions
overlap_eval <- evaluate_overlapping_predictions(global_overlap, predicted_col = "predicted_global")
cat("Global level overlapping cells accuracy:", overlap_eval$accuracy, "\n")

# 5. Data Preparation for Subset Model Training

#Isoform count matrix
isoform_subset_train_data <- split_isoform_matrix$train_data
isoform_subset_test_data <- global_non_overlap_isoforms

# Extract expression data and labels
isoform_subset_expr_data <- extract_expression_data(isoform_subset_train_data, isoform_subset_test_data ,features = RIU, nfeatures = 500)
isoform_subset_train_expr <- isoform_subset_expr_data$train_expr
isoform_subset_test_expr <- isoform_subset_expr_data$test_expr
isoform_subset_train_labels <- isoform_subset_expr_data$train_labels
isoform_subset_test_labels <- isoform_subset_expr_data$test_labels

#Gene count matrix
gene_subset_train_data <- gene_train_data
gene_subset_test_data <- global_non_overlap_genes

# Extract expression data and labels
gene_subset_expr_data <- extract_expression_data(gene_subset_train_data, gene_subset_test_data ,features = RGU, nfeatures = 500)
gene_subset_train_expr <- gene_subset_expr_data$train_expr
gene_subset_test_expr <- gene_subset_expr_data$test_expr
gene_subset_train_labels <- gene_subset_expr_data$train_labels
gene_subset_test_labels <- gene_subset_expr_data$test_labels

# 6. Train subset Random Forest model

#Subset isoform Rf model
subset_isoform_rf_result <- train_random_forest(isoform_subset_train_expr, isoform_subset_train_labels, ntree = 100, seed = 241)
subset_isoform_rf_model <- subset_isoform_rf_result$model

# Predict cell types
predictions <- predict_cell_types(subset_isoform_rf_model, isoform_subset_test_expr, isoform_subset_test_data, 
                                  column_name = "Subset_isoform_predicted_cell_type")
isoform_subset_test_data <- predictions$test_data

# Evaluate predictions
evaluation <- evaluate_predictions(test_data = isoform_subset_test_data,
                                   actual_col = "cell_type",
                                   predicted_col = "Subset_isoform_predicted_cell_type")
cat("Accuracy subset isoform RF:", evaluation$accuracy, "\n")
cat("Accuracy subset isoform RF manual calculation:", evaluation$accuracy_manual, "\n")

#Subset gene Rf model
subset_gene_rf_result <- train_random_forest(gene_subset_train_expr, gene_subset_train_labels, ntree = 100, seed = 241)
subset_gene_rf_model <- subset_gene_rf_result$model

# Predict cell types
predictions <- predict_cell_types(subset_gene_rf_model, gene_subset_test_expr, gene_subset_test_data, 
                                  column_name = "Subset_gene_predicted_cell_type")
gene_subset_test_data <- predictions$test_data

# Evaluate predictions
evaluation <- evaluate_predictions(test_data = gene_subset_test_data,
                                   actual_col = "cell_type",
                                   predicted_col = "Subset_gene_predicted_cell_type")
cat("Accuracy subset gene RF:", evaluation$accuracy, "\n")
cat("Accuracy subset gene RF manual calculation:", evaluation$accuracy_manual, "\n")

# 7. Compare Isoform and Gene Predictions in subset level

# Merge isoform and gene predictions
merged_meta_subset <- merge_isoform_gene_predictions(isoform_subset_test_data, gene_subset_test_data,
                                                     isoform_pred_col = "Subset_isoform_predicted_cell_type",
                                                     gene_pred_col = "Subset_gene_predicted_cell_type")

# Identify overlapping cells (where both models agree)
overlap_result_subset <- identify_overlapping_cells(isoform_subset_test_data, gene_subset_test_data,
                                                    isoform_pred_col = "Subset_isoform_predicted_cell_type",
                                                    gene_pred_col = "Subset_gene_predicted_cell_type",
                                                    predicted_global_col = "predicted_subset")
subset_overlap <- overlap_result_subset$overlap_data
subset_non_overlap_isoforms <- overlap_result_subset$non_overlap_isoforms
subset_non_overlap_genes <- overlap_result_subset$non_overlap_genes

# Evaluate overlapping predictions
subset_overlap_eval <- evaluate_overlapping_predictions(subset_overlap,  predicted_col = "predicted_subset",)
cat("Subset overlapping cells accuracy:", subset_overlap_eval$accuracy, "\n")

# 8. Combine final predictions from all sources (global, subset overlap, subset non-overlap)
final_eval <- combined_final_predictions(global_overlap,
                                         subset_overlap,
                                         subset_non_overlap_isoforms,
                                         global_pred_col = "predicted_global",
                                         subset_pred_col = "predicted_subset",
                                         subset_non_overlap_pred_col = "Subset_isoform_predicted_cell_type")
cat("Final Combined Accuracy:", final_eval$accuracy, "\n")
cat("Final Combined Accuracy manually calculated:", final_eval$accuracy_manual, "\n")
Final_Combined_model_accuracy <- final_eval$accuracy_manual

# Save the dataframe with predicted cell types as a CSV
write.csv(
  final_eval$final_predictions_df,
  file = "Predicted cell types.csv",
  row.names = FALSE
)

print(Global_isform_RF_accuracy)
print(Global_gene_RF_accuracy)
print(Final_Combined_model_accuracy)




