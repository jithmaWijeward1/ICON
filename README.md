# IsoHRF

This repository contains an R-based workflow and example datasets for cell type annotation in single-cell RNA-seq data using a hierarchical Random Forest framework that integrates isoform- and gene-based predictions.

## Files
- `scRNA_analysis_functions.R` - Contains all the analysis functions
- `example.R` - Minimal script for quick start
- `README.md` - This file

## Installation
First, install required R packages:

```r
install.packages("Seurat")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("Matrix")
install.packages("randomForest")
```

## Functions

### 1. `import_count_matrix()`
Validates and imports a Seurat object.

**Note:** This function only accepts Seurat objects. To create a Seurat object from files, use Seurat's native functions first (e.g., `readRDS()`, `CreateSeuratObject()`, `Read10X()`, etc.).

**Example:**
```r
# Load from RDS file
seurat_obj <- readRDS("path/to/seurat_object.rds")
seurat_obj <- import_count_matrix(seurat_obj)

# Or create from count matrix
count_matrix <- read.csv("data.csv", row.names = 1)
seurat_obj <- CreateSeuratObject(counts = count_matrix)
seurat_obj <- import_count_matrix(seurat_obj)
```

### 2. `preprocess_data()`
Performs data normalization and finds variable features.

**Includes:**
- `NormalizeData()` - Normalizes the count data
- `FindVariableFeatures()` - Identifies highly variable genes

**Example:**
```r
seurat_obj <- preprocess_data(seurat_obj, nfeatures = 2000)
```

### 3. `assign_cell_types()`
Assigns cell types by loading annotations from a TSV/CSV/TXT file and merging with Seurat object metadata.

**Parameters:**
- `annotation_file` - Path to TSV file containing cell type annotations (optional if `annotation_df` is provided)
- `annotation_df` - Data frame containing cell type annotations (optional if `annotation_file` is provided)
- `barcode_col` - Name of the column containing cell barcodes in annotation file (default: "cell")
- `cell_type_col` - Name of the column containing cell types in annotation file (default: "donor_id")
- `column_name` - Name of the metadata column to store cell types (default: "cell_type")
- `sep` - Separator for TSV file (default: "\t")

**Example:**
```r
# Load from TSV file
seurat_obj <- assign_cell_types(seurat_obj,
                                annotation_file = "path/to/donor_ids.tsv",
                                barcode_col = "cell",
                                cell_type_col = "donor_id")

# Or use a pre-loaded data frame
donor_ids <- read.csv("path/to/donor_ids.tsv", header = TRUE, sep = "\t")
seurat_obj <- assign_cell_types(seurat_obj,
                                annotation_df = donor_ids,
                                barcode_col = "cell",
                                cell_type_col = "donor_id")
```

### 4. `check_unique_cell_types()`
Checks and summarizes unique cell types in the dataset.

**Example:**
```r
summary <- check_unique_cell_types(seurat_obj, plot = TRUE)
```

### 5. `calculate_relative_feature_usage()`
Calculates top highly expressed features (isoforms/genes) per cell type and returns a unique set of features across all cell types. This is useful for calculating Relative Feature Usage Index (Isoforms: RIU and Genes: RGU).

**Parameters:**
- `column_name` - Name of the metadata column containing cell types (default: "cell_type")
- `top_n` - Number of top features to select per cell type (default: 100)
- `verbose` - Print progress messages (default: TRUE)

**Returns:**
- `top_features_list` - Named list of top features per cell type
- `unique_top_features` - Unique set of features across all cell types
- `n_unique_features` - Number of unique features

**Example:**
```r
rui_result <- calculate_relative_feature_usage(seurat_obj, top_n = 100)
top_features <- rui_result$unique_top_features
cat("Unique features:", rui_result$n_unique_features, "\n")
```

### 6. `scale_and_reduce_dimensions()`
Scales all features and performs PCA for dimensional reduction.

**Parameters:**
- `features` - Features to scale. If NULL, scales all features (default: NULL)
- `use_variable_features` - Whether to use variable features for PCA (default: TRUE)
- `npcs` - Number of principal components to compute (default: 50)
- `verbose` - Print progress messages (default: TRUE)

**Example:**
```r
seurat_obj <- scale_and_reduce_dimensions(seurat_obj)
# Visualize elbow plot to determine dimensionality:
ElbowPlot(seurat_obj)
```

### 7. `train_test_split_cells()`
Splits cells into training (70%) and testing (30%) sets, stratified by cell type. The split is performed separately for each cell type to ensure balanced representation.

**Parameters:**
- `column_name` - Name of the metadata column containing cell types (default: "cell_type")
- `train_proportion` - Proportion of cells to use for training (default: 0.7)
- `seed` - Random seed for reproducibility (default: 241)
- `verbose` - Print progress messages (default: TRUE)

**Returns:**
- `train_data` - Seurat object with training cells
- `test_data` - Seurat object with testing cells
- `train_cells` - Vector of training cell names
- `test_cells` - Vector of testing cell names
- `distribution` - List with train and test cell type distributions

**Example:**
```r
split_result <- train_test_split_cells(seurat_obj, seed = 241)
train_data <- split_result$train_data
test_data <- split_result$test_data
# View distributions
print(split_result$distribution$train_df)
print(split_result$distribution$test_df)
```

### 8. `extract_expression_data()`
Extracts normalized expression matrices and cell type labels for training and testing sets. Optionally finds variable features in the training set.

**Parameters:**
- `train_data` - Seurat object for training set
- `test_data` - Seurat object for testing set
- `column_name` - Name of the metadata column containing cell types (default: "cell_type")
- `find_variable_features` - Logical, whether to find variable features in training set (default: TRUE)
- `nfeatures` - Number of variable features to select (default: 500)
- `selection_method` - Method for selecting variable features (default: "vst")
- `use_variable_features` - Logical, whether to extract only variable features (default: TRUE)
- `layer` - Layer to extract from (default: "data" for normalized data)
- `verbose` - Print progress messages (default: TRUE)

**Returns:**
- `train_expr` - Expression matrix for training set (data frame)
- `test_expr` - Expression matrix for testing set (data frame)
- `train_labels` - Cell type labels for training set (factor)
- `test_labels` - Cell type labels for testing set (factor)
- `variable_features` - Vector of variable feature names (if found)

**Example:**
```r
expr_data <- extract_expression_data(train_data, test_data, nfeatures = 500)
train_expr <- expr_data$train_expr
test_expr <- expr_data$test_expr
train_labels <- expr_data$train_labels
test_labels <- expr_data$test_labels
```

### 9. `train_random_forest()`
Trains a Random Forest model for cell type classification using expression data.

**Parameters:**
- `train_expr` - Expression matrix for training set (features as rows, cells as columns)
- `train_labels` - Cell type labels for training set (factor or character vector)
- `ntree` - Number of trees to grow (default: 100)
- `importance` - Logical, whether to calculate variable importance (default: TRUE)
- `seed` - Random seed for reproducibility (default: 241)
- `verbose` - Print progress messages (default: TRUE)
- `...` - Additional arguments passed to randomForest()

**Returns:**
- `model` - Trained Random Forest model
- `training_time` - Training time in seconds

**Example:**
```r
rf_result <- train_random_forest(train_expr, train_labels, ntree = 100)
rf_model <- rf_result$model
```

### 10. `predict_cell_types()`
Predicts cell types for test data using a trained Random Forest model. Optionally stores predictions in Seurat object metadata.

**Parameters:**
- `rf_model` - Trained Random Forest model (from train_random_forest)
- `test_expr` - Expression matrix for test set (features as rows, cells as columns)
- `test_data` - Seurat object for test set (optional, for storing predictions)
- `column_name` - Name of the metadata column to store predictions (default: "predicted_cell_type")
- `return_probs` - Logical, whether to return class probabilities (default: FALSE)
- `verbose` - Print progress messages (default: TRUE)

**Returns:**
- `predicted_labels` - Predicted cell type labels (factor)
- `predicted_probs` - Class probabilities matrix (if return_probs = TRUE)
- `test_data` - Seurat object with predictions in metadata (if test_data provided)

**Example:**
```r
predictions <- predict_cell_types(rf_model, test_expr, test_data)
predicted_labels <- predictions$predicted_labels
test_data <- predictions$test_data
```

### 11. `evaluate_predictions()`
Evaluates predicted cell types by comparing with actual cell types. Calculates confusion matrix and accuracy metrics.

**Parameters:**
- `actual_labels` - Actual cell type labels (optional if test_data provided)
- `predicted_labels` - Predicted cell type labels (optional if test_data provided)
- `test_data` - Seurat object with predictions (optional, for extracting labels)
- `actual_col` - Name of metadata column containing actual labels (default: "cell_type")
- `predicted_col` - Name of metadata column containing predicted labels (default: "predicted_cell_type")
- `verbose` - Print evaluation results (default: TRUE)

**Returns:**
- `confusion_matrix` - Confusion matrix table
- `accuracy` - Overall accuracy (from confusion matrix)
- `accuracy_manual` - Manual accuracy calculation (matching cells / total cells)
- `accurate_cells_count` - Number of correctly predicted cells
- `total_cells` - Total number of cells

**Example:**
```r
# Using Seurat object
evaluation <- evaluate_predictions(test_data = test_data,
                                   actual_col = "cell_type",
                                   predicted_col = "predicted_cell_type")
print(evaluation$confusion_matrix)
cat("Accuracy:", evaluation$accuracy, "\n")

# Using extracted labels
evaluation <- evaluate_predictions(actual_labels = test_labels,
                                   predicted_labels = predicted_labels)
```

### 12. `merge_isoform_gene_predictions()`
Merges predictions from isoform-based and gene-based Seurat objects into a single data frame.

**Parameters:**
- `test_data_isoforms` - Seurat object with isoform-based predictions
- `test_data_genes` - Seurat object with gene-based predictions
- `isoform_pred_col` - Name of the column containing isoform predictions (default: "predicted_label_global_isoforms")
- `gene_pred_col` - Name of the column containing gene predictions (default: "predicted_label_global_genes")
- `true_labels_col` - Name of the column containing true labels (default: "cell_type")
- `use_isoform_labels` - Logical, whether to use isoform object for true labels (default: TRUE)
- `verbose` - Print progress messages (default: TRUE)

**Returns:**
- A data frame with merged metadata containing true labels and both predictions

**Example:**
```r
merged_meta <- merge_isoform_gene_predictions(test_data, test_data_genes)
```

### 13. `identify_overlapping_cells()`
Identifies cells where isoform and gene predictions agree (overlapping) or disagree (non-overlapping). Returns subset Seurat objects for each category.

**Parameters:**
- `test_data_isoforms` - Seurat object with isoform-based predictions
- `test_data_genes` - Seurat object with gene-based predictions
- `merged_meta` - Merged metadata data frame (from merge_isoform_gene_predictions). If NULL, will be created automatically.
- `isoform_pred_col` - Name of the column containing isoform predictions (default: "predicted_label_global_isoforms")
- `gene_pred_col` - Name of the column containing gene predictions (default: "predicted_label_global_genes")
- `predicted_global_col` - Name of the column to store agreed predictions in overlapping cells (default: "predicted_global")
- `verbose` - Print progress messages (default: TRUE)

**Returns:**
- `overlap_data` - Seurat object with overlapping cells (same prediction from both)
- `non_overlap_isoforms` - Seurat object with non-overlapping cells from isoform object
- `non_overlap_genes` - Seurat object with non-overlapping cells from gene object
- `overlap_barcodes` - Vector of overlapping cell barcodes
- `non_overlap_barcodes` - Vector of non-overlapping cell barcodes
- `merged_meta` - Merged metadata data frame

**Example:**
```r
overlap_result <- identify_overlapping_cells(test_data, test_data_genes)
global_overlap <- overlap_result$overlap_data
global_non_overlap_isoforms <- overlap_result$non_overlap_isoforms
global_non_overlap_genes <- overlap_result$non_overlap_genes
```

### 14. `evaluate_overlapping_predictions()`
Evaluates prediction accuracy for cells where isoform and gene predictions agree (overlapping cells).

**Parameters:**
- `overlap_data` - Seurat object with overlapping cells (from identify_overlapping_cells)
- `actual_col` - Name of metadata column containing actual labels (default: "cell_type")
- `predicted_col` - Name of metadata column containing predicted labels (default: "predicted_global")
- `verbose` - Print evaluation results (default: TRUE)

**Returns:**
- `confusion_matrix` - Confusion matrix table
- `accuracy` - Overall accuracy (from confusion matrix)
- `accuracy_manual` - Manual accuracy calculation (matching cells / total cells)
- `accurate_cells_count` - Number of correctly predicted cells
- `total_cells` - Total number of cells

**Example:**
```r
overlap_eval <- evaluate_overlapping_predictions(global_overlap)
print(overlap_eval$confusion_matrix)
cat("Accuracy:", overlap_eval$accuracy, "\n")
```

### 15. `combined_final_predictions()`
Combines predictions from global overlap, subset overlap, and subset non-overlap Seurat objects to create a final combined evaluation. Calculates overall accuracy using both manual calculation and confusion matrix.

**Parameters:**
- `global_overlap` - Seurat object with global overlapping cells (from identify_overlapping_cells)
- `subset_overlap` - Seurat object with subset overlapping cells
- `subset_non_overlap_isoforms` - Seurat object with subset non-overlapping cells (isoform-based)
- `global_pred_col` - Name of the column containing global predictions (default: "predicted_global")
- `subset_pred_col` - Name of the column containing subset overlap predictions (default: "predicted_subset")
- `subset_non_overlap_pred_col` - Name of the column containing subset non-overlap predictions (default: "predicted_label_subset_isoforms")
- `actual_col` - Name of the metadata column containing actual labels (default: "cell_type")
- `verbose` - Print progress messages and results (default: TRUE)

**Returns:**
- `final_predictions_df` - Combined data frame with all predictions 
- `confusion_matrix` - Confusion matrix table
- `accuracy` - Overall accuracy (from confusion matrix)
- `accuracy_manual` - Manual accuracy calculation (matching cells / total cells)
- `accurate_cells_count` - Number of correctly predicted cells
- `total_cells` - Total number of cells

**Example:**
```r
final_eval <- combined_final_predictions(global_overlap,
                                        subset_overlap,
                                        subset_non_overlap_isoforms)
cat("Final Combined Accuracy:", final_eval$accuracy, "\n")
print(final_eval$confusion_matrix)

# Save the dataframe with predicted cell types as a CSV
write.csv(
  final_eval$final_predictions_df,
  file = "Predicted cell types.csv",
  row.names = FALSE
)
```

## Quick Start

Usage note: If you have the isoform- and gene-level Seurat objects and the required cell type annotation information (for use with assign_cell_types(), including annotation_file, barcode_col, and cell_type_col), you can use the Example.R script to quickly run the complete analysis workflow using the functions provided in this repository.

## Notes

- Make sure your count matrix has genes as rows and cells as columns
- Cell barcodes should match between your count matrix and cell type assignments
- The annotation file should contain columns for cell barcodes and cell types
- The function automatically handles merging and preserves row names in the Seurat object metadata
