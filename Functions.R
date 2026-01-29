# Single-Cell RNA-seq Analysis Functions
# This script contains functions for importing, preprocessing, and analyzing single-cell data

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(randomForest)

#' Import Seurat Object
#' 
#' @param file_path A Seurat object
#' @return A Seurat object
import_count_matrix <- function(file_path) {
  
  # Check if file_path is already a Seurat object
  if (inherits(file_path, "Seurat")) {
    cat("Input is already a Seurat object. Returning as-is.\n")
    cat("  Cells:", ncol(file_path), "\n")
    cat("  Features:", nrow(file_path), "\n")
    return(file_path)
  }
  
  stop("file_path must be a Seurat object")
}


#' Pre-process Seurat Object
#' 
#' Performs normalization and finds variable features
#' 
#' @param seurat_obj A Seurat object
#' @param normalization_method Method for normalization: "LogNormalize" (default), "CLR", or "RC"
#' @param scale_factor Scale factor for normalization (default: 10000)
#' @param nfeatures Number of variable features to select (default: 2000)
#' @param selection_method Method for selecting variable features: "vst" (default), "mean.var.plot", or "dispersion"
#' @param verbose Print progress messages (default: TRUE)
#' @return A preprocessed Seurat object
preprocess_data <- function(seurat_obj, 
                            normalization_method = "LogNormalize",
                            scale_factor = 10000,
                            nfeatures = 2000,
                            selection_method = "vst",
                            verbose = TRUE) {
  
  if (verbose) cat("Starting preprocessing...\n")
  
  # Step 1: Normalize data
  if (verbose) cat("Normalizing data...\n")
  seurat_obj <- NormalizeData(seurat_obj,
                              normalization.method = normalization_method,
                              scale.factor = scale_factor,
                              verbose = verbose)
  
  # Step 2: Find variable features
  if (verbose) cat("Finding variable features...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                     selection.method = selection_method,
                                     nfeatures = nfeatures,
                                     verbose = verbose)
  
  if (verbose) {
    cat("Preprocessing complete!\n")
    cat("  Variable features found:", length(VariableFeatures(seurat_obj)), "\n")
  }
  
  return(seurat_obj)
}


#' Assign Cell Types to Data
#' 
#' Assigns cell types by loading annotations from a TSV file and merging with Seurat object metadata
#' 
#' @param seurat_obj A Seurat object
#' @param annotation_file Path to TSV file containing cell type annotations (optional if annotation_df is provided)
#' @param annotation_df Data frame containing cell type annotations (optional if annotation_file is provided)
#' @param barcode_col Name of the column containing cell barcodes in annotation file (default: "cell")
#' @param cell_type_col Name of the column containing cell types in annotation file (default: "donor_id")
#' @param column_name Name of the metadata column to store cell types (default: "cell_type")
#' @param sep Separator for TSV file (default: "\t")
#' @param verbose Print progress messages (default: TRUE)
#' @return A Seurat object with cell type assignments in metadata
#' @examples
#' # Using file path
#' seurat_obj <- assign_cell_types(seurat_obj, 
#'                                 annotation_file = "path/to/donor_ids.tsv")
#' 
#' # Using data frame
#' annotations <- read.csv("donor_ids.tsv", sep = "\t", header = TRUE)
#' seurat_obj <- assign_cell_types(seurat_obj, annotation_df = annotations)
assign_cell_types <- function(seurat_obj,
                              annotation_file = NULL,
                              annotation_df = NULL,
                              barcode_col = "cell",
                              cell_type_col = "donor_id",
                              column_name = "cell_type",
                              #sep = "\t",
                              verbose = TRUE) {
  
  if (verbose) cat("Assigning cell types...\n")
  
  # Load cell type annotations
  if (!is.null(annotation_file)) {
    if (verbose) cat("Loading annotations from file:", annotation_file, "\n")
    
    # Detect file extension
    file_ext <- tolower(tools::file_ext(annotation_file))
    
    # Set separator based on file type
    sep <- switch(
      file_ext,
      "tsv" = "\t",
      "txt" = "\t",
      "csv" = ",",
      stop("Unsupported file type: .", file_ext,
           ". Supported types are: .tsv, .csv, .txt")
    )
    
    donor_ids <- read.csv(annotation_file, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  } else if (!is.null(annotation_df)) {
    if (verbose) cat("Using provided annotation data frame\n")
    donor_ids <- annotation_df
  } else {
    stop("Either 'annotation_file' or 'annotation_df' must be provided")
  }
  
  # Check if required columns exist
  if (!barcode_col %in% colnames(donor_ids)) {
    stop(paste("Column", barcode_col, "not found in annotation file. Available columns:",
               paste(colnames(donor_ids), collapse = ", ")))
  }
  if (!cell_type_col %in% colnames(donor_ids)) {
    stop(paste("Column", cell_type_col, "not found in annotation file. Available columns:",
               paste(colnames(donor_ids), collapse = ", ")))
  }
  
  # Subset the two columns
  if (verbose) cat("Subsetting columns:", barcode_col, "and", cell_type_col, "\n")
  cell_types <- donor_ids[, c(barcode_col, cell_type_col)]
  
  # Rename columns to barcode and cell_type
  colnames(cell_types) <- c("barcode", "cell_type")
  
  # Step 1: Convert row names (barcodes) into a regular column
  if (verbose) cat("Adding barcode column to metadata...\n")
  seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
  
  # Step 2: Perform the left_join to add the 'cell_type' column
  if (verbose) cat("Merging cell type annotations...\n")
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    left_join(cell_types, by = "barcode")
  
  # Step 3: Reassign the 'barcode' column back to row names
  rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$barcode
  
  # Step 4: Remove the barcode column after merging
  seurat_obj@meta.data$barcode <- NULL
  
  # Rename the merged column if different from default
  if (column_name != "cell_type" && "cell_type" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data[[column_name]] <- seurat_obj@meta.data$cell_type
    seurat_obj@meta.data$cell_type <- NULL
  }
  
  if (verbose) {
    cat("Cell type assignment complete!\n")
    cat("  Cell types assigned to column:", column_name, "\n")
    # Count how many cells got annotations
    n_annotated <- sum(!is.na(seurat_obj@meta.data[[column_name]]))
    cat("  Cells with annotations:", n_annotated, "out of", ncol(seurat_obj), "\n")
  }
  
  return(seurat_obj)
}


#' Check Unique Cell Types in Dataset
#' 
#' Returns information about unique cell types in the dataset
#' 
#' @param seurat_obj A Seurat object with cell type assignments
#' @param column_name Name of the metadata column containing cell types (default: "cell_type")
#' @param return_table Logical, whether to return a table with counts (default: TRUE)
#' @param plot Logical, whether to create a bar plot (default: FALSE)
#' @return A data frame with cell type counts (if return_table = TRUE) or a plot (if plot = TRUE)
check_unique_cell_types <- function(seurat_obj,
                                    column_name = "cell_type",
                                    return_table = TRUE,
                                    plot = FALSE) {
  
  # Check if column exists
  if (!column_name %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", column_name, "not found in metadata. Available columns:",
               paste(colnames(seurat_obj@meta.data), collapse = ", ")))
  }
  
  # Get cell types
  cell_types <- seurat_obj@meta.data[[column_name]]
  
  # Remove NA values
  cell_types <- cell_types[!is.na(cell_types)]
  
  # Get unique cell types and counts
  unique_types <- unique(cell_types)
  type_counts <- table(cell_types)
  
  # Create summary data frame
  summary_df <- data.frame(
    Cell_Type = names(type_counts),
    Count = as.numeric(type_counts),
    Percentage = round(as.numeric(type_counts) / length(cell_types) * 100, 2)
  )
  summary_df <- summary_df[order(-summary_df$Count), ]
  
  # Print summary
  cat("Cell Type Summary:\n")
  cat("==================\n")
  cat("Total unique cell types:", length(unique_types), "\n")
  cat("Total cells:", length(cell_types), "\n\n")
  print(summary_df)
  
  # Create plot if requested
  if (plot) {
    p <- ggplot(summary_df, aes(x = reorder(Cell_Type, Count), y = Count, fill = Cell_Type)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = "Cell Type Distribution",
           x = "Cell Type",
           y = "Number of Cells") +
      theme_minimal() +
      theme(legend.position = "none")
    print(p)
  }
  
  if (return_table) {
    return(summary_df)
  } else {
    return(invisible(NULL))
  }
}


#' Calculate Relative Feature Usage (RUI)
#' 
#' Calculates top highly expressed features (isoforms/genes) per cell type
#' and returns a unique set of features across all cell types
#' 
#' @param seurat_obj A Seurat object with cell type assignments
#' @param column_name Name of the metadata column containing cell types (default: "cell_type")
#' @param top_n Number of top features to select per cell type (default: 100)
#' @param verbose Print progress messages (default: TRUE)
#' @return A list containing:
#'   - `top_features_list`: Named list of top features per cell type
#'   - `unique_top_features`: Unique set of features across all cell types
#'   - `n_unique_features`: Number of unique features
#' @examples
#' result <- calculate_relative_feature_usage(seurat_obj, top_n = 100)
#' top_features <- result$unique_top_features
calculate_relative_feature_usage <- function(seurat_obj,
                                             column_name = "cell_type",
                                             top_n = 100,
                                             verbose = TRUE) {
  
  if (verbose) cat("Calculating relative feature usage...\n")
  
  # Check if column exists
  if (!column_name %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", column_name, "not found in metadata. Available columns:",
               paste(colnames(seurat_obj@meta.data), collapse = ", ")))
  }
  
  # Get unique cell types
  cell_types <- unique(seurat_obj@meta.data[[column_name]])
  cell_types <- cell_types[!is.na(cell_types)]
  
  if (verbose) cat("  Found", length(cell_types), "cell types\n")
  
  # Create a named list of Seurat objects, each corresponding to a single cell type
  if (verbose) cat("  Creating cell type subsets...\n")
  isoform_subsets <- lapply(cell_types, function(ct) {
    # Get cells matching this cell type
    cells_to_keep <- seurat_obj@meta.data[[column_name]] == ct
    cells_to_keep[is.na(cells_to_keep)] <- FALSE
    cell_names <- colnames(seurat_obj)[cells_to_keep]
    subset(seurat_obj, cells = cell_names)
  })
  names(isoform_subsets) <- cell_types
  
  # Create an empty list to store top features per cell type
  top_isoform_list <- list()
  
  # Loop through each cell type subset
  if (verbose) cat("  Calculating top features per cell type...\n")
  for (ct in names(isoform_subsets)) {
    # Extract Seurat object
    seurat_subset <- isoform_subsets[[ct]]
    
    # Get normalized expression matrix
    norm_expr <- GetAssayData(seurat_subset, layer = "data")
    
    # Sum expression per feature across all cells
    isoform_totals <- Matrix::rowSums(norm_expr)
    
    # Get top N highly expressed features
    top_isoforms <- names(sort(isoform_totals, decreasing = TRUE))[1:min(top_n, length(isoform_totals))]
    
    # Store in the list
    top_isoform_list[[ct]] <- top_isoforms
    
    if (verbose) cat("    ", ct, ":", length(top_isoforms), "top features\n")
  }
  
  # Combine all top features across cell types into one vector
  combined_top_isoforms <- unlist(top_isoform_list)
  
  # Get the unique set (i.e., remove duplicates)
  unique_top_isoforms <- unique(combined_top_isoforms)
  
  if (verbose) {
    cat("  Unique features across all cell types:", length(unique_top_isoforms), "\n")
    cat("Relative feature usage calculation complete!\n")
  }
  
  return(list(
    top_features_list = top_isoform_list,
    unique_top_features = unique_top_isoforms,
    n_unique_features = length(unique_top_isoforms)
  ))
}


#' Scale Data and Perform Dimensional Reduction
#' 
#' Scales all features and performs PCA for dimensional reduction
#' 
#' @param seurat_obj A Seurat object (should be normalized and have variable features)
#' @param features Features to scale. If NULL, scales all features (default: NULL)
#' @param use_variable_features Logical, whether to use variable features for PCA (default: TRUE)
#' @param npcs Number of principal components to compute (default: 50)
#' @param verbose Print progress messages (default: TRUE)
#' @return A Seurat object with scaled data and PCA results
#' @examples
#' seurat_obj <- scale_and_reduce_dimensions(seurat_obj)
#' # To visualize elbow plot:
#' ElbowPlot(seurat_obj)
scale_and_reduce_dimensions <- function(seurat_obj,
                                        features = NULL,
                                        use_variable_features = TRUE,
                                        npcs = 50,
                                        verbose = TRUE) {
  
  if (verbose) cat("Scaling data and performing dimensional reduction...\n")
  
  # Scale the data
  if (is.null(features)) {
    if (verbose) cat("  Scaling all features...\n")
    all_features <- rownames(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features = all_features, verbose = verbose)
  } else {
    if (verbose) cat("  Scaling specified features...\n")
    seurat_obj <- ScaleData(seurat_obj, features = features, verbose = verbose)
  }
  
  # Determine which features to use for PCA
  if (use_variable_features) {
    pca_features <- VariableFeatures(object = seurat_obj)
    if (verbose) cat("  Using", length(pca_features), "variable features for PCA\n")
  } else {
    pca_features <- rownames(seurat_obj)
    if (verbose) cat("  Using all", length(pca_features), "features for PCA\n")
  }
  
  # Perform linear dimensional reduction
  if (verbose) cat("  Running PCA...\n")
  seurat_obj <- RunPCA(seurat_obj, features = pca_features, npcs = npcs, verbose = verbose)
  
  if (verbose) {
    cat("  PCA complete! Computed", npcs, "principal components\n")
    cat("  Use ElbowPlot(seurat_obj) to determine the dimensionality of the dataset\n")
    cat("Scaling and dimensional reduction complete!\n")
  }
  
  return(seurat_obj)
}


#' Split Cells into Training and Testing Sets
#' 
#' Splits cells into training (70%) and testing (30%) sets, stratified by cell type.
#' The split is performed separately for each cell type to ensure balanced representation.
#' 
#' @param seurat_obj A Seurat object with cell type assignments
#' @param column_name Name of the metadata column containing cell types (default: "cell_type")
#' @param train_proportion Proportion of cells to use for training (default: 0.7)
#' @param seed Random seed for reproducibility (default: 241)
#' @param verbose Print progress messages (default: TRUE)
#' @return A list containing:
#'   - `train_data`: Seurat object with training cells
#'   - `test_data`: Seurat object with testing cells
#'   - `train_cells`: Vector of training cell names
#'   - `test_cells`: Vector of testing cell names
#'   - `distribution`: List with train and test cell type distributions
#' @examples
#' split_result <- train_test_split_cells(seurat_obj, seed = 241)
#' train_data <- split_result$train_data
#' test_data <- split_result$test_data
train_test_split_cells <- function(seurat_obj,
                                   column_name = "cell_type",
                                   train_proportion = 0.7,
                                   seed = 241,
                                   verbose = TRUE) {
  
  if (verbose) cat("Splitting cells into training and testing sets...\n")
  
  # Check if column exists
  if (!column_name %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", column_name, "not found in metadata. Available columns:",
               paste(colnames(seurat_obj@meta.data), collapse = ", ")))
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialize vectors to hold training and testing cell names
  train_cells <- c()
  test_cells <- c()
  
  # Get unique cell types
  cell_types <- unique(seurat_obj@meta.data[[column_name]])
  cell_types <- cell_types[!is.na(cell_types)]
  
  if (verbose) cat("  Found", length(cell_types), "cell types\n")
  
  # Loop through each unique cell type and split into training and testing sets
  for (cell_type in cell_types) {
    # Get cells belonging to the current cell type
    cell_type_mask <- seurat_obj@meta.data[[column_name]] == cell_type
    cell_type_mask[is.na(cell_type_mask)] <- FALSE
    cell_type_cells <- colnames(seurat_obj)[cell_type_mask]
    
    # Calculate number of cells for training and testing
    num_train_cells <- round(train_proportion * length(cell_type_cells))
    num_test_cells <- length(cell_type_cells) - num_train_cells
    
    # Sample cells for training and testing
    train_cell_type_cells <- sample(cell_type_cells, size = num_train_cells)
    test_cell_type_cells <- setdiff(cell_type_cells, train_cell_type_cells)
    
    # Append to the training and testing vectors
    train_cells <- c(train_cells, train_cell_type_cells)
    test_cells <- c(test_cells, test_cell_type_cells)
    
    if (verbose) {
      cat("    ", cell_type, ": ", num_train_cells, " train, ", num_test_cells, " test\n", sep = "")
    }
  }
  
  # Create Seurat objects for training and testing sets
  if (verbose) cat("  Creating training and testing Seurat objects...\n")
  train_data <- subset(seurat_obj, cells = train_cells)
  test_data <- subset(seurat_obj, cells = test_cells)
  
  # Get cell type distributions
  train_cell_type_counts <- table(train_data@meta.data[[column_name]])
  test_cell_type_counts <- table(test_data@meta.data[[column_name]])
  
  # Create distribution data frames
  train_distribution_df <- as.data.frame(train_cell_type_counts)
  colnames(train_distribution_df) <- c("Cell_Type", "Cell_Count")
  test_distribution_df <- as.data.frame(test_cell_type_counts)
  colnames(test_distribution_df) <- c("Cell_Type", "Cell_Count")
  
  if (verbose) {
    cat("\nSplit Summary:\n")
    cat("  Training set: ", length(train_cells), " cells\n", sep = "")
    cat("  Testing set: ", length(test_cells), " cells\n", sep = "")
    cat("\nTraining set cell type distribution:\n")
    print(train_cell_type_counts)
    cat("\nTesting set cell type distribution:\n")
    print(test_cell_type_counts)
  }
  
  return(list(
    train_data = train_data,
    test_data = test_data,
    train_cells = train_cells,
    test_cells = test_cells,
    distribution = list(
      train_counts = train_cell_type_counts,
      test_counts = test_cell_type_counts,
      train_df = train_distribution_df,
      test_df = test_distribution_df
    )
  ))
}


#' Get Cell Type Distribution Summary
#' 
#' Returns a summary of cell type distribution for a Seurat object
#' 
#' @param seurat_obj A Seurat object with cell type assignments
#' @param column_name Name of the metadata column containing cell types (default: "cell_type")
#' @param return_df Logical, whether to return data frame (default: TRUE)
#' @param verbose Print summary (default: TRUE)
#' @return A data frame or table with cell type counts
#' @examples
#' dist_summary <- get_cell_type_distribution(seurat_obj)
get_cell_type_distribution <- function(seurat_obj,
                                       column_name = "cell_type",
                                       return_df = TRUE,
                                       verbose = TRUE) {
  
  # Check if column exists
  if (!column_name %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", column_name, "not found in metadata. Available columns:",
               paste(colnames(seurat_obj@meta.data), collapse = ", ")))
  }
  
  # Get cell type counts
  cell_type_counts <- table(seurat_obj@meta.data[[column_name]])
  
  if (verbose) {
    cat("Cell Type Distribution:\n")
    cat("======================\n")
    print(cell_type_counts)
    cat("\nTotal cells:", sum(cell_type_counts), "\n")
  }
  
  if (return_df) {
    distribution_df <- as.data.frame(cell_type_counts)
    colnames(distribution_df) <- c("Cell_Type", "Cell_Count")
    return(distribution_df)
  } else {
    return(cell_type_counts)
  }
}


#' Extract Expression Data and Labels
#' 
#' Extracts normalized expression matrices and cell type labels for training and testing sets.
#' Optionally finds variable features in the training set.
#' 
#' @param train_data Seurat object for training set
#' @param test_data Seurat object for testing set
#' @param column_name Name of the metadata column containing cell types (default: "cell_type")
#' @param find_variable_features Logical, whether to find variable features in training set (default: TRUE)
#' @param nfeatures Number of variable features to select (default: 500)
#' @param selection_method Method for selecting variable features (default: "vst")
#' @param use_variable_features Logical, whether to extract only variable features (default: TRUE)
#' @param layer Layer to extract from (default: "data" for normalized data)
#' @param verbose Print progress messages (default: TRUE)
#' @return A list containing:
#'   - `train_expr`: Expression matrix for training set (data frame)
#'   - `test_expr`: Expression matrix for testing set (data frame)
#'   - `train_labels`: Cell type labels for training set (factor)
#'   - `test_labels`: Cell type labels for testing set (factor)
#'   - `variable_features`: Vector of variable feature names (if found)
#' @examples
#' expr_data <- extract_expression_data(train_data, test_data, nfeatures = 500)
#' train_expr <- expr_data$train_expr
#' train_labels <- expr_data$train_labels
extract_expression_data <- function(train_data,
                                    test_data,
                                    column_name = "cell_type",
                                    features = NULL,
                                    find_variable_features = TRUE,
                                    nfeatures = 500,
                                    selection_method = "vst",
                                    use_variable_features = TRUE,
                                    layer = "data",
                                    verbose = TRUE) {
  
  if (verbose) cat("Extracting expression data and labels...\n")
  
  # Check if column exists in both objects
  for (obj_name in c("train_data", "test_data")) {
    obj <- get(obj_name)
    if (!column_name %in% colnames(obj@meta.data)) {
      stop(paste("Column", column_name, "not found in", obj_name, "metadata"))
    }
  }
  
  # Find variable features in training set if requested
  variable_features <- NULL
  
  if (!is.null(features)) {
    if (verbose) cat("  Using user-specified features...\n")
    variable_features <- features
  } else if (find_variable_features) {
    if (verbose) cat("  Finding variable features in training set...\n")
    train_data <- FindVariableFeatures(train_data,
                                       selection.method = selection_method,
                                       nfeatures = nfeatures,
                                       verbose = verbose)
    variable_features <- VariableFeatures(train_data)
    if (verbose) cat("  Found", length(variable_features), "variable features\n")
  }
  
  # Determine which features to extract
  if (use_variable_features && !is.null(variable_features)) {
    features_to_extract <- variable_features
    if (verbose) cat("  Extracting", length(features_to_extract), "variable features\n")
  } else {
    # Get common features between train and test
    common_features <- intersect(rownames(train_data), rownames(test_data))
    features_to_extract <- common_features
    if (verbose) cat("  Extracting all", length(features_to_extract), "common features\n")
  }
  
  # Extract normalized expression data
  if (verbose) cat("  Extracting expression matrices...\n")
  train_expr <- as.data.frame(GetAssayData(train_data, layer = layer)[features_to_extract, ])
  test_expr <- as.data.frame(GetAssayData(test_data, layer = layer)[features_to_extract, ])
  
  # Extract cell type labels
  if (verbose) cat("  Extracting cell type labels...\n")
  train_labels <- as.factor(train_data@meta.data[[column_name]])
  test_labels <- as.factor(test_data@meta.data[[column_name]])
  
  if (verbose) {
    cat("  Training set: ", ncol(train_expr), " cells, ", nrow(train_expr), " features\n", sep = "")
    cat("  Testing set: ", ncol(test_expr), " cells, ", nrow(test_expr), " features\n", sep = "")
    cat("Expression data extraction complete!\n")
  }
  
  return(list(
    train_expr = train_expr,
    test_expr = test_expr,
    train_labels = train_labels,
    test_labels = test_labels,
    features_used = variable_features
  ))
}

#' Train Random Forest Model
#' 
#' Trains a Random Forest model for cell type classification using expression data.
#' 
#' @param train_expr Expression matrix for training set (features as rows, cells as columns)
#' @param train_labels Cell type labels for training set (factor or character vector)
#' @param ntree Number of trees to grow (default: 100)
#' @param importance Logical, whether to calculate variable importance (default: TRUE)
#' @param seed Random seed for reproducibility (default: 241)
#' @param verbose Print progress messages (default: TRUE)
#' @param ... Additional arguments passed to randomForest()
#' @return A list containing:
#'   - `model`: Trained Random Forest model
#'   - `training_time`: Training time in seconds
#' @examples
#' rf_result <- train_random_forest(train_expr, train_labels, ntree = 100)
#' rf_model <- rf_result$model
train_random_forest <- function(train_expr,
                                train_labels,
                                ntree = 100,
                                importance = TRUE,
                                seed = 241,
                                verbose = TRUE,
                                ...) {
  
  if (verbose) cat("Training Random Forest model...\n")
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Ensure train_labels is a factor
  if (!is.factor(train_labels)) {
    train_labels <- as.factor(train_labels)
  }
  
  # Check dimensions
  if (ncol(train_expr) != length(train_labels)) {
    stop("Number of columns in train_expr must match length of train_labels")
  }
  
  # Transpose expression matrix (randomForest expects samples as rows, features as columns)
  if (verbose) cat("  Transposing expression matrix...\n")
  train_expr_transposed <- t(train_expr)
  
  if (verbose) {
    cat("  Training set: ", nrow(train_expr_transposed), " cells, ", ncol(train_expr_transposed), " features\n", sep = "")
    cat("  Number of trees: ", ntree, "\n", sep = "")
  }
  
  # Train the Random Forest model
  if (verbose) cat("  Training model (this may take a while)...\n")
  training_time <- system.time({
    rf_model <- randomForest(x = train_expr_transposed,
                             y = train_labels,
                             ntree = ntree,
                             importance = importance,
                             ...)
  })
  
  if (verbose) {
    cat("  Training complete!\n")
    cat("  Training time (in seconds):", round(training_time["elapsed"], 2), "\n")
  }
  
  return(list(
    model = rf_model,
    training_time = training_time["elapsed"]
  ))
}


#' Predict Cell Types Using Random Forest Model
#' 
#' Predicts cell types for test data using a trained Random Forest model.
#' Optionally stores predictions in Seurat object metadata.
#' 
#' @param rf_model Trained Random Forest model (from train_random_forest)
#' @param test_expr Expression matrix for test set (features as rows, cells as columns)
#' @param test_data Seurat object for test set (optional, for storing predictions)
#' @param column_name Name of the metadata column to store predictions (default: "predicted_cell_type")
#' @param return_probs Logical, whether to return class probabilities (default: FALSE)
#' @param verbose Print progress messages (default: TRUE)
#' @return A list containing:
#'   - `predicted_labels`: Predicted cell type labels (factor)
#'   - `predicted_probs`: Class probabilities matrix (if return_probs = TRUE)
#'   - `test_data`: Seurat object with predictions in metadata (if test_data provided)
#' @examples
#' predictions <- predict_cell_types(rf_model, test_expr, test_data)
#' predicted_labels <- predictions$predicted_labels
predict_cell_types <- function(rf_model,
                               test_expr,
                               test_data = NULL,
                               column_name = "predicted_cell_type",
                               return_probs = FALSE,
                               verbose = TRUE) {
  
  if (verbose) cat("Predicting cell types...\n")
  
  # Transpose expression matrix
  if (verbose) cat("  Transposing expression matrix...\n")
  test_expr_transposed <- t(test_expr)
  
  if (verbose) {
    cat("  Test set: ", nrow(test_expr_transposed), " cells, ", ncol(test_expr_transposed), " features\n", sep = "")
  }
  
  # Get class probabilities for each test cell
  if (verbose) cat("  Generating predictions...\n")
  predicted_probs <- predict(rf_model, newdata = test_expr_transposed, type = "prob")
  
  # For each test cell, get the class with the highest probability
  predicted_labels <- colnames(predicted_probs)[apply(predicted_probs, 1, which.max)]
  
  # Get factor levels from the model
  # Always define factor levels using TRAINING labels
  predicted_labels <- factor(
    predicted_labels,
    levels = levels(rf_model$y)
  )
  
  # Assign names to predicted labels (cell barcodes)
  names(predicted_labels) <- colnames(test_expr)
  
  # Store in Seurat metadata if test_data is provided
  if (!is.null(test_data)) {
    if (verbose) cat("  Storing predictions in Seurat metadata...\n")
    test_data@meta.data[[column_name]] <- predicted_labels[colnames(test_data)]
    if (verbose) cat("  Predictions stored in column:", column_name, "\n")
  }
  
  if (verbose) cat("Prediction complete!\n")
  
  # Prepare return list
  result <- list(
    predicted_labels = predicted_labels
  )
  
  if (return_probs) {
    result$predicted_probs <- predicted_probs
  }
  
  if (!is.null(test_data)) {
    result$test_data <- test_data
  }
  
  return(result)
}


#' Evaluate Cell Type Predictions
#' 
#' Evaluates predicted cell types by comparing with actual cell types.
#' Calculates confusion matrix and accuracy metrics.
#' 
#' @param actual_labels Actual cell type labels (factor or character vector)
#' @param predicted_labels Predicted cell type labels (factor or character vector)
#' @param test_data Seurat object with predictions (optional, for extracting labels)
#' @param actual_col Name of metadata column containing actual labels (default: "cell_type")
#' @param predicted_col Name of metadata column containing predicted labels (default: "predicted_cell_type")
#' @param verbose Print evaluation results (default: TRUE)
#' @return A list containing:
#'   - `confusion_matrix`: Confusion matrix table
#'   - `accuracy`: Overall accuracy (from confusion matrix)
#'   - `accuracy_manual`: Manual accuracy calculation (matching cells / total cells)
#'   - `accurate_cells_count`: Number of correctly predicted cells
#'   - `total_cells`: Total number of cells
#' @examples
#' evaluation <- evaluate_predictions(actual_labels, predicted_labels)
#' print(evaluation$confusion_matrix)
#' cat("Accuracy:", evaluation$accuracy, "\n")
evaluate_predictions <- function(actual_labels = NULL,
                                 predicted_labels = NULL,
                                 test_data = NULL,
                                 actual_col = "cell_type",
                                 predicted_col = "predicted_cell_type",
                                 verbose = TRUE) {
  
  if (verbose) cat("Evaluating predictions...\n")
  
  # Extract labels from test_data if provided
  if (!is.null(test_data)) {
    if (!actual_col %in% colnames(test_data@meta.data)) {
      stop(paste("Column", actual_col, "not found in test_data metadata"))
    }
    if (!predicted_col %in% colnames(test_data@meta.data)) {
      stop(paste("Column", predicted_col, "not found in test_data metadata"))
    }
    
    actual_labels <- test_data@meta.data[[actual_col]]
    predicted_labels <- test_data@meta.data[[predicted_col]]
    
    # Ensure they're factors with same levels
    all_levels <- unique(c(levels(as.factor(actual_labels)), levels(as.factor(predicted_labels))))
    actual_labels <- factor(actual_labels, levels = all_levels)
    predicted_labels <- factor(predicted_labels, levels = all_levels)
  } else {
    if (is.null(actual_labels) || is.null(predicted_labels)) {
      stop("Either provide test_data or both actual_labels and predicted_labels")
    }
    
    # Ensure they're factors with same levels
    all_levels <- unique(c(levels(as.factor(actual_labels)), levels(as.factor(predicted_labels))))
    actual_labels <- factor(actual_labels, levels = all_levels)
    predicted_labels <- factor(predicted_labels, levels = all_levels)
  }
  
  # Check lengths match
  if (length(actual_labels) != length(predicted_labels)) {
    stop("Length of actual_labels and predicted_labels must match")
  }
  
  # Create confusion matrix
  if (verbose) cat("  Creating confusion matrix...\n")
  confusion_matrix <- table(Predicted = predicted_labels, Actual = actual_labels)
  
  # Compute accuracy from confusion matrix
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  # Manual accuracy calculation (comparing each cell individually)
  accurate_cells_count <- sum(actual_labels == predicted_labels, na.rm = TRUE)
  total_cells <- length(actual_labels)
  accuracy_manual <- accurate_cells_count / total_cells
  
  if (verbose) {
    cat("\nEvaluation Results:\n")
    cat("===================\n")
    cat("Confusion Matrix:\n")
    print(confusion_matrix)
    cat("\nAccuracy Metrics:\n")
    cat("  Accuracy (from confusion matrix):", round(accuracy, 6), "\n")
    cat("  Accuracy (manual calculation):", round(accuracy_manual, 6), "\n")
    cat("  Accurate cells:", accurate_cells_count, "\n")
    cat("  Total cells:", total_cells, "\n")
    cat("\n")
  }
  
  return(list(
    confusion_matrix = confusion_matrix,
    accuracy = accuracy,
    accuracy_manual = accuracy_manual,
    accurate_cells_count = accurate_cells_count,
    total_cells = total_cells
  ))
}


#' Merge Isoform and Gene Predictions
#' 
#' Merges predictions from isoform-based and gene-based Seurat objects into a single data frame.
#' 
#' @param test_data_isoforms Seurat object with isoform-based predictions
#' @param test_data_genes Seurat object with gene-based predictions
#' @param isoform_pred_col Name of the column containing isoform predictions (default: "predicted_label_global_isoforms")
#' @param gene_pred_col Name of the column containing gene predictions (default: "predicted_label_global_genes")
#' @param true_labels_col Name of the column containing true labels (default: "cell_type")
#' @param use_isoform_labels Logical, whether to use isoform object for true labels (default: TRUE)
#' @param verbose Print progress messages (default: TRUE)
#' @return A data frame with merged metadata containing true labels and both predictions
#' @examples
#' merged_meta <- merge_isoform_gene_predictions(test_data, test_data_genes)
merge_isoform_gene_predictions <- function(test_data_isoforms,
                                           test_data_genes,
                                           isoform_pred_col = "predicted_label_global_isoforms",
                                           gene_pred_col = "predicted_label_global_genes",
                                           true_labels_col = "cell_type",
                                           use_isoform_labels = TRUE,
                                           verbose = TRUE) {
  
  if (verbose) cat("Merging isoform and gene predictions...\n")
  
  # Check if prediction columns exist
  if (!isoform_pred_col %in% colnames(test_data_isoforms@meta.data)) {
    stop(paste("Column", isoform_pred_col, "not found in test_data_isoforms metadata"))
  }
  if (!gene_pred_col %in% colnames(test_data_genes@meta.data)) {
    stop(paste("Column", gene_pred_col, "not found in test_data_genes metadata"))
  }
  
  # ---- Align by cell barcodes ----
  common_barcodes <- intersect(
    rownames(test_data_isoforms@meta.data),
    rownames(test_data_genes@meta.data)
  )
  
  if (length(common_barcodes) == 0) {
    stop("No overlapping cell barcodes found between isoform and gene objects")
  }
  
  if (verbose) {
    cat("  Common cells:", length(common_barcodes), "\n")
  }
  
  # Extract predictions
  if (verbose) cat("  Extracting predictions from isoform object...\n")
  isoform_pred <- test_data_isoforms@meta.data[common_barcodes, isoform_pred_col, drop = FALSE]
  
  if (verbose) cat("  Extracting predictions from gene object...\n")
  gene_pred <- test_data_genes@meta.data[common_barcodes, gene_pred_col, drop = FALSE]
  
  # Extract true labels from both objects
  if (!true_labels_col %in% colnames(test_data_isoforms@meta.data)) {
    stop(paste("Column", true_labels_col, "not found in test_data_isoforms metadata"))
  }
  if (!true_labels_col %in% colnames(test_data_genes@meta.data)) {
    stop(paste("Column", true_labels_col, "not found in test_data_genes metadata"))
  }
  
  if (verbose) cat("  Extracting true labels from both isoform and gene objects...\n")
  
  true_labels_isoform <- test_data_isoforms@meta.data[
    common_barcodes, true_labels_col, drop = FALSE
  ]
  colnames(true_labels_isoform) <- "cell_type_from_isoform"
  
  true_labels_gene <- test_data_genes@meta.data[
    common_barcodes, true_labels_col, drop = FALSE
  ]
  colnames(true_labels_gene) <- "cell_type_from_gene"
  
  # Merge metadata into one dataframe with cell barcodes as rownames
  if (verbose) cat("  Merging metadata...\n")
  merged_meta <- cbind(
    true_labels_isoform,
    true_labels_gene,
    isoform_pred,
    gene_pred
  )
  
  if (verbose) {
    cat("  Merged metadata contains", nrow(merged_meta), "cells\n")
    cat("Merging complete!\n")
  }
  
  return(merged_meta)
}


#' Identify Overlapping Cells
#' 
#' Identifies cells where isoform and gene predictions agree (overlapping) or disagree (non-overlapping).
#' Returns subset Seurat objects for each category.
#' 
#' @param test_data_isoforms Seurat object with isoform-based predictions
#' @param test_data_genes Seurat object with gene-based predictions
#' @param merged_meta Merged metadata data frame (from merge_isoform_gene_predictions). If NULL, will be created automatically.
#' @param isoform_pred_col Name of the column containing isoform predictions (default: "predicted_label_global_isoforms")
#' @param gene_pred_col Name of the column containing gene predictions (default: "predicted_label_global_genes")
#' @param predicted_global_col Name of the column to store agreed predictions in overlapping cells (default: "predicted_global")
#' @param verbose Print progress messages (default: TRUE)
#' @return A list containing:
#'   - `overlap_data`: Seurat object with overlapping cells (same prediction from both)
#'   - `non_overlap_isoforms`: Seurat object with non-overlapping cells from isoform object
#'   - `non_overlap_genes`: Seurat object with non-overlapping cells from gene object
#'   - `overlap_barcodes`: Vector of overlapping cell barcodes
#'   - `non_overlap_barcodes`: Vector of non-overlapping cell barcodes
#'   - `merged_meta`: Merged metadata data frame
#' @examples
#' overlap_result <- identify_overlapping_cells(test_data, test_data_genes)
#' global_overlap <- overlap_result$overlap_data
#' global_non_overlap_isoforms <- overlap_result$non_overlap_isoforms
identify_overlapping_cells <- function(test_data_isoforms,
                                       test_data_genes,
                                       merged_meta = NULL,
                                       isoform_pred_col = "predicted_label_global_isoforms",
                                       gene_pred_col = "predicted_label_global_genes",
                                       predicted_global_col = "predicted_subset",
                                       verbose = TRUE) {
  
  if (verbose) cat("Identifying overlapping and non-overlapping cells...\n")
  
  # Merge metadata if not provided
  if (is.null(merged_meta)) {
    if (verbose) cat("  Merging predictions (merged_meta not provided)...\n")
    merged_meta <- merge_isoform_gene_predictions(test_data_isoforms,
                                                  test_data_genes,
                                                  isoform_pred_col = isoform_pred_col,
                                                  gene_pred_col = gene_pred_col,
                                                  verbose = verbose)
  }
  
  # Check if required columns exist in merged_meta
  if (!isoform_pred_col %in% colnames(merged_meta)) {
    stop(paste("Column", isoform_pred_col, "not found in merged_meta"))
  }
  if (!gene_pred_col %in% colnames(merged_meta)) {
    stop(paste("Column", gene_pred_col, "not found in merged_meta"))
  }
  
  # Identify overlapping and non-overlapping barcodes
  if (verbose) cat("  Identifying overlapping cells (same prediction from both)...\n")
  overlap_barcodes <- rownames(merged_meta)[merged_meta[[isoform_pred_col]] == merged_meta[[gene_pred_col]]]
  non_overlap_barcodes <- rownames(merged_meta)[merged_meta[[isoform_pred_col]] != merged_meta[[gene_pred_col]]]
  
  # Remove NA comparisons
  overlap_barcodes <- overlap_barcodes[!is.na(overlap_barcodes)]
  non_overlap_barcodes <- non_overlap_barcodes[!is.na(non_overlap_barcodes)]
  
  if (verbose) {
    cat("    Overlapping cells:", length(overlap_barcodes), "\n")
    cat("    Non-overlapping cells:", length(non_overlap_barcodes), "\n")
  }
  
  # Overlapping cells ----
  if (length(overlap_barcodes) > 0) {
    if (verbose) cat("  Creating Seurat object with overlapping cells...\n")
    
    overlap_data <- subset(test_data_isoforms, cells = overlap_barcodes)
    overlap_data@meta.data[[predicted_global_col]] <-
      merged_meta[overlap_barcodes, gene_pred_col]
  } else {
    if (verbose) cat("  No overlapping cells found.\n")
    overlap_data <- NULL
  }
  
  # Non-overlapping cells ----
  if (length(non_overlap_barcodes) > 0) {
    if (verbose) cat("  Creating Seurat objects with non-overlapping cells...\n")
    
    non_overlap_isoforms <- subset(test_data_isoforms, cells = non_overlap_barcodes)
    non_overlap_genes <- subset(test_data_genes, cells = non_overlap_barcodes)
  } else {
    if (verbose) cat("  No non-overlapping cells found. Skipping subset.\n")
    
    non_overlap_isoforms <- NULL
    non_overlap_genes <- NULL
  }
  
  if (verbose) {
    if (!is.null(overlap_data))
      cat("    Overlap Seurat object:", ncol(overlap_data), "cells\n")
    if (!is.null(non_overlap_isoforms))
      cat("    Non-overlap isoforms Seurat object:", ncol(non_overlap_isoforms), "cells\n")
    if (!is.null(non_overlap_genes))
      cat("    Non-overlap genes Seurat object:", ncol(non_overlap_genes), "cells\n")
    cat("Identification complete!\n")
  }
  
  return(list(
    overlap_data = overlap_data,
    non_overlap_isoforms = non_overlap_isoforms,
    non_overlap_genes = non_overlap_genes,
    overlap_barcodes = overlap_barcodes,
    non_overlap_barcodes = non_overlap_barcodes,
    merged_meta = merged_meta
  ))
}


#' Evaluate Overlapping Predictions
#' 
#' Evaluates prediction accuracy for cells where isoform and gene predictions agree (overlapping cells).
#' 
#' @param overlap_data Seurat object with overlapping cells (from identify_overlapping_cells)
#' @param actual_col Name of metadata column containing actual labels (default: "cell_type")
#' @param predicted_col Name of metadata column containing predicted labels (default: "predicted_global")
#' @param verbose Print evaluation results (default: TRUE)
#' @return A list containing:
#'   - `confusion_matrix`: Confusion matrix table
#'   - `accuracy`: Overall accuracy (from confusion matrix)
#'   - `accuracy_manual`: Manual accuracy calculation (matching cells / total cells)
#'   - `accurate_cells_count`: Number of correctly predicted cells
#'   - `total_cells`: Total number of cells
#' @examples
#' overlap_eval <- evaluate_overlapping_predictions(global_overlap)
#' print(overlap_eval$confusion_matrix)
#' cat("Accuracy:", overlap_eval$accuracy, "\n")
evaluate_overlapping_predictions <- function(overlap_data,
                                             actual_col = "cell_type",
                                             predicted_col = "predicted_global",
                                             verbose = TRUE) {
  
  if (verbose) cat("Evaluating overlapping predictions...\n")
  
  # Check if columns exist
  if (!actual_col %in% colnames(overlap_data@meta.data)) {
    stop(paste("Column", actual_col, "not found in overlap_data metadata"))
  }
  if (!predicted_col %in% colnames(overlap_data@meta.data)) {
    stop(paste("Column", predicted_col, "not found in overlap_data metadata"))
  }
  
  # Extract actual and predicted labels from metadata
  if (verbose) cat("  Extracting labels...\n")
  actual_labels <- overlap_data@meta.data[[actual_col]]
  predicted_labels <- overlap_data@meta.data[[predicted_col]]
  
  # Ensure they're factors with same levels
  all_levels <- unique(c(levels(as.factor(actual_labels)), levels(as.factor(predicted_labels))))
  actual_labels <- factor(actual_labels, levels = all_levels)
  predicted_labels <- factor(predicted_labels, levels = all_levels)
  
  # Create confusion matrix
  if (verbose) cat("  Creating confusion matrix...\n")
  confusion_matrix <- table(Predicted = predicted_labels, Actual = actual_labels)
  
  # Compute accuracy from confusion matrix
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  # Manual accuracy calculation (comparing each cell individually)
  accurate_cells_count <- sum(actual_labels == predicted_labels, na.rm = TRUE)
  total_cells <- length(actual_labels)
  accuracy_manual <- accurate_cells_count / total_cells
  
  if (verbose) {
    cat("\nEvaluation Results for Overlapping Cells:\n")
    cat("==========================================\n")
    cat("Confusion Matrix:\n")
    print(confusion_matrix)
    cat("\nAccuracy Metrics:\n")
    cat("  Accuracy (from confusion matrix):", round(accuracy, 6), "\n")
    cat("  Accuracy (manual calculation):", round(accuracy_manual, 6), "\n")
    cat("  Accurate cells:", accurate_cells_count, "\n")
    cat("  Total cells:", total_cells, "\n")
    cat("\n")
  }
  
  return(list(
    confusion_matrix = confusion_matrix,
    accuracy = accuracy,
    accuracy_manual = accuracy_manual,
    accurate_cells_count = accurate_cells_count,
    total_cells = total_cells
  ))
}

#' Combine Final Predictions from Multiple Sources
#' 
#' Combines predictions from global overlap, subset overlap, and subset non-overlap
#' Seurat objects to create a final combined evaluation. Calculates overall accuracy
#' using both manual calculation and confusion matrix.
#' 
#' @param global_overlap Seurat object with global overlapping cells (from identify_overlapping_cells)
#' @param subset_overlap Seurat object with subset overlapping cells
#' @param subset_non_overlap_isoforms Seurat object with subset non-overlapping cells (isoform-based)
#' @param global_pred_col Name of the column containing global predictions (default: "predicted_global")
#' @param subset_pred_col Name of the column containing subset overlap predictions (default: "predicted_subset")
#' @param subset_non_overlap_pred_col Name of the column containing subset non-overlap predictions (default: "predicted_label_subset_isoforms")
#' @param actual_col Name of the metadata column containing actual labels (default: "cell_type")
#' @param verbose Print progress messages and results (default: TRUE)
#' @return A list containing:
#'   - `final_predictions_df`: Combined data frame with all predictions
#'   - `confusion_matrix`: Confusion matrix table
#'   - `accuracy`: Overall accuracy (from confusion matrix)
#'   - `accuracy_manual`: Manual accuracy calculation (matching cells / total cells)
#'   - `accurate_cells_count`: Number of correctly predicted cells
#'   - `total_cells`: Total number of cells
#' @examples
#' final_eval <- combined_final_predictions(global_overlap,
#'                                          subset_overlap,
#'                                          subset_non_overlap_isoforms)
#' cat("Final Combined Accuracy:", final_eval$accuracy, "\n")
combined_final_predictions <- function(global_overlap,
                                       subset_overlap,
                                       subset_non_overlap_isoforms,
                                       global_pred_col = "predicted_global",
                                       subset_pred_col = "predicted_subset",
                                       subset_non_overlap_pred_col = "predicted_label_subset_isoforms",
                                       actual_col = "cell_type",
                                       verbose = TRUE) {
  
  if (verbose) cat("Combining final predictions from multiple sources...\n")
  
  final_dfs <- list()
  
  ## ---- Global overlap ----
  if (!is.null(global_overlap)) {
    
    if (!global_pred_col %in% colnames(global_overlap@meta.data)) {
      stop(paste("Column", global_pred_col, "not found in global_overlap metadata"))
    }
    if (!actual_col %in% colnames(global_overlap@meta.data)) {
      stop(paste("Column", actual_col, "not found in global_overlap metadata"))
    }
    
    global_overlap@meta.data$predicted_label <-
      global_overlap@meta.data[[global_pred_col]]
    
    final_dfs[["global"]] <- data.frame(
      barcode = colnames(global_overlap),
      cell_type = global_overlap@meta.data[[actual_col]],
      predicted_label = global_overlap@meta.data$predicted_label,
      stringsAsFactors = FALSE
    )
    
    if (verbose) cat("  Added global overlap:", nrow(final_dfs[["global"]]), "cells\n")
  }
  
  ## ---- Subset overlap ----
  if (!is.null(subset_overlap)) {
    
    if (!subset_pred_col %in% colnames(subset_overlap@meta.data)) {
      stop(paste("Column", subset_pred_col, "not found in subset_overlap metadata"))
    }
    if (!actual_col %in% colnames(subset_overlap@meta.data)) {
      stop(paste("Column", actual_col, "not found in subset_overlap metadata"))
    }
    
    subset_overlap@meta.data$predicted_label <-
      subset_overlap@meta.data[[subset_pred_col]]
    
    final_dfs[["subset_overlap"]] <- data.frame(
      barcode = colnames(subset_overlap),
      cell_type = subset_overlap@meta.data[[actual_col]],
      predicted_label = subset_overlap@meta.data$predicted_label,
      stringsAsFactors = FALSE
    )
    
    if (verbose) cat("  Added subset overlap:", nrow(final_dfs[["subset_overlap"]]), "cells\n")
  }
  
  ## ---- Subset non-overlap (isoforms) ----
  if (!is.null(subset_non_overlap_isoforms)) {
    
    if (!subset_non_overlap_pred_col %in% colnames(subset_non_overlap_isoforms@meta.data)) {
      stop(paste("Column", subset_non_overlap_pred_col,
                 "not found in subset_non_overlap_isoforms metadata"))
    }
    if (!actual_col %in% colnames(subset_non_overlap_isoforms@meta.data)) {
      stop(paste("Column", actual_col,
                 "not found in subset_non_overlap_isoforms metadata"))
    }
    
    subset_non_overlap_isoforms@meta.data$predicted_label <-
      subset_non_overlap_isoforms@meta.data[[subset_non_overlap_pred_col]]
    
    final_dfs[["subset_non_overlap"]] <- data.frame(
      barcode = colnames(subset_non_overlap_isoforms),
      cell_type = subset_non_overlap_isoforms@meta.data[[actual_col]],
      predicted_label = subset_non_overlap_isoforms@meta.data$predicted_label,
      stringsAsFactors = FALSE
    )
    
    if (verbose) cat("  Added subset non-overlap:", 
                     nrow(final_dfs[["subset_non_overlap"]]), "cells\n")
  }
  
  ## ---- Safety check ----
  if (length(final_dfs) == 0) {
    stop("No valid prediction objects provided. All inputs are NULL.")
  }
  
  ## ---- Combine all predictions ----
  final_predictions_df <- do.call(rbind, final_dfs)
  
  if (verbose) {
    cat("  Total combined cells:", nrow(final_predictions_df), "\n")
    cat("  Missing actual labels:",
        sum(is.na(final_predictions_df$cell_type)), "\n")
    cat("  Missing predicted labels:",
        sum(is.na(final_predictions_df$predicted_label)), "\n")
  }
  
  ## ---- Remove incomplete rows ----
  final_predictions_df_complete <- final_predictions_df[
    !is.na(final_predictions_df$cell_type) &
      !is.na(final_predictions_df$predicted_label), ]
  
  # Manually calculate accuracy
  if (verbose) cat("  Calculating accuracy...\n")
  accurate_cells_count <- sum(final_predictions_df_complete$cell_type == final_predictions_df_complete$predicted_label, na.rm = TRUE)
  total_cells <- nrow(final_predictions_df_complete)
  accuracy_manual <- accurate_cells_count / total_cells
  
  # Ensure all levels match for confusion matrix
  all_labels <- sort(unique(c(final_predictions_df_complete$cell_type, 
                              final_predictions_df_complete$predicted_label)))
  
  # Create confusion matrix
  if (verbose) cat("  Creating confusion matrix...\n")
  confusion_matrix_final <- table(
    Predicted = factor(final_predictions_df_complete$predicted_label, levels = all_labels),
    Actual = factor(final_predictions_df_complete$cell_type, levels = all_labels)
  )
  
  # Calculate accuracy from confusion matrix
  accuracy_final <- sum(diag(confusion_matrix_final)) / sum(confusion_matrix_final)
  
  if (verbose) {
    cat("\nFinal Combined Evaluation Results:\n")
    cat("===================================\n")
    cat("Total cells:", total_cells, "\n")
    cat("Accurate cells:", accurate_cells_count, "\n")
    cat("Accuracy (manual calculation):", round(accuracy_manual, 6), "\n")
    cat("Accuracy (from confusion matrix):", round(accuracy_final, 6), "\n")
    cat("\nConfusion Matrix:\n")
    print(confusion_matrix_final)
    cat("\n")
  }
  
  return(list(
    final_predictions_df = final_predictions_df,
    confusion_matrix = confusion_matrix_final,
    accuracy = accuracy_final,
    accuracy_manual = accuracy_manual,
    accurate_cells_count = accurate_cells_count,
    total_cells = total_cells
  ))
}



