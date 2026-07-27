# Single-Cell RNA-seq Analysis Functions
# Supporting functions for the ICON pipeline.
# Note: library() calls are NOT used here — dependencies are declared in DESCRIPTION.

#' Import Seurat Object
#'
#' Validates that the input is a Seurat object and prints basic summary information.
#'
#' @param file_path A Seurat object.
#' @return A Seurat object.
#' @export
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
#' Performs normalization and finds variable features.
#'
#' @param seurat_obj A Seurat object.
#' @param normalization_method Method for normalization: \code{"LogNormalize"}
#'   (default), \code{"CLR"}, or \code{"RC"}.
#' @param scale_factor Scale factor for normalization (default: \code{10000}).
#' @param nfeatures Number of variable features to select (default: \code{2000}).
#' @param selection_method Method for selecting variable features: \code{"vst"}
#'   (default), \code{"mean.var.plot"}, or \code{"dispersion"}.
#' @param verbose Print progress messages (default: \code{TRUE}).
#' @return A preprocessed Seurat object.
#' @export
preprocess_data <- function(seurat_obj,
                            normalization_method = "LogNormalize",
                            scale_factor = 10000,
                            nfeatures = 2000,
                            selection_method = "vst",
                            verbose = TRUE) {

  if (verbose) cat("Starting preprocessing...\n")

  # Step 1: Normalize data
  if (verbose) cat("Normalizing data...\n")
  seurat_obj <- Seurat::NormalizeData(seurat_obj,
                                      normalization.method = normalization_method,
                                      scale.factor = scale_factor,
                                      verbose = verbose)

  # Step 2: Find variable features
  if (verbose) cat("Finding variable features...\n")
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj,
                                             selection.method = selection_method,
                                             nfeatures = nfeatures,
                                             verbose = verbose)

  if (verbose) {
    cat("Preprocessing complete!\n")
    cat("  Variable features found:", length(Seurat::VariableFeatures(seurat_obj)), "\n")
  }

  return(seurat_obj)
}


#' Assign Cell Types to Data
#'
#' Assigns cell types by loading annotations from a file and merging with Seurat
#' object metadata.
#'
#' @param seurat_obj A Seurat object.
#' @param annotation_file Path to a CSV or TSV file containing cell type
#'   annotations (optional if \code{annotation_df} is provided).
#' @param annotation_df Data frame containing cell type annotations (optional if
#'   \code{annotation_file} is provided).
#' @param barcode_col Name of the column containing cell barcodes in the
#'   annotation file (default: \code{"cell"}).
#' @param cell_type_col Name of the column containing cell types in the
#'   annotation file (default: \code{"donor_id"}).
#' @param column_name Name of the metadata column to store cell types
#'   (default: \code{"cell_type"}).
#' @param verbose Print progress messages (default: \code{TRUE}).
#' @return A Seurat object with cell type assignments in metadata.
#' @examples
#' \dontrun{
#' # Using file path
#' seurat_obj <- assign_cell_types(seurat_obj,
#'                                 annotation_file = "path/to/donor_ids.tsv")
#'
#' # Using data frame
#' annotations <- read.csv("donor_ids.tsv", sep = "\t", header = TRUE)
#' seurat_obj <- assign_cell_types(seurat_obj, annotation_df = annotations)
#' }
#' @export
assign_cell_types <- function(seurat_obj,
                              annotation_file = NULL,
                              annotation_df = NULL,
                              barcode_col = "cell",
                              cell_type_col = "donor_id",
                              column_name = "cell_type",
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

    donor_ids <- read.csv(annotation_file, header = TRUE, sep = sep,
                          stringsAsFactors = FALSE)
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
  seurat_obj@meta.data <- dplyr::left_join(seurat_obj@meta.data, cell_types,
                                           by = "barcode")

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
    n_annotated <- sum(!is.na(seurat_obj@meta.data[[column_name]]))
    cat("  Cells with annotations:", n_annotated, "out of", ncol(seurat_obj), "\n")
  }

  return(seurat_obj)
}


#' Check Unique Cell Types in Dataset
#'
#' Returns information about unique cell types in the dataset.
#'
#' @param seurat_obj A Seurat object with cell type assignments.
#' @param column_name Name of the metadata column containing cell types
#'   (default: \code{"cell_type"}).
#' @param return_table Logical, whether to return a table with counts
#'   (default: \code{TRUE}).
#' @param plot Logical, whether to create a bar plot (default: \code{FALSE}).
#' @return A data frame with cell type counts (if \code{return_table = TRUE}).
#' @export
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
    p <- ggplot2::ggplot(summary_df,
                         ggplot2::aes(x = reorder(Cell_Type, Count),
                                      y = Count, fill = Cell_Type)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Cell Type Distribution",
                    x = "Cell Type",
                    y = "Number of Cells") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
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
#' Calculates the top highly expressed features (isoforms/genes) per cell type
#' and returns a unique set of features across all cell types.
#'
#' @param seurat_obj A Seurat object with cell type assignments.
#' @param column_name Name of the metadata column containing cell types
#'   (default: \code{"cell_type"}).
#' @param top_n Number of top features to select per cell type
#'   (default: \code{100}).
#' @param verbose Print progress messages (default: \code{TRUE}).
#' @return A list containing:
#'   \describe{
#'     \item{top_features_list}{Named list of top features per cell type.}
#'     \item{unique_top_features}{Unique set of features across all cell types.}
#'     \item{n_unique_features}{Number of unique features.}
#'   }
#' @examples
#' \dontrun{
#' result <- calculate_relative_feature_usage(seurat_obj, top_n = 100)
#' top_features <- result$unique_top_features
#' }
#' @export
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
    seurat_subset <- isoform_subsets[[ct]]

    # Get normalized expression matrix
    norm_expr <- Seurat::GetAssayData(seurat_subset, layer = "data")

    # Sum expression per feature across all cells
    isoform_totals <- Matrix::rowSums(norm_expr)

    # Get top N highly expressed features
    top_isoforms <- names(sort(isoform_totals, decreasing = TRUE))[
      1:min(top_n, length(isoform_totals))
    ]

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
#' Scales all features and performs PCA for dimensional reduction.
#'
#' @param seurat_obj A Seurat object (should be normalized and have variable
#'   features).
#' @param features Features to scale. If \code{NULL}, scales all features
#'   (default: \code{NULL}).
#' @param use_variable_features Logical, whether to use variable features for
#'   PCA (default: \code{TRUE}).
#' @param npcs Number of principal components to compute (default: \code{50}).
#' @param verbose Print progress messages (default: \code{TRUE}).
#' @return A Seurat object with scaled data and PCA results.
#' @examples
#' \dontrun{
#' seurat_obj <- scale_and_reduce_dimensions(seurat_obj)
#' # To visualize elbow plot:
#' Seurat::ElbowPlot(seurat_obj)
#' }
#' @export
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
    seurat_obj <- Seurat::ScaleData(seurat_obj, features = all_features,
                                    verbose = verbose)
  } else {
    if (verbose) cat("  Scaling specified features...\n")
    seurat_obj <- Seurat::ScaleData(seurat_obj, features = features,
                                    verbose = verbose)
  }

  # Determine which features to use for PCA
  if (use_variable_features) {
    pca_features <- Seurat::VariableFeatures(object = seurat_obj)
    if (verbose) cat("  Using", length(pca_features), "variable features for PCA\n")
  } else {
    pca_features <- rownames(seurat_obj)
    if (verbose) cat("  Using all", length(pca_features), "features for PCA\n")
  }

  # Perform linear dimensional reduction
  if (verbose) cat("  Running PCA...\n")
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = pca_features, npcs = npcs,
                               verbose = verbose)

  if (verbose) {
    cat("  PCA complete! Computed", npcs, "principal components\n")
    cat("  Plotting PCA...\n")
    print(Seurat::DimPlot(seurat_obj, reduction = "pca"))
    cat("  Plotting Elbow Plot...\n")
    print(Seurat::ElbowPlot(seurat_obj))
    cat("  Use ElbowPlot(seurat_obj) to determine the dimensionality of the dataset\n")
    cat("Scaling and dimensional reduction complete!\n")
  }

  return(seurat_obj)
}

#' Prepare a Seurat Dataset for ICON
#'
#' Imports a count matrix, preprocesses the data, optionally assigns cell type
#' annotations, and performs dimensionality reduction.
#'
#' This function can be used for both training and testing datasets. Cell type
#' annotations are assigned only when \code{assign_cell_types = TRUE}.
#'
#' @param seurat_obj A Seurat object or a character string giving the path to
#'   an RDS file containing a Seurat object.
#' @param annotation_file Path to a CSV or TSV file containing cell type
#'   annotations. Required only if \code{assign_cell_types = TRUE}.
#' @param barcode_col Name of the column in the annotation file containing cell
#'   barcodes (default: \code{"barcode"}).
#' @param cell_type_col Name of the column in the annotation file containing
#'   cell type labels (default: \code{"cell_type"}).
#' @param nfeatures Number of highly variable features to identify during
#'   preprocessing (default: \code{2000}).
#' @param assign_cell_types Logical; whether to assign cell type annotations to
#'   the Seurat object (default: \code{TRUE}).
#' @param verbose Logical; whether to print progress messages (default:
#'   \code{TRUE}).
#'
#' @return A processed Seurat object.
#'
#' @examples
#' \dontrun{
#' # Training dataset
#' train_obj <- prepare_seurat_dataset(
#'   seurat_obj = train_obj,
#'   annotation_file = "celltypes.csv",
#'   assign_cell_types = TRUE
#' )
#'
#' # Testing dataset without annotations
#' test_obj <- prepare_seurat_dataset(
#'   seurat_obj = test_obj,
#'   assign_cell_types = FALSE
#' )
#' }
#'
#' @export
prepare_seurat_dataset <- function(
    seurat_obj,
    annotation_file = NULL,
    barcode_col = "barcode",
    cell_type_col = "cell_type",
    nfeatures = 2000,
    assign_cell_types = TRUE,
    verbose = TRUE) {
  
  if (is.character(seurat_obj)) {
    seurat_obj <- readRDS(seurat_obj)
  }
  
  seurat_obj <- import_count_matrix(seurat_obj)
  
  seurat_obj <- preprocess_data(
    seurat_obj,
    nfeatures = nfeatures,
    verbose = verbose
  )
  
  if (assign_cell_types) {
    
    if (is.null(annotation_file)) {
      stop(
        "'annotation_file' must be provided when ",
        "'assign_cell_types = TRUE'."
      )
    }
    
    seurat_obj <- assign_cell_types(
      seurat_obj,
      annotation_file = annotation_file,
      barcode_col = barcode_col,
      cell_type_col = cell_type_col,
      verbose = verbose
    )
  }
  
  seurat_obj <- scale_and_reduce_dimensions(
    seurat_obj,
    verbose = verbose
  )
  
  return(seurat_obj)
}

#' Validate Paired Isoform and Gene Datasets
#'
#' Checks whether paired isoform and gene Seurat objects contain identical
#' cell barcodes. This validation ensures that isoform-level and gene-level
#' datasets correspond to the same set of cells before downstream integration
#' and model training.
#'
#' @param isoform_data A Seurat object containing isoform-level expression data.
#' @param gene_data A Seurat object containing gene-level expression data.
#' @param dataset_name Character string describing the dataset being validated
#'   (e.g., \code{"training"} or \code{"testing"}).
#'   Default is \code{"dataset"}.
#' @param verbose Logical; whether to print validation messages.
#'   Default is \code{TRUE}.
#'
#' @return Invisibly returns \code{TRUE} if the cell barcodes match.
#'   Throws an error if the isoform and gene datasets contain different
#'   cell barcodes.
#'
#' @examples
#' \dontrun{
#' # Validate training datasets
#' validate_paired_datasets(
#'   isoform_data = isoform_train,
#'   gene_data = gene_train,
#'   dataset_name = "training"
#' )
#'
#' # Validate testing datasets
#' validate_paired_datasets(
#'   isoform_data = isoform_test,
#'   gene_data = gene_test,
#'   dataset_name = "testing"
#' )
#' }
#'
#' @export
validate_paired_datasets <- function(
    isoform_data,
    gene_data,
    dataset_name = "dataset",
    verbose = TRUE) {
  
  # Extract cell barcodes
  iso_cells <- colnames(isoform_data)
  gene_cells <- colnames(gene_data)
  
  # Check matching cell identities
  if (!setequal(iso_cells, gene_cells)) {
    
    missing_in_gene <- setdiff(iso_cells, gene_cells)
    missing_in_isoform <- setdiff(gene_cells, iso_cells)
    
    stop(
      paste0(
        "Isoform and gene ", dataset_name,
        " datasets contain different cell barcodes.\n",
        "Missing in gene dataset: ",
        paste(missing_in_gene, collapse = ", "),
        "\nMissing in isoform dataset: ",
        paste(missing_in_isoform, collapse = ", ")
      )
    )
  }
  
  if (verbose) {
    cat(
      dataset_name,
      "datasets validated:",
      length(iso_cells),
      "cells matched.\n"
    )
  }
  
  #invisible(TRUE)
}

#' Extract Expression Data and Labels
#'
#' Extracts normalized expression matrices and cell type labels for training
#' and testing datasets. Variable features are identified using only the
#' training dataset to avoid information leakage. Cell type labels are required
#' for training data but are optional for testing data, allowing prediction on
#' unlabeled datasets.
#'
#' @param train_data Seurat object containing the training dataset.
#'   The metadata must contain the cell type annotation column.
#' @param test_data Seurat object containing the testing dataset.
#'   Cell type annotations are optional and will only be extracted if available.
#' @param column_name Name of the metadata column containing cell types
#'   (default: \code{"cell_type"}).
#' @param features Optional character vector of specific features to extract.
#'   If provided, overrides variable feature selection.
#' @param find_variable_features Logical; whether to identify variable features
#'   using the training dataset (default: \code{TRUE}).
#' @param nfeatures Number of variable features to select from the training
#'   dataset (default: \code{500}).
#' @param selection_method Method used for variable feature selection
#'   (default: \code{"vst"}).
#' @param use_variable_features Logical; whether to extract only selected
#'   variable features. If \code{FALSE}, all common features between training
#'   and testing datasets are extracted (default: \code{TRUE}).
#' @param layer Assay layer used for expression extraction
#'   (default: \code{"data"}).
#' @param verbose Logical; whether to print progress messages
#'   (default: \code{TRUE}).
#'
#' @return A list containing:
#'   \describe{
#'     \item{train_expr}{Expression matrix for training dataset as a data frame.}
#'     \item{test_expr}{Expression matrix for testing dataset as a data frame.}
#'     \item{train_labels}{Cell type labels for training dataset as a factor.}
#'     \item{test_labels}{Cell type labels for testing dataset as a factor if
#'     available, otherwise \code{NULL}.}
#'     \item{features_used}{Character vector containing selected features.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Extract expression from labelled training data and unlabelled test data
#' expr_data <- extract_expression_data(
#'   train_data = train_seurat,
#'   test_data = test_seurat,
#'   nfeatures = 500
#' )
#'
#' train_expr <- expr_data$train_expr
#' train_labels <- expr_data$train_labels
#' predictions_input <- expr_data$test_expr
#' }
#'
#' @export
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
  
  if (verbose) {
    cat("Extracting expression data and labels...\n")
  }
  
  # Validate training labels (mandatory)
  
  if (!column_name %in% colnames(train_data@meta.data)) {
    
    stop(
      paste(
        "Column", column_name,
        "not found in training dataset metadata."
      )
    )
  }
  

  # Check testing labels (optional)
  
  test_labels_available <- column_name %in%
    colnames(test_data@meta.data)
  
  if (verbose) {
    
    if (test_labels_available) {
      cat("  Testing cell type labels detected.\n")
    } else {
      cat("  Testing cell type labels not detected. Prediction mode enabled.\n")
    }
  }
  

  # Feature selection using training data only
  
  variable_features <- NULL
  
  if (!is.null(features)) {
    
    if (verbose) {
      cat("  Using user-specified features...\n")
    }
    
    variable_features <- features
    
  } else if (find_variable_features) {
    
    if (verbose) {
      cat("  Finding variable features in training set...\n")
    }
    
    train_data <- Seurat::FindVariableFeatures(
      train_data,
      selection.method = selection_method,
      nfeatures = nfeatures,
      verbose = verbose
    )
    
    variable_features <- Seurat::VariableFeatures(train_data)
    
    if (verbose) {
      cat(
        "  Found",
        length(variable_features),
        "variable features\n"
      )
    }
  }
  
  # Determine features for extraction
  
  if (use_variable_features && !is.null(variable_features)) {
    
    features_to_extract <- variable_features
    
    if (verbose) {
      cat(
        "  Extracting",
        length(features_to_extract),
        "variable features\n"
      )
    }
    
  } else {
    
    common_features <- intersect(
      rownames(train_data),
      rownames(test_data)
    )
    
    features_to_extract <- common_features
    
    if (verbose) {
      cat(
        "  Extracting",
        length(features_to_extract),
        "common features\n"
      )
    }
  }
  

  # Extract expression matrices
  
  if (verbose) {
    cat("  Extracting expression matrices...\n")
  }
  
  train_expr <- as.data.frame(
    Seurat::GetAssayData(
      train_data,
      layer = layer
    )[features_to_extract, ]
  )
  
  test_expr <- as.data.frame(
    Seurat::GetAssayData(
      test_data,
      layer = layer
    )[features_to_extract, ]
  )
  

  # Extract labels
  
  if (verbose) {
    cat("  Extracting cell type labels...\n")
  }
  
  # Training labels are always required
  train_labels <- as.factor(
    train_data@meta.data[[column_name]]
  )
  
  
  # Testing labels are optional
  if (test_labels_available) {
    
    test_labels <- as.factor(
      test_data@meta.data[[column_name]]
    )
    
  } else {
    
    test_labels <- NULL
    
  }

  # Summary output
  
  if (verbose) {
    
    cat(
      "  Training set: ",
      ncol(train_expr),
      " cells, ",
      nrow(train_expr),
      " features\n",
      sep = ""
    )
    
    cat(
      "  Testing set: ",
      ncol(test_expr),
      " cells, ",
      nrow(test_expr),
      " features\n",
      sep = ""
    )
    
    if (!is.null(test_labels)) {
      cat(
        "  Testing labels available:",
        length(test_labels),
        "\n"
      )
    }
    
    cat("Expression data extraction complete!\n")
  }
  
  
  return(
    list(
      train_expr = train_expr,
      test_expr = test_expr,
      train_labels = train_labels,
      test_labels = test_labels,
      features_used = variable_features
    )
  )
}


#' Train Random Forest Model
#'
#' Trains a Random Forest model for cell type classification using expression
#' data.
#'
#' @param train_expr Expression matrix for the training set (features as rows,
#'   cells as columns).
#' @param train_labels Cell type labels for the training set (factor or
#'   character vector).
#' @param ntree Number of trees to grow (default: \code{100}).
#' @param importance Logical, whether to calculate variable importance
#'   (default: \code{TRUE}).
#' @param seed Random seed for reproducibility (default: \code{241}).
#' @param verbose Print progress messages (default: \code{TRUE}).
#' @param ... Additional arguments passed to \code{randomForest::randomForest()}.
#' @return A list containing:
#'   \describe{
#'     \item{model}{Trained Random Forest model.}
#'     \item{training_time}{Training time in seconds.}
#'   }
#' @examples
#' \dontrun{
#' rf_result <- train_random_forest(train_expr, train_labels, ntree = 100)
#' rf_model <- rf_result$model
#' }
#' @export
train_random_forest <- function(train_expr,
                                train_labels,
                                ntree = 100,
                                importance = TRUE,
                                seed = 241,
                                verbose = TRUE,
                                ...) {

  if (verbose) cat("Training Random Forest model...\n")

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
    cat("  Training set: ", nrow(train_expr_transposed), " cells, ",
        ncol(train_expr_transposed), " features\n", sep = "")
    cat("  Number of trees: ", ntree, "\n", sep = "")
  }

  if (verbose) cat("  Training model (this may take a while)...\n")
  training_time <- system.time({
    rf_model <- randomForest::randomForest(x = train_expr_transposed,
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
#' @param rf_model Trained Random Forest model (from \code{train_random_forest}).
#' @param test_expr Expression matrix for the test set (features as rows,
#'   cells as columns).
#' @param test_data Seurat object for the test set (optional, for storing
#'   predictions).
#' @param column_name Name of the metadata column to store predictions
#'   (default: \code{"predicted_cell_type"}).
#' @param return_probs Logical, whether to return class probabilities
#'   (default: \code{FALSE}).
#' @param verbose Print progress messages (default: \code{TRUE}).
#' @return A list containing:
#'   \describe{
#'     \item{predicted_labels}{Predicted cell type labels (factor).}
#'     \item{predicted_probs}{Class probabilities matrix (if
#'       \code{return_probs = TRUE}).}
#'     \item{test_data}{Seurat object with predictions in metadata (if
#'       \code{test_data} provided).}
#'   }
#' @examples
#' \dontrun{
#' predictions <- predict_cell_types(rf_model, test_expr, test_data)
#' predicted_labels <- predictions$predicted_labels
#' }
#' @export
predict_cell_types <- function(rf_model,
                               test_expr,
                               test_data = NULL,
                               column_name = "predicted_cell_type",
                               return_probs = FALSE,
                               verbose = TRUE) {

  if (verbose) cat("Predicting cell types...\n")

  if (verbose) cat("  Transposing expression matrix...\n")
  test_expr_transposed <- t(test_expr)

  if (verbose) {
    cat("  Test set: ", nrow(test_expr_transposed), " cells, ",
        ncol(test_expr_transposed), " features\n", sep = "")
  }

  if (verbose) cat("  Generating predictions...\n")
  predicted_probs <- stats::predict(rf_model, newdata = test_expr_transposed,
                                    type = "prob")

  # For each test cell, get the class with the highest probability
  predicted_labels <- colnames(predicted_probs)[apply(predicted_probs, 1, which.max)]

  # Always define factor levels using TRAINING labels
  predicted_labels <- factor(predicted_labels, levels = levels(rf_model$y))

  # Assign names to predicted labels (cell barcodes)
  names(predicted_labels) <- colnames(test_expr)

  # Store in Seurat metadata if test_data is provided
  if (!is.null(test_data)) {
    if (verbose) cat("  Storing predictions in Seurat metadata...\n")
    test_data@meta.data[[column_name]] <- predicted_labels[colnames(test_data)]
    if (verbose) cat("  Predictions stored in column:", column_name, "\n")
  }

  if (verbose) cat("Prediction complete!\n")

  result <- list(predicted_labels = predicted_labels)

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
#' @param actual_labels Actual cell type labels (factor or character vector).
#' @param predicted_labels Predicted cell type labels (factor or character
#'   vector).
#' @param test_data Seurat object with predictions (optional, for extracting
#'   labels).
#' @param actual_col Name of metadata column containing actual labels
#'   (default: \code{"cell_type"}).
#' @param predicted_col Name of metadata column containing predicted labels
#'   (default: \code{"predicted_cell_type"}).
#' @param verbose Print evaluation results (default: \code{TRUE}).
#' @return A list containing:
#'   \describe{
#'     \item{confusion_matrix}{Confusion matrix table.}
#'     \item{accuracy}{Overall accuracy (from confusion matrix).}
#'     \item{accuracy_manual}{Manual accuracy calculation.}
#'     \item{accurate_cells_count}{Number of correctly predicted cells.}
#'     \item{total_cells}{Total number of cells.}
#'   }
#' @examples
#' \dontrun{
#' evaluation <- evaluate_predictions(actual_labels, predicted_labels)
#' print(evaluation$confusion_matrix)
#' cat("Accuracy:", evaluation$accuracy, "\n")
#' }
#' @export
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

    all_levels <- unique(c(levels(as.factor(actual_labels)),
                           levels(as.factor(predicted_labels))))
    actual_labels <- factor(actual_labels, levels = all_levels)
    predicted_labels <- factor(predicted_labels, levels = all_levels)
  } else {
    if (is.null(actual_labels) || is.null(predicted_labels)) {
      stop("Either provide test_data or both actual_labels and predicted_labels")
    }

    all_levels <- unique(c(levels(as.factor(actual_labels)),
                           levels(as.factor(predicted_labels))))
    actual_labels <- factor(actual_labels, levels = all_levels)
    predicted_labels <- factor(predicted_labels, levels = all_levels)
  }

  if (length(actual_labels) != length(predicted_labels)) {
    stop("Length of actual_labels and predicted_labels must match")
  }

  if (verbose) cat("  Creating confusion matrix...\n")
  confusion_matrix <- table(Predicted = predicted_labels, Actual = actual_labels)

  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

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
    cat("  Total cells:", total_cells, "\n\n")
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
#' Merges cell type predictions from isoform-based and gene-based Seurat
#' objects into a single metadata data frame. If true cell type labels are
#' available, they are also included for evaluation purposes.
#'
#' @param test_data_isoforms Seurat object containing isoform-based
#'   predictions in its metadata.
#' @param test_data_genes Seurat object containing gene-based predictions
#'   in its metadata.
#' @param isoform_pred_col Name of the metadata column containing isoform-based
#'   predictions (default: \code{"predicted_label_global_isoforms"}).
#' @param gene_pred_col Name of the metadata column containing gene-based
#'   predictions (default: \code{"predicted_label_global_genes"}).
#' @param true_labels_col Optional name of the metadata column containing true
#'   cell type labels. Used only when evaluating datasets with known
#'   annotations (default: \code{"cell_type"}).
#' @param verbose Logical; whether to print progress messages
#'   (default: \code{TRUE}).
#'
#' @return A data frame containing:
#' \describe{
#'   \item{isoform predictions}{Predicted cell types from the isoform model.}
#'   \item{gene predictions}{Predicted cell types from the gene model.}
#'   \item{cell_type_from_isoform}{True labels from isoform object, if available.}
#'   \item{cell_type_from_gene}{True labels from gene object, if available.}
#' }
#'
#' @details
#' Cells are matched between isoform and gene datasets using shared cell
#' barcodes. True cell type labels are not required for merging predictions
#' and are included only when available.
#'
#' @examples
#' \dontrun{
#' # Merge predictions from datasets without true labels
#' merged_meta <- merge_isoform_gene_predictions(
#'   test_data_isoforms,
#'   test_data_genes,
#'   isoform_pred_col = "Global_isoform_predicted_cell_type",
#'   gene_pred_col = "Global_gene_predicted_cell_type"
#' )
#'
#' # Merge predictions with known cell types for benchmarking
#' merged_meta <- merge_isoform_gene_predictions(
#'   test_data_isoforms,
#'   test_data_genes,
#'   isoform_pred_col = "Global_isoform_predicted_cell_type",
#'   gene_pred_col = "Global_gene_predicted_cell_type",
#'   true_labels_col = "cell_type"
#' )
#' }
#'
#' @export
merge_isoform_gene_predictions <- function(test_data_isoforms,
                                           test_data_genes,
                                           isoform_pred_col = "predicted_label_global_isoforms",
                                           gene_pred_col = "predicted_label_global_genes",
                                           true_labels_col = "cell_type",
                                           verbose = TRUE) {
  
  if (verbose) cat("Merging isoform and gene predictions...\n")
  
  # Check prediction columns exist
  if (!isoform_pred_col %in% colnames(test_data_isoforms@meta.data)) {
    stop(paste(
      "Column", isoform_pred_col,
      "not found in test_data_isoforms metadata"
    ))
  }
  
  if (!gene_pred_col %in% colnames(test_data_genes@meta.data)) {
    stop(paste(
      "Column", gene_pred_col,
      "not found in test_data_genes metadata"
    ))
  }
  
  
  # Find common cells
  common_barcodes <- intersect(
    rownames(test_data_isoforms@meta.data),
    rownames(test_data_genes@meta.data)
  )
  
  if (length(common_barcodes) == 0) {
    stop("No overlapping cell barcodes found between isoform and gene objects")
  }
  
  if (verbose)
    cat("  Common cells:", length(common_barcodes), "\n")
  
  
  # Extract predictions
  if (verbose)
    cat("  Extracting isoform predictions...\n")
  
  isoform_pred <- test_data_isoforms@meta.data[
    common_barcodes,
    isoform_pred_col,
    drop = FALSE
  ]
  
  
  if (verbose)
    cat("  Extracting gene predictions...\n")
  
  gene_pred <- test_data_genes@meta.data[
    common_barcodes,
    gene_pred_col,
    drop = FALSE
  ]
  
  
  # Check whether true labels are available
  has_true_labels <- 
    true_labels_col %in% colnames(test_data_isoforms@meta.data) &&
    true_labels_col %in% colnames(test_data_genes@meta.data)
  
  
  if (has_true_labels) {
    
    if (verbose)
      cat("  True labels detected. Adding evaluation labels...\n")
    
    true_labels_isoform <- test_data_isoforms@meta.data[
      common_barcodes,
      true_labels_col,
      drop = FALSE
    ]
    
    colnames(true_labels_isoform) <- "cell_type_from_isoform"
    
    
    true_labels_gene <- test_data_genes@meta.data[
      common_barcodes,
      true_labels_col,
      drop = FALSE
    ]
    
    colnames(true_labels_gene) <- "cell_type_from_gene"
    
    
    merged_meta <- cbind(
      true_labels_isoform,
      true_labels_gene,
      isoform_pred,
      gene_pred
    )
    
  } else {
    
    if (verbose)
      cat("  True labels not available. Merging predictions only...\n")
    
    merged_meta <- cbind(
      isoform_pred,
      gene_pred
    )
  }
  
  
  if (verbose) {
    cat("  Merged metadata contains",
        nrow(merged_meta),
        "cells\n")
    cat("Merging complete!\n")
  }
  
  
  return(merged_meta)
}


#' Identify Overlapping Cells
#'
#' Identifies cells where isoform and gene predictions agree (overlapping) or
#' disagree (non-overlapping). Returns subset Seurat objects for each category.
#'
#' @param test_data_isoforms Seurat object with isoform-based predictions.
#' @param test_data_genes Seurat object with gene-based predictions.
#' @param merged_meta Merged metadata data frame (from
#'   \code{merge_isoform_gene_predictions}). If \code{NULL}, will be created
#'   automatically.
#' @param isoform_pred_col Name of the column containing isoform predictions
#'   (default: \code{"predicted_label_global_isoforms"}).
#' @param gene_pred_col Name of the column containing gene predictions
#'   (default: \code{"predicted_label_global_genes"}).
#' @param predicted_global_col Name of the column to store agreed predictions
#'   in overlapping cells (default: \code{"predicted_subset"}).
#' @param verbose Print progress messages (default: \code{TRUE}).
#' @return A list containing:
#'   \describe{
#'     \item{overlap_data}{Seurat object with overlapping cells.}
#'     \item{non_overlap_isoforms}{Seurat object with non-overlapping isoform cells.}
#'     \item{non_overlap_genes}{Seurat object with non-overlapping gene cells.}
#'     \item{overlap_barcodes}{Vector of overlapping cell barcodes.}
#'     \item{non_overlap_barcodes}{Vector of non-overlapping cell barcodes.}
#'     \item{merged_meta}{Merged metadata data frame.}
#'   }
#' @examples
#' \dontrun{
#' overlap_result <- identify_overlapping_cells(test_data, test_data_genes)
#' global_overlap <- overlap_result$overlap_data
#' global_non_overlap_isoforms <- overlap_result$non_overlap_isoforms
#' }
#' @export
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
    merged_meta <- merge_isoform_gene_predictions(
      test_data_isoforms,
      test_data_genes,
      isoform_pred_col = isoform_pred_col,
      gene_pred_col = gene_pred_col,
      verbose = verbose
    )
  }

  if (!isoform_pred_col %in% colnames(merged_meta)) {
    stop(paste("Column", isoform_pred_col, "not found in merged_meta"))
  }
  if (!gene_pred_col %in% colnames(merged_meta)) {
    stop(paste("Column", gene_pred_col, "not found in merged_meta"))
  }

  if (verbose) cat("  Identifying overlapping cells (same prediction from both)...\n")
  overlap_barcodes <- rownames(merged_meta)[
    merged_meta[[isoform_pred_col]] == merged_meta[[gene_pred_col]]
  ]
  non_overlap_barcodes <- rownames(merged_meta)[
    merged_meta[[isoform_pred_col]] != merged_meta[[gene_pred_col]]
  ]

  overlap_barcodes <- overlap_barcodes[!is.na(overlap_barcodes)]
  non_overlap_barcodes <- non_overlap_barcodes[!is.na(non_overlap_barcodes)]

  if (verbose) {
    cat("    Overlapping cells:", length(overlap_barcodes), "\n")
    cat("    Non-overlapping cells:", length(non_overlap_barcodes), "\n")
  }

  # Overlapping cells
  if (length(overlap_barcodes) > 0) {
    if (verbose) cat("  Creating Seurat object with overlapping cells...\n")
    overlap_data <- subset(test_data_isoforms, cells = overlap_barcodes)
    overlap_data@meta.data[[predicted_global_col]] <-
      merged_meta[overlap_barcodes, gene_pred_col]
  } else {
    if (verbose) cat("  No overlapping cells found.\n")
    overlap_data <- NULL
  }

  # Non-overlapping cells
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
      cat("    Non-overlap isoforms Seurat object:",
          ncol(non_overlap_isoforms), "cells\n")
    if (!is.null(non_overlap_genes))
      cat("    Non-overlap genes Seurat object:",
          ncol(non_overlap_genes), "cells\n")
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
#' Evaluates prediction accuracy for cells where isoform and gene predictions
#' agree (overlapping cells).
#'
#' @param overlap_data Seurat object with overlapping cells (from
#'   \code{identify_overlapping_cells}).
#' @param actual_col Name of metadata column containing actual labels
#'   (default: \code{"cell_type"}).
#' @param predicted_col Name of metadata column containing predicted labels
#'   (default: \code{"predicted_global"}).
#' @param verbose Print evaluation results (default: \code{TRUE}).
#' @return A list containing:
#'   \describe{
#'     \item{confusion_matrix}{Confusion matrix table.}
#'     \item{accuracy}{Overall accuracy (from confusion matrix).}
#'     \item{accuracy_manual}{Manual accuracy calculation.}
#'     \item{accurate_cells_count}{Number of correctly predicted cells.}
#'     \item{total_cells}{Total number of cells.}
#'   }
#' @examples
#' \dontrun{
#' overlap_eval <- evaluate_overlapping_predictions(global_overlap)
#' print(overlap_eval$confusion_matrix)
#' cat("Accuracy:", overlap_eval$accuracy, "\n")
#' }
#' @export
evaluate_overlapping_predictions <- function(overlap_data,
                                             actual_col = "cell_type",
                                             predicted_col = "predicted_global",
                                             verbose = TRUE) {

  if (verbose) cat("Evaluating overlapping predictions...\n")

  if (!actual_col %in% colnames(overlap_data@meta.data)) {
    stop(paste("Column", actual_col, "not found in overlap_data metadata"))
  }
  if (!predicted_col %in% colnames(overlap_data@meta.data)) {
    stop(paste("Column", predicted_col, "not found in overlap_data metadata"))
  }

  if (verbose) cat("  Extracting labels...\n")
  actual_labels <- overlap_data@meta.data[[actual_col]]
  predicted_labels <- overlap_data@meta.data[[predicted_col]]

  all_levels <- unique(c(levels(as.factor(actual_labels)),
                         levels(as.factor(predicted_labels))))
  actual_labels <- factor(actual_labels, levels = all_levels)
  predicted_labels <- factor(predicted_labels, levels = all_levels)

  if (verbose) cat("  Creating confusion matrix...\n")
  confusion_matrix <- table(Predicted = predicted_labels, Actual = actual_labels)

  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

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
    cat("  Total cells:", total_cells, "\n\n")
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
#' Combines predictions from global overlap, subset overlap, and subset
#' non-overlap Seurat objects to generate final cell type predictions.
#' If true cell type labels are available, the function calculates accuracy
#' and confusion matrix. For prediction-only datasets without true labels,
#' evaluation metrics are skipped.
#'
#' @param global_overlap Seurat object containing cells where global isoform
#'   and gene predictions agree.
#' @param subset_overlap Seurat object containing cells where subset isoform
#'   and gene predictions agree.
#' @param subset_non_overlap_isoforms Seurat object containing cells requiring
#'   subset isoform-based predictions.
#' @param global_pred_col Name of metadata column containing global predictions
#'   (default: \code{"predicted_global"}).
#' @param subset_pred_col Name of metadata column containing subset overlap
#'   predictions (default: \code{"predicted_subset"}).
#' @param subset_non_overlap_pred_col Name of metadata column containing
#'   subset non-overlap isoform predictions
#'   (default: \code{"predicted_label_subset_isoforms"}).
#' @param actual_col Name of metadata column containing true cell type labels
#'   (default: \code{"cell_type"}).
#' @param export_file Optional output Excel file path for saving final
#'   predictions (default: \code{NULL}).
#' @param verbose Print progress messages (default: \code{TRUE}).
#'
#' @return A list containing:
#'   \describe{
#'     \item{final_predictions_df}{Combined prediction data frame containing
#'     barcode, true labels (if available), and predicted labels.}
#'     \item{confusion_matrix}{Confusion matrix if true labels are available,
#'     otherwise \code{NULL}.}
#'     \item{accuracy}{Accuracy calculated from confusion matrix, or
#'     \code{NA} if true labels are unavailable.}
#'     \item{accuracy_manual}{Manually calculated accuracy, or
#'     \code{NA} if true labels are unavailable.}
#'     \item{accurate_cells_count}{Number of correctly predicted cells, or
#'     \code{NA} if labels are unavailable.}
#'     \item{total_cells}{Number of evaluated cells.}
#'     \item{prediction_file}{Exported prediction file path if provided.}
#'   }
#'
#' @examples
#' \dontrun{
#' final_eval <- combined_final_predictions(
#'   global_overlap,
#'   subset_overlap,
#'   subset_non_overlap_isoforms
#' )
#'
#' final_predictions <- final_eval$final_predictions_df
#' }
#'
#' @export
combined_final_predictions <- function(global_overlap,
                                       subset_overlap,
                                       subset_non_overlap_isoforms,
                                       global_pred_col = "predicted_global",
                                       subset_pred_col = "predicted_subset",
                                       subset_non_overlap_pred_col =
                                         "predicted_label_subset_isoforms",
                                       actual_col = "cell_type",
                                       export_file = NULL,
                                       verbose = TRUE) {
  
  
  if (verbose)
    cat("Combining final predictions from multiple sources...\n")
  
  
  final_dfs <- list()
  
  
  ## -------------------------------
  ## Global overlap predictions
  ## -------------------------------
  
  if (!is.null(global_overlap)) {
    
    if (!global_pred_col %in%
        colnames(global_overlap@meta.data)) {
      
      stop(paste("Column", global_pred_col,
                 "not found in global_overlap metadata"))
    }
    
    
    final_dfs[["global"]] <- data.frame(
      
      barcode = colnames(global_overlap),
      
      cell_type =
        if (actual_col %in%
            colnames(global_overlap@meta.data)) {
          
          global_overlap@meta.data[[actual_col]]
          
        } else {
          NA
        },
      
      predicted_label =
        global_overlap@meta.data[[global_pred_col]],
      
      stringsAsFactors = FALSE
    )
    
    
    if (verbose)
      cat("  Added global overlap:",
          nrow(final_dfs[["global"]]),
          "cells\n")
  }
  
  
  
  ## -------------------------------
  ## Subset overlap predictions
  ## -------------------------------
  
  if (!is.null(subset_overlap)) {
    
    
    if (!subset_pred_col %in%
        colnames(subset_overlap@meta.data)) {
      
      stop(paste("Column", subset_pred_col,
                 "not found in subset_overlap metadata"))
    }
    
    
    final_dfs[["subset_overlap"]] <- data.frame(
      
      barcode = colnames(subset_overlap),
      
      cell_type =
        if (actual_col %in%
            colnames(subset_overlap@meta.data)) {
          
          subset_overlap@meta.data[[actual_col]]
          
        } else {
          NA
        },
      
      
      predicted_label =
        subset_overlap@meta.data[[subset_pred_col]],
      
      stringsAsFactors = FALSE
    )
    
    
    if (verbose)
      cat("  Added subset overlap:",
          nrow(final_dfs[["subset_overlap"]]),
          "cells\n")
    
  }
  
  
  
  ## -------------------------------
  ## Subset non-overlap isoform
  ## -------------------------------
  
  if (!is.null(subset_non_overlap_isoforms)) {
    
    
    if (!subset_non_overlap_pred_col %in%
        colnames(subset_non_overlap_isoforms@meta.data)) {
      
      stop(paste(
        "Column",
        subset_non_overlap_pred_col,
        "not found in subset_non_overlap_isoforms metadata"
      ))
    }
    
    
    final_dfs[["subset_non_overlap"]] <- data.frame(
      
      barcode = colnames(subset_non_overlap_isoforms),
      
      cell_type =
        if (actual_col %in%
            colnames(subset_non_overlap_isoforms@meta.data)) {
          
          subset_non_overlap_isoforms@meta.data[[actual_col]]
          
        } else {
          
          NA
        },
      
      predicted_label =
        subset_non_overlap_isoforms@meta.data[[subset_non_overlap_pred_col]],
      
      stringsAsFactors = FALSE
    )
    
    
    if (verbose)
      cat("  Added subset non-overlap:",
          nrow(final_dfs[["subset_non_overlap"]]),
          "cells\n")
    
  }
  
  
  
  if (length(final_dfs) == 0) {
    stop("No prediction objects provided.")
  }
  
  
  
  ## Combine predictions
  
  final_predictions_df <- do.call(
    rbind,
    final_dfs
  )
  
  
  if (verbose) {
    
    cat("\nTotal combined cells:",
        nrow(final_predictions_df),
        "\n")
    
    cat("Missing true labels:",
        sum(is.na(final_predictions_df$cell_type)),
        "\n")
    
  }
  
  
  
  ## -------------------------------
  ## Evaluation
  ## -------------------------------
  
  if (all(is.na(final_predictions_df$cell_type))) {
    
    
    if (verbose) {
      
      cat("\nNo true cell type labels available.\n")
      cat("Skipping accuracy and confusion matrix calculation.\n")
      
    }
    
    
    confusion_matrix_final <- NULL
    accuracy_final <- NA
    accuracy_manual <- NA
    accurate_cells_count <- NA
    total_cells <- nrow(final_predictions_df)
    
    
  } else {
    
    
    evaluation_df <- final_predictions_df[
      !is.na(final_predictions_df$cell_type) &
        !is.na(final_predictions_df$predicted_label),
    ]
    
    
    accurate_cells_count <-
      sum(
        evaluation_df$cell_type ==
          evaluation_df$predicted_label
      )
    
    
    total_cells <- nrow(evaluation_df)
    
    
    accuracy_manual <-
      accurate_cells_count /
      total_cells
    
    
    labels <- sort(
      unique(
        c(
          evaluation_df$cell_type,
          evaluation_df$predicted_label
        )
      )
    )
    
    
    confusion_matrix_final <-
      table(
        Predicted =
          factor(
            evaluation_df$predicted_label,
            levels = labels
          ),
        
        Actual =
          factor(
            evaluation_df$cell_type,
            levels = labels
          )
      )
    
    
    accuracy_final <-
      sum(diag(confusion_matrix_final)) /
      sum(confusion_matrix_final)
    
    
    
    if (verbose) {
      
      cat("\nFinal Evaluation Results\n")
      cat("=======================\n")
      cat("Total evaluated cells:",
          total_cells,
          "\n")
      
      cat("Accuracy:",
          round(accuracy_final,6),
          "\n")
      
      cat("\nConfusion Matrix:\n")
      
      print(confusion_matrix_final)
      
    }
    
  }
  
  
  
  ## -------------------------------
  ## Export
  ## -------------------------------
  
  if (!is.null(export_file)) {
    
    
    if (verbose)
      cat("\nExporting predictions to:",
          export_file,
          "\n")
    
    
    wb <- openxlsx::createWorkbook()
    
    openxlsx::addWorksheet(
      wb,
      "Final_Predictions"
    )
    
    openxlsx::writeData(
      wb,
      "Final_Predictions",
      final_predictions_df
    )
    
    
    openxlsx::saveWorkbook(
      wb,
      export_file,
      overwrite = TRUE
    )
    
  }
  
  
  
  return(
    list(
      
      final_predictions_df =
        final_predictions_df,
      
      confusion_matrix =
        confusion_matrix_final,
      
      accuracy =
        accuracy_final,
      
      accuracy_manual =
        accuracy_manual,
      
      accurate_cells_count =
        accurate_cells_count,
      
      total_cells =
        total_cells,
      
      prediction_file =
        export_file
      
    )
  )
  
}

#' Generate UMAP Plot of Final Predictions
#'
#' Combines predictions from multiple sources (global overlap, subset overlap,
#' subset non-overlap) into a base Seurat object and generates a UMAP plot
#' coloured by predicted cell type.
#'
#' @param base_seurat A Seurat object to use as the base for the UMAP. Should
#'   contain all test cells.
#' @param global_overlap Seurat object with global overlapping cells (optional).
#' @param subset_overlap Seurat object with subset overlapping cells (optional).
#' @param subset_non_overlap_isoforms Seurat object with subset non-overlapping
#'   isoform cells (optional).
#' @param global_pred_col Column name for global predictions
#'   (default: \code{"predicted_global"}).
#' @param subset_pred_col Column name for subset predictions
#'   (default: \code{"predicted_subset"}).
#' @param subset_non_overlap_pred_col Column name for subset non-overlap
#'   predictions
#'   (default: \code{"Subset_isoform_predicted_cell_type"}).
#' @param reduction Dimensionality reduction to use for the plot
#'   (default: \code{"umap"}).
#' @param dims PCA dimensions to use if UMAP needs to be computed
#'   (default: \code{1:10}).
#' @param cols Optional named character vector of colours for cell types.
#' @param label Logical, whether to label clusters on the plot
#'   (default: \code{TRUE}).
#' @param verbose Print progress messages (default: \code{TRUE}).
#' @return A list containing:
#'   \describe{
#'     \item{seurat_object}{The base Seurat object with predictions added.}
#'     \item{umap_plot}{A ggplot2 UMAP plot coloured by predicted cell type.}
#'   }
#' @examples
#' \dontrun{
#' umap_result <- Final_prediction_UMAP(isoform_seurat_obj,
#'                                      global_overlap = global_overlap,
#'                                      subset_overlap = subset_overlap,
#'                                      subset_non_overlap_isoforms = subset_non_overlap_isoforms)
#' print(umap_result$umap_plot)
#' }
#' @export
Final_prediction_UMAP <- function(base_seurat,
                                  global_overlap = NULL,
                                  subset_overlap = NULL,
                                  subset_non_overlap_isoforms = NULL,
                                  global_pred_col = "predicted_global",
                                  subset_pred_col = "predicted_subset",
                                  subset_non_overlap_pred_col = "Subset_isoform_predicted_cell_type",
                                  reduction = "umap",
                                  dims = 1:10,
                                  cols = NULL,
                                  label = TRUE,
                                  verbose = TRUE) {

  ## Ensure UMAP exists
  if (!reduction %in% Seurat::Reductions(base_seurat)) {
    if (verbose) cat("UMAP not found. Running RunUMAP() using PCA...\n")
    base_seurat <- Seurat::RunUMAP(base_seurat, dims = dims)
  }

  ## Initialise prediction vector
  predicted_cell_type <- rep(NA_character_, ncol(base_seurat))
  names(predicted_cell_type) <- colnames(base_seurat)

  ## Level 1: Global
  if (!is.null(global_overlap)) {
    global_cells <- colnames(global_overlap)
    predicted_cell_type[global_cells] <-
      as.character(global_overlap@meta.data[[global_pred_col]])
  }

  ## Level 2: Subset overlap
  if (!is.null(subset_overlap)) {
    subset_cells <- colnames(subset_overlap)
    predicted_cell_type[subset_cells] <-
      as.character(subset_overlap@meta.data[[subset_pred_col]])
  }

  ## Level 3: Isoform subset
  if (!is.null(subset_non_overlap_isoforms)) {
    subset_iso_cells <- colnames(subset_non_overlap_isoforms)
    predicted_cell_type[subset_iso_cells] <-
      as.character(subset_non_overlap_isoforms@meta.data[[subset_non_overlap_pred_col]])
  }

  ## Add predictions to Seurat
  base_seurat$predicted_cell_type <- as.character(predicted_cell_type)

  if (verbose) {
    cat("Total cells:", ncol(base_seurat), "\n")
    cat("Cells with predictions:",
        sum(!is.na(base_seurat$predicted_cell_type)), "\n")
    cat("Cells without predictions:",
        sum(is.na(base_seurat$predicted_cell_type)), "\n")
    cat("Unique predicted cell types:",
        length(unique(stats::na.omit(base_seurat$predicted_cell_type))), "\n")
  }

  ## UMAP plot
  umap_plot <- Seurat::DimPlot(
    base_seurat,
    reduction = reduction,
    group.by = "predicted_cell_type",
    label = label,
    cols = cols
  )

  return(list(
    seurat_object = base_seurat,
    umap_plot = umap_plot
  ))
}

#' Extract Top Random Forest Features
#'
#' Extracts the top-ranking features for each cell type from a trained
#' Random Forest model based on variable importance scores.
#'
#' For each cell type, the function selects the top \code{top_n} features
#' ranked by their importance values and returns them as a data frame.
#'
#' @param rf_model A trained Random Forest model generated using
#'   \code{randomForest::randomForest()}.
#' @param top_n Number of top features to extract for each cell type
#'   (default: \code{10}).
#'
#' @return A data frame where each column corresponds to a cell type and
#'   each row contains one of the top-ranked features formatted as
#'   \code{"feature_name (importance_score)"}.
#'
#' @examples
#' \dontrun{
#' top_features <- get_top_rf_features(rf_model, top_n = 10)
#' head(top_features)
#' }
#'
#' @export
get_top_rf_features <- function(rf_model, top_n = 10) {
  
  # Check that top_n is valid
  if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1) {
    stop("'top_n' must be a positive integer.")
  }
  
  # Extract variable importance
  imp <- randomForest::importance(rf_model)
  
  # Remove overall importance columns
  cell_type_cols <- setdiff(
    colnames(imp),
    c("MeanDecreaseAccuracy", "MeanDecreaseGini")
  )
  
  # Helper function to extract top features for one cell type
  format_top_features <- function(cell_type) {
    
    sorted <- sort(imp[, cell_type], decreasing = TRUE)
    
    top_features <- head(sorted, top_n)
    
    paste0(
      names(top_features),
      " (",
      round(top_features, 3),
      ")"
    )
  }
  
  # Create summary table
  top_table <- data.frame(
    lapply(cell_type_cols, format_top_features),
    check.names = FALSE
  )
  
  colnames(top_table) <- cell_type_cols
  
  return(top_table)
}

#' Export Random Forest Top Features
#'
#' Combines available Random Forest top feature tables and exports them
#' into separate sheets of an Excel workbook.
#'
#' @param top_feature_list Named list containing top feature data frames
#'   generated by \code{get_top_rf_features()}.
#' @param file Output Excel filename.
#' @param verbose Print progress messages.
#'
#' @return Invisibly returns the output filename.
#'
#' @export
export_rf_top_features <- function(top_feature_list,
                                   file = "RF_Top_Features.xlsx",
                                   verbose = TRUE) {
  
  if (verbose)
    cat("Exporting Random Forest top features...\n")
  
  # Validate input
  if (!is.list(top_feature_list) || length(top_feature_list) == 0) {
    stop("'top_feature_list' must be a non-empty named list.")
  }
  
  if (is.null(names(top_feature_list)) ||
      any(names(top_feature_list) == "")) {
    stop("'top_feature_list' must contain named elements.")
  }
  
  wb <- openxlsx::createWorkbook()
  
  for (model_name in names(top_feature_list)) {
    
    top_features <- top_feature_list[[model_name]]
    
    # Skip unavailable models
    if (is.null(top_features)) {
      if (verbose)
        cat("Skipping:", model_name, "- not available\n")
      next
    }
    
    if (verbose)
      cat("Adding sheet:", model_name, "\n")
    
    # Excel worksheet names cannot exceed 31 characters
    sheet_name <- substr(model_name, 1, 31)
    
    openxlsx::addWorksheet(
      wb,
      sheet_name
    )
    
    openxlsx::writeData(
      wb,
      sheet_name,
      top_features
    )
  }
  
  openxlsx::saveWorkbook(
    wb,
    file,
    overwrite = TRUE
  )
  
  if (verbose)
    cat("Saved to:", file, "\n")
  
  invisible(file)
}