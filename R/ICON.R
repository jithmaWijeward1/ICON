#' ICON: Isoform-Cell-type classification via Overlapping Networks
#'
#' Runs the full hierarchical Random Forest pipeline for single-cell RNA-seq
#' cell type classification. Takes pre-split isoform and gene Seurat objects
#' for training and testing, along with a cell type annotation file, and
#' returns predicted cell types and Random Forest top feature lists. When test
#' cells contain cell type labels, accuracy metrics are also computed.
#'
#' @param isoform_train Seurat object (or path to an RDS file) containing the
#'   isoform count matrix for training cells.
#' @param isoform_test Seurat object (or path to an RDS file) containing the
#'   isoform count matrix for test cells.
#' @param gene_train Seurat object (or path to an RDS file) containing the
#'   gene count matrix for training cells.
#' @param gene_test Seurat object (or path to an RDS file) containing the
#'   gene count matrix for test cells.
#' @param annotation_file Path to a CSV, TSV, or TXT file containing cell type
#'   annotations. Must contain barcode and cell type columns.
#' @param barcode_col Name of the column in \code{annotation_file} containing
#'   cell barcodes (default: \code{"barcode"}).
#' @param cell_type_col Name of the column in \code{annotation_file} containing
#'   cell type labels (default: \code{"cell_type"}).
#' @param nfeatures Number of variable features for preprocessing
#'   (default: \code{2000}).
#' @param top_n Number of top features per cell type to select for relative
#'   feature usage (RIU/RGU) calculation (default: \code{100}).
#' @param nfeatures_model Number of variable features to use during model
#'   training (default: \code{500}).
#' @param ntree Number of trees in each Random Forest model (default: \code{100}).
#' @param seed Random seed for reproducibility (default: \code{241}).
#' @param export_predictions_file Optional path for exporting final predictions
#'   as an Excel workbook. Set to \code{NULL} to skip export.
#' @param export_features_file Optional path for exporting top RF features as
#'   an Excel workbook. Set to \code{NULL} to skip export.
#' @param assign_test_cell_types Logical; whether to assign cell type labels to
#'   test datasets from \code{annotation_file} (default: \code{FALSE}). Set to
#'   \code{TRUE} when test cells are annotated and accuracy should be calculated.
#' @param verbose Logical; whether to print progress messages (default: \code{TRUE}).
#'
#' @return A list containing:
#' \describe{
#'   \item{final_predictions_df}{Data frame of all test cells with predicted cell types.}
#'   \item{rf_top_features}{Named list of top Random Forest feature tables.}
#'   \item{global_isoform_accuracy}{Accuracy of the global isoform RF model, or
#'     \code{NA} when test labels are unavailable.}
#'   \item{global_gene_accuracy}{Accuracy of the global gene RF model, or
#'     \code{NA} when test labels are unavailable.}
#'   \item{global_overlapping_accuracy}{Accuracy on globally overlapping cells, or
#'     \code{NA} when test labels are unavailable.}
#'   \item{subset_overlapping_accuracy}{Accuracy on subset overlapping cells, or
#'     \code{NA} when test labels are unavailable.}
#'   \item{final_combined_accuracy}{Accuracy of the final hierarchical combined model, or
#'     \code{NA} when test labels are unavailable.}
#'   \item{confusion_matrix}{Confusion matrix for the final combined predictions, or
#'     \code{NULL} when test labels are unavailable.}
#'   \item{rf_top_features_file}{Path to the exported top-features workbook, or
#'     \code{NULL} if export was skipped.}
#' }
#'
#' @examples
#' \dontrun{
#' isoform_train <- readRDS("Iso_train.rds")
#' isoform_test <- readRDS("Iso_test.rds")
#' gene_train <- readRDS("gene_train.rds")
#' gene_test <- readRDS("gene_test.rds")
#'
#' result <- ICON(
#'   isoform_train   = isoform_train,
#'   isoform_test    = isoform_test,
#'   gene_train      = gene_train,
#'   gene_test       = gene_test,
#'   annotation_file = "Annotation_file.txt",
#'   barcode_col     = "sample",
#'   cell_type_col   = "stage",
#'   nfeatures       = 2000,
#'   top_n           = 100,
#'   nfeatures_model = 500,
#'   ntree           = 100,
#'   seed            = 241,
#'   assign_test_cell_types = TRUE,
#'   verbose         = TRUE
#' )
#'
#' head(result$final_predictions_df)
#' result$rf_top_features$Global_Isoform
#' }
#'
#' @export
ICON <- function(isoform_train,
                 isoform_test,
                 gene_train,
                 gene_test,
                 annotation_file,
                 barcode_col             = "barcode",
                 cell_type_col           = "cell_type",
                 nfeatures               = 2000,
                 top_n                   = 100,
                 nfeatures_model         = 500,
                 ntree                   = 100,
                 seed                    = 241,
                 export_predictions_file = "ICON_Final_Predictions.xlsx",
                 export_features_file    = "RF_Top_Features.xlsx",
                 assign_test_cell_types  = FALSE,
                 verbose                 = TRUE) {

  has_test_labels <- function(seurat_obj) {
    "cell_type" %in% colnames(seurat_obj@meta.data) &&
      any(!is.na(seurat_obj@meta.data$cell_type))
  }

  remove_existing_cell_type <- function(seurat_obj) {
    if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
      seurat_obj$cell_type <- NULL
    }
    seurat_obj
  }

  if (verbose)
    cat("\n=== Step 1: Preparing training datasets ===\n")

  isoform_train <- remove_existing_cell_type(isoform_train)
  gene_train <- remove_existing_cell_type(gene_train)

  if (assign_test_cell_types) {
    isoform_test <- remove_existing_cell_type(isoform_test)
    gene_test <- remove_existing_cell_type(gene_test)
  }

  isoform_train <- prepare_seurat_dataset(
    seurat_obj = isoform_train,
    annotation_file = annotation_file,
    barcode_col = barcode_col,
    cell_type_col = cell_type_col,
    nfeatures = nfeatures,
    assign_cell_types = TRUE,
    verbose = verbose
  )

  gene_train <- prepare_seurat_dataset(
    seurat_obj = gene_train,
    annotation_file = annotation_file,
    barcode_col = barcode_col,
    cell_type_col = cell_type_col,
    nfeatures = nfeatures,
    assign_cell_types = TRUE,
    verbose = verbose
  )

  if (verbose)
    cat("\n=== Step 2: Preparing testing datasets ===\n")

  isoform_test <- prepare_seurat_dataset(
    seurat_obj = isoform_test,
    annotation_file = annotation_file,
    barcode_col = barcode_col,
    cell_type_col = cell_type_col,
    nfeatures = nfeatures,
    assign_cell_types = assign_test_cell_types,
    verbose = verbose
  )

  gene_test <- prepare_seurat_dataset(
    seurat_obj = gene_test,
    annotation_file = annotation_file,
    barcode_col = barcode_col,
    cell_type_col = cell_type_col,
    nfeatures = nfeatures,
    assign_cell_types = assign_test_cell_types,
    verbose = verbose
  )

  if (verbose)
    cat("\n=== Step 3: Validating paired datasets ===\n")

  validate_paired_datasets(
    isoform_data = isoform_train,
    gene_data = gene_train,
    dataset_name = "training",
    verbose = verbose
  )

  validate_paired_datasets(
    isoform_data = isoform_test,
    gene_data = gene_test,
    dataset_name = "testing",
    verbose = verbose
  )

  evaluate_test_labels <- has_test_labels(isoform_test)

  global_isoform_accuracy <- NA_real_
  global_gene_accuracy <- NA_real_
  global_overlapping_accuracy <- NA_real_
  subset_overlapping_accuracy <- NA_real_
  final_combined_accuracy <- NA_real_
  confusion_matrix <- NULL

  RIU_result <- calculate_relative_feature_usage(
    isoform_train,
    top_n = top_n,
    verbose = verbose
  )
  RIU <- RIU_result$unique_top_features

  RGU_result <- calculate_relative_feature_usage(
    gene_train,
    top_n = top_n,
    verbose = verbose
  )
  RGU <- RGU_result$unique_top_features

  isoform_expr <- extract_expression_data(
    isoform_train,
    isoform_test,
    nfeatures = nfeatures_model,
    verbose = verbose
  )

  isoform_train_expr <- isoform_expr$train_expr
  isoform_test_expr <- isoform_expr$test_expr
  isoform_train_labels <- isoform_expr$train_labels

  gene_expr <- extract_expression_data(
    gene_train,
    gene_test,
    nfeatures = nfeatures_model,
    verbose = verbose
  )

  gene_train_expr <- gene_expr$train_expr
  gene_test_expr <- gene_expr$test_expr
  gene_train_labels <- gene_expr$train_labels

  if (verbose)
    cat("\n=== Step 4: Training global Random Forest models ===\n")

  global_iso_rf <- train_random_forest(
    isoform_train_expr,
    isoform_train_labels,
    ntree = ntree,
    seed = seed,
    verbose = verbose
  )
  global_iso_top_features <- get_top_rf_features(global_iso_rf$model)

  iso_preds <- predict_cell_types(
    global_iso_rf$model,
    isoform_test_expr,
    isoform_test,
    column_name = "Global_isoform_predicted_cell_type",
    verbose = verbose
  )
  isoform_test_data <- iso_preds$test_data

  if (evaluate_test_labels) {
    iso_eval <- evaluate_predictions(
      test_data = isoform_test_data,
      actual_col = "cell_type",
      predicted_col = "Global_isoform_predicted_cell_type",
      verbose = verbose
    )
    global_isoform_accuracy <- iso_eval$accuracy
    if (verbose) {
      cat("  Global isoform RF accuracy:",
          round(global_isoform_accuracy, 4), "\n")
    }
  } else if (verbose) {
    cat("  Test cell types not found. Skipping global isoform accuracy.\n")
  }

  global_gene_rf <- train_random_forest(
    gene_train_expr,
    gene_train_labels,
    ntree = ntree,
    seed = seed,
    verbose = verbose
  )
  global_gene_top_features <- get_top_rf_features(global_gene_rf$model)

  gene_preds <- predict_cell_types(
    global_gene_rf$model,
    gene_test_expr,
    gene_test,
    column_name = "Global_gene_predicted_cell_type",
    verbose = verbose
  )
  gene_test_data <- gene_preds$test_data

  if (evaluate_test_labels) {
    gene_eval <- evaluate_predictions(
      test_data = gene_test_data,
      actual_col = "cell_type",
      predicted_col = "Global_gene_predicted_cell_type",
      verbose = verbose
    )
    global_gene_accuracy <- gene_eval$accuracy
    if (verbose) {
      cat("  Global gene RF accuracy:",
          round(global_gene_accuracy, 4), "\n")
    }
  } else if (verbose) {
    cat("  Test cell types not found. Skipping global gene accuracy.\n")
  }

  if (verbose)
    cat("\n=== Step 5: Comparing global predictions ===\n")

  global_overlap_result <- identify_overlapping_cells(
    isoform_test_data,
    gene_test_data,
    isoform_pred_col = "Global_isoform_predicted_cell_type",
    gene_pred_col = "Global_gene_predicted_cell_type",
    predicted_global_col = "predicted_global",
    verbose = verbose
  )
  global_overlap <- global_overlap_result$overlap_data
  global_non_overlap_isoforms <- global_overlap_result$non_overlap_isoforms
  global_non_overlap_genes <- global_overlap_result$non_overlap_genes

  if (evaluate_test_labels && !is.null(global_overlap)) {
    overlap_eval <- evaluate_overlapping_predictions(
      global_overlap,
      predicted_col = "predicted_global",
      verbose = verbose
    )
    global_overlapping_accuracy <- overlap_eval$accuracy
    if (verbose) {
      cat("  Global overlapping cells accuracy:",
          round(global_overlapping_accuracy, 4), "\n")
    }
  } else if (verbose) {
    cat("  Skipping global overlapping accuracy calculation.\n")
  }

  if (verbose)
    cat("\n=== Step 6: Preparing subset model data ===\n")

  isoform_subset_expr <- extract_expression_data(
    isoform_train,
    global_non_overlap_isoforms,
    features = RIU,
    nfeatures = nfeatures_model,
    verbose = verbose
  )
  isoform_subset_train_expr <- isoform_subset_expr$train_expr
  isoform_subset_test_expr <- isoform_subset_expr$test_expr
  isoform_subset_train_labels <- isoform_subset_expr$train_labels

  gene_subset_expr <- extract_expression_data(
    gene_train,
    global_non_overlap_genes,
    features = RGU,
    nfeatures = nfeatures_model,
    verbose = verbose
  )
  gene_subset_train_expr <- gene_subset_expr$train_expr
  gene_subset_test_expr <- gene_subset_expr$test_expr
  gene_subset_train_labels <- gene_subset_expr$train_labels

  if (verbose)
    cat("\n=== Step 7: Training subset Random Forest models ===\n")

  subset_iso_top_features <- NULL
  subset_gene_top_features <- NULL

  subset_iso_rf <- train_random_forest(
    isoform_subset_train_expr,
    isoform_subset_train_labels,
    ntree = ntree,
    seed = seed,
    verbose = verbose
  )
  if (!is.null(subset_iso_rf)) {
    subset_iso_top_features <- get_top_rf_features(subset_iso_rf$model)
  }

  iso_sub_preds <- predict_cell_types(
    subset_iso_rf$model,
    isoform_subset_test_expr,
    global_non_overlap_isoforms,
    column_name = "Subset_isoform_predicted_cell_type",
    verbose = verbose
  )
  isoform_subset_test_data <- iso_sub_preds$test_data

  subset_gene_rf <- train_random_forest(
    gene_subset_train_expr,
    gene_subset_train_labels,
    ntree = ntree,
    seed = seed,
    verbose = verbose
  )
  if (!is.null(subset_gene_rf)) {
    subset_gene_top_features <- get_top_rf_features(subset_gene_rf$model)
  }

  gene_sub_preds <- predict_cell_types(
    subset_gene_rf$model,
    gene_subset_test_expr,
    global_non_overlap_genes,
    column_name = "Subset_gene_predicted_cell_type",
    verbose = verbose
  )
  gene_subset_test_data <- gene_sub_preds$test_data

  if (verbose)
    cat("\n=== Step 8: Comparing subset predictions ===\n")

  subset_overlap_result <- identify_overlapping_cells(
    isoform_subset_test_data,
    gene_subset_test_data,
    isoform_pred_col = "Subset_isoform_predicted_cell_type",
    gene_pred_col = "Subset_gene_predicted_cell_type",
    predicted_global_col = "predicted_subset",
    verbose = verbose
  )
  subset_overlap <- subset_overlap_result$overlap_data
  subset_non_overlap_isoforms <- subset_overlap_result$non_overlap_isoforms

  if (evaluate_test_labels && !is.null(subset_overlap)) {
    subset_overlap_eval <- evaluate_overlapping_predictions(
      subset_overlap,
      predicted_col = "predicted_subset",
      verbose = verbose
    )
    subset_overlapping_accuracy <- subset_overlap_eval$accuracy
    if (verbose) {
      cat("  Subset overlapping cells accuracy:",
          round(subset_overlapping_accuracy, 4), "\n")
    }
  } else if (verbose) {
    cat("  Skipping subset overlapping accuracy calculation.\n")
  }

  if (verbose)
    cat("\n=== Step 9: Combining final predictions ===\n")

  final_eval <- combined_final_predictions(
    global_overlap,
    subset_overlap,
    subset_non_overlap_isoforms,
    global_pred_col = "predicted_global",
    subset_pred_col = "predicted_subset",
    subset_non_overlap_pred_col = "Subset_isoform_predicted_cell_type",
    export_file = export_predictions_file,
    verbose = verbose
  )

  if (evaluate_test_labels) {
    final_combined_accuracy <- final_eval$accuracy_manual
    confusion_matrix <- final_eval$confusion_matrix
  }

  rf_top_features <- list(
    Global_Isoform = global_iso_top_features,
    Global_Gene = global_gene_top_features
  )

  if (!is.null(subset_iso_top_features)) {
    rf_top_features$Subset_Isoform <- subset_iso_top_features
  }

  if (!is.null(subset_gene_top_features)) {
    rf_top_features$Subset_Gene <- subset_gene_top_features
  }

  rf_top_features_file <- NULL
  if (!is.null(export_features_file)) {
    rf_top_features_file <- export_rf_top_features(
      top_feature_list = rf_top_features,
      file = export_features_file,
      verbose = verbose
    )
  }

  if (verbose && evaluate_test_labels) {
    cat("\n============================================================\n")
    cat("                    ICON Results Summary                    \n")
    cat("============================================================\n")
    cat("  Global isoform RF accuracy:  ", round(global_isoform_accuracy, 4), "\n")
    cat("  Global gene RF accuracy:     ", round(global_gene_accuracy, 4), "\n")
    cat("  Final combined accuracy:     ", round(final_combined_accuracy, 4), "\n")
    cat("============================================================\n")
  }

  list(
    final_predictions_df = final_eval$final_predictions_df,
    rf_top_features = rf_top_features,
    global_isoform_accuracy = global_isoform_accuracy,
    global_gene_accuracy = global_gene_accuracy,
    global_overlapping_accuracy = global_overlapping_accuracy,
    subset_overlapping_accuracy = subset_overlapping_accuracy,
    final_combined_accuracy = final_combined_accuracy,
    confusion_matrix = confusion_matrix,
    rf_top_features_file = rf_top_features_file
  )
}
