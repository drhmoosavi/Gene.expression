select_top_genes <- function(gene_expression_matrix, 
                             method = "sd", percentage = NULL,
                             top_n = NULL) {
  if (!is.matrix(gene_expression_matrix)) {
    stop("Input must be a matrix.")
  }
  
  if (is.null(percentage) && is.null(top_n)) {
    stop("Either percentage or top_n must be provided.")
  }
  
  if (!is.null(top_n) && (top_n <= 0 || top_n > nrow(gene_expression_matrix))) {
    stop("Invalid value for top_n. It must be between 1 and the number of genes.")
  }
  
  if (!is.null(percentage) && (percentage <= 0 || percentage > 1)) {
    stop("Invalid value for percentage. It must be between 0 and 1.")
  }
  
  # Calculate the chosen variation for each gene
  variation <- NULL
  if (method == "sd") {
    variation <- rowSds(gene_expression_matrix)
  } else if (method == "mad") {
    variation <- rowMads(gene_expression_matrix)
  } else if (method == "iqr") {
    variation <- rowIQRs(gene_expression_matrix)
  } else {
    stop("Invalid method. Choose 'sd', 'mad', or 'iqr'.")
  }
  
  # Determine the threshold based on percentage or top_n
  threshold <- NULL
  if (!is.null(top_n)) {
    threshold <- sort(variation, decreasing = TRUE)[top_n]
  } else {
    threshold <- quantile(variation, 1 - percentage)
  }
  
  # Select the rows (genes) that are above or equal to the threshold
  selected_genes <- gene_expression_matrix[variation >= threshold,]
  
  return(selected_genes)
}
