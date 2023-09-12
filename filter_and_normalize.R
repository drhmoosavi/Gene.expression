filter_and_normalize <- function(expression_set, 
                                 subset = NULL) {
  
  # Ensure the required libraries are loaded
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("`Biobase` package must be installed and loaded.")
  }
  
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("`crayon` package must be installed and loaded.")
  }
  
  expression_matrix <- Biobase::exprs(expression_set)
  
  # Remove "_at" suffix from row names
  rownames(expression_matrix) <- gsub("\\_at$", "", rownames(expression_matrix))
  
  # Check if dimension of expression_matrix is zero
  if (all(dim(expression_matrix) == 0)) {
    stop(crayon::red("Dimension of expression matrix is zero. Please provide a valid expression matrix."))
  }
  
  # Check if subset is provided
  if (!is.null(subset)) {
    cat(bold(crayon::green("Length of subset vector:")), paste(yellow(length(subset))) ,"\n")
    
    # Filter rows based on subset (e.g., Entrez IDs)
    expression_matrix <- expression_matrix[rownames(expression_matrix) %in% subset,]
    cat("Dimensions after filtering:", paste(yellow(dim(expression_matrix))))
    
    # Check intersection between row names and subset
    intersection <- intersect(rownames(expression_matrix), subset)
    cat(bold(green(" Intersection between row names and the subset")))
  } else {
    cat("Subset not provided; proceeding without filtering.")
  }
  
  # Ensure the required library is loaded
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("`limma` package must be installed and loaded.")
  }
  
  # Normalize between arrays using limma package
  cat(crayon::bold("\n`Limma` Quantile Normalizating between arrays"), "\n")
  expression_matrix <- limma::normalizeBetweenArrays(expression_matrix, method = "quantile")
  
  return(expression_matrix)
}
