bg_adjustment <- function(expression_data, 
                          gene_sets_list,
                          subetN = NULL, 
                          use_subset = TRUE,
                          type.lm = "simple") {
  
  # Initialize a matrix to store Pi values for each gene set
  if (use_subset == TRUE && !is.null(subetN)) {
    col_orders <- colnames(expression_data)
    subset_expression_data <- expression_data[, subetN]
    Pi_matrix <- matrix(0, nrow = length(gene_sets_list), ncol = ncol(subset_expression_data))
  } else {
    Pi_matrix <- matrix(0, nrow = length(gene_sets_list), ncol = ncol(expression_data))
  }
  
  for (i in seq_along(gene_sets_list)) {
    # print(paste("Iterating for gene set:", i))

    if (use_subset == TRUE && !is.null(subetN)) {
      row_centered_mat <- sweep(subset_expression_data, 1, rowMeans(subset_expression_data, na.rm = TRUE))
    # row_centered_mat <- scale(subset_expression_data, center = rowMeans(subset_expression_data, na.rm = TRUE), scale = FALSE)
    } else {
      row_centered_mat <- sweep(expression_data, 1, rowMeans(expression_data, na.rm = TRUE))
    }
    
    Pi <- row_centered_mat[rownames(row_centered_mat) %in% gene_sets_list[[i]], ] %>% 
      colSums(na.rm = TRUE) %>% {(.-min(.)) / (max(.)-min(.))}
      # range01()
    
    if(length(Pi) == ncol(Pi_matrix)) {
      Pi_matrix[i, ] <- Pi
    } else {
      stop("Length of Pi does not match the number of columns in Pi_matrix.")
    }
  }
  
  # Convert Pi_matrix to a data frame
  Pi_df <- as.data.frame(t(Pi_matrix))
  cat("Number of `Independent Variable`: ", nrow(t(Pi_df)),"\n")
  
  # Create formula function
  make_lm_formula <- function(Invar, dependent_var = "y", type = type.lm) {
    if (type == "simple") {
      return(formula(paste(dependent_var, "~", Invar[1])))
    } else if (type == "multiple") {
      return(formula(paste(dependent_var, "~", paste(Invar, collapse = " + "))))
    } else {
      stop("Invalid type. Use 'simple' or 'multiple'.")
    }
  }
  
  # Fit a linear model with multiple predictors
  if (use_subset == TRUE && !is.null(subetN)) {
  cat("Dimention of expression matrix(response variable): ",dim(t(subset_expression_data)))
    
    # Create lm formula using make_lm_formula function
    lm_formula <- make_lm_formula(Invar, dependent_var = "subset_expression_data" , type = type.lm)
    
    fit <- lm(t(subset_expression_data) ~ ., data = Pi_df)
  } else {
    cat("Dimention of expression matrix:\n",print(dim(t(expression_data))))
    
    # Create lm formula using make_lm_formula function
    lm_formula <- make_lm_formula(Invar, dependent_var = "expression_data" , type = type.lm)
    
    fit <- lm(lm_formula, data = Pi_df)
  }
  
  Adjusted_data <- t(residuals(fit))
  
  # Apply lnn function to Adjusted gene expression data
  lnn <- function(x) x + abs(min(x))
  row_centered_mat <- lnn(Adjusted_data)
  
  if (use_subset == TRUE && !is.null(subetN)) {
    combined_result <- BiocGenerics::combine(subset_expression_data, expression_data[,-subetN])
    # print(dim(combined_result))
    combined_result <- combined_result[, match(col_orders, colnames(combined_result))]
    return(invisible(combined_result))
  } else {
    return(invisible(row_centered_mat))
  }
  
}

# bg_adjustment(exp_liv.pri, 
#               gene_sets_list = all_gene_sets[grepl("MacP_",names(all_gene_sets))],
#               subetN = grepl("met",info_df_p2$type2),
#               type.lm = "simple")
