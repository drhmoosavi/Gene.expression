limma_camera_test <- function(y, 
                              class, 
                              gene_set_collection,
                              minGeneSetSize = 5,
                              maxGeneSetSize = 1000,
                              use.ranks = T, 
                              inter.gene.cor = 0.01, 
                              trend.var = FALSE, 
                              sort = TRUE,
                              comparisonType = "1vRest") {
  
  # Validate input length
  if (length(class) != ncol(y)) {
    stop("Length of class must match number of columns in y")
  }
  
  cat(yellow$bold(paste0("Min gene set size is set to: ", minGeneSetSize, "\n")))
  cat(yellow$bold(paste0("Max gene set size is set to: ", maxGeneSetSize, "\n")))
  
  # Ensure class is a factor and drop samples with unknown class label
  class <- as.factor(class)
  keepN <- !is.na(class)
  class <- droplevels(class[keepN])
  y <- y[, keepN]
  
  # Create index for CAMERA from gene_set_collection
  index <- lapply(names(gene_set_collection), function(n) {
    which(rownames(y) %in% gene_set_collection[[n]])
  })
  names(index) <- names(gene_set_collection)
  index <- index[sapply(index, function(genes) {
    length(genes) >= minGeneSetSize & length(genes) <= maxGeneSetSize
  })]
  
  # Create design matrix
  design <- model.matrix(~ 0 + class, data = class)
  colnames(design) <- levels(class)
  
  # Initialize results list
  camera_results_list <- list()
  
  # Create contrasts and run CAMERA for each
  if (comparisonType == "1vRest") {
    cat(green$bold("Building contrast strings...\n\n"))
    
    contrast_strings <- sapply(levels(class), function(this_level) {
      cat(yellow$bold(paste0("Processing ", this_level, "...\n\n")))
      
      other_levels <- setdiff(levels(class), this_level)
      cat(magenta$bold(paste0("Other levels: ", paste(other_levels, collapse = ", "))),"\n")
      
      n_other_levels <- length(other_levels)
      
      other_part <- paste(other_levels, collapse = " + ")
      cat(magenta$bold(paste0("Other part: ", other_part)),"\n\n")
      
      contrast_str <- paste0(this_level, " - (", other_part, ")/", n_other_levels)
      cat(yellow$bold(paste0("Contrast string for ", this_level, ": ", contrast_str)),"\n")
      
      return(contrast_str)
    })
    
    cat(green$bold("Generated contrast strings:"),"\n")
    print(contrast_strings)
    
    contrast.matrix <- makeContrasts(contrasts = contrast_strings, levels = design)
    colnames(contrast.matrix) <- names(contrast_strings)
    
    cat(green$bold("\nGenerated contrast matrix:\n"))
    # Print the matrix itself. Assuming it doesn't need coloring.
    print(contrast.matrix)
    
    for (level in levels(class)) {
      if (!level %in% colnames(contrast.matrix)) {
        stop(paste0("No matching contrast found for level: ", level))
      }
      res <- camera(y = y, index = index, 
                    design = design,
                    contrast = contrast.matrix[, level],
                    use.ranks = use.ranks, 
                    inter.gene.cor = inter.gene.cor, 
                    trend.var = trend.var, 
                    sort = sort)
      camera_results_list[[level]] <- res
    }
  }
  
  return(camera_results_list)
}
