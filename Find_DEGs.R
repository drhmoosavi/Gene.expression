Find_DEGs <- function(expression_matrix, 
                      class, 
                      lfc = log2(1.5),
                      padj= 0.05, 
                      comparisonType="1vRest", 
                      sortBy="adj.P.Val",
                      convert.ids = TRUE,
                      volcano.plot = TRUE) {
  
  if (is.matrix(expression_matrix) == FALSE) {
    cat(crayon::yellow("Data is not a matrix. Coercing to a matrix.\n\n"))
    expression_matrix <- as.matrix(expression_matrix)
  } else {
    cat(crayon::green("Data isa matrix.\n\n"))
  }
  
  # Check for required libraries
  if (!requireNamespace("limma", quietly=TRUE)) {
    stop(crayon::red("Function requires the limma package available from Bioconductor.\n\n"))
  }
  
  # Check for valid inputs
  if (length(class) != ncol(expression_matrix)) {
    stop(crayon::red("Length of class must be equal to the number of columns in the expression matrix.\n\n"))
  }
  
  # Data Preprocessing
  cat(crayon::green("Preprocessing data...\n\n"))
  
  # Ensure class is a factor
  class <- as.factor(class)
  
  # Drop samples with unknown class label
  keepN <- !is.na(class)
  class <- droplevels(class[keepN])
  expression_matrix <- expression_matrix[, keepN]
  
  # Differential Expression Analysis
  cat(crayon::green("Performing differential expression analysis...\n\n"))
  
  # Design matrix creation
  design <- model.matrix(~0+class)
  colnames(design) <- levels(class)
  
  # Linear modeling with limma
  fit <- lmFit(expression_matrix, design, method="ls")
  cat(crayon::green("lmFit done! Fitting method; `least squares`.\n\n"))
  
  # Empirical Bayes statistics
  fit <- eBayes(fit, robust= TRUE)
  cat(crayon::green("eBayes done!\n\n"))
  
  # Perform 1vRest or pairwise comparisons
  if (comparisonType == "1vRest") {
    # each group is individually compared to the average of all other groups
    num_levels <- length(levels(class))
    contrast_strings <- sapply(levels(class), function(level) {
      other_levels <- setdiff(levels(class), level)
      contrast_formula <- paste(c(rep(1, length(level)), rep(-1/length(other_levels), length(other_levels))), 
                                c(level, other_levels), sep = "*", collapse =" + ")
      return(contrast_formula)
    })
    
    contrasts <- makeContrasts(contrasts = contrast_strings, levels = design)
    colnames(contrasts) <- levels(class)
    rownames(contrasts) <- levels(class)
    
    cat(crayon::green(paste0("Performing 1 versus rest comparisons for ", num_levels, " levels.\n\n")))
    
  } else if (comparisonType == "pairwise") { 
    
    contrasts <- makeContrasts(
      contrasts = combn(levels(class), 2, function(x) paste(x, collapse=" - ")),
      levels = design
    )
    cat(crayon::green(paste0("Performing pairwise comparisons for ", length(levels(class)), " levels, resulting in ", choose(length(levels(class)), 2), " comparisons.\\n")))
  } else {
    
    stop("Invalid comparisonType. Valid options are '1vRest', or 'pairwise'.")
    
  }
  print(contrasts)
  
  fit2 <- contrasts.fit(fit, contrasts)
  cat(crayon::green("Estimated coefficients were computed using `contrasts.fit()` function!\n\n"))
  fit2 <- eBayes(fit2, robust= TRUE)
  cat(crayon::green("`eBayes()` with `robust= TRUE` done!\n\n"))
  
  # Initialize a list to store results for all contrasts
  results_list <- list()
  
  # Iterate through the contrasts and apply topTable to each one
  for (i in seq_len(ncol(contrasts))) {
    contrast_name <- colnames(contrasts)[i]
    results <- topTable(fit2, adjust.method= "BH", number= Inf, coef= i, lfc = 0, p.value = padj , sort.by = "p")
    
    # Converting to a data frame if necessary
    if (is.matrix(results)) {
      results <- as.data.frame(results)
    }
    
    results_list[[contrast_name]] <- results
    
  }
  
  # Convert Gene Ids to Gene SYMBOL
  #if(all(sapply(x, function(df) all(!is.na(as.numeric(row.names(df)))))) | !missing(convert.ids) && convert.ids == TRUE )
  if(!missing(convert.ids) && convert.ids == TRUE ) { 
    
    # Sorting by adj.P.Val or logFC, as specified by the user
    results <- results %>% dplyr::arrange(!!rlang::sym(sortBy))
    
    if (!requireNamespace("org.Hs.eg.db", quietly=TRUE)) {
      stop("The package org.Hs.eg.db is required for this function to work. Please install it. OR
            set `convert.ids = FALSE`")
    }
    
    cat(bold(crayon::yellow("Convert Gene IDs to SYMBOLs using `org.Hs.eg.db` package!\n\n")))
    results_list <- purrr::map(results_list, function(df) {
      mutate(df, gene = unname(annotate::getSYMBOL(row.names(df), "org.Hs.eg.db")))
    })
  }
  
  results_list_volcano <- results_list
  
  results_list <- purrr::map(results_list, ~dplyr::filter(.x, logFC > lfc | logFC < -lfc))
  results_list <- purrr::map(results_list, ~dplyr::arrange(.x, -adj.P.Val))
  
  # Print a message indicating the contrasts and results
  cat(crayon::green(paste0("Results obtained for ", bold(red(length(results_list))), " contrasts.\n\n")))
  
  cat(bold(yellow("Number of DEGs for " %+% red(paste(levels(class), collapse =" & ")) %+% " in comparison\n")))
  
  sum_fit2 <- decideTests(fit2, p.value = padj, adjust = "fdr", lfc=lfc)
  sum_fit2_summary <- as.data.frame(summary(sum_fit2))
  print(sum_fit2)
  sum_fit2_summary <- stats::setNames(sum_fit2_summary, c("Direction","Comparison","Number"))
  
  print(knitr::kable(sum_fit2_summary, "simple"))
  cat("\n")
  
  
  # Check if volcano.plot is TRUE and then produce the plots
  if(volcano.plot == TRUE) {
    
    # Check for required libraries
    if (!requireNamespace("EnhancedVolcano", quietly=TRUE)) {
      stop("Function requires the EnhancedVolcano package.")
    }
    if (!requireNamespace("cowplot", quietly=TRUE)) {
      stop("Function requires the cowplot package.")
    }
    cat(bold(yellow("Making `volcano` plots: " %+% red(paste(length(results_list_volcano))),"\n")))
    # Initialize an empty list to store ggplot objects
    volcano_plot_list <- list()
    
    # Loop through results_list list and create volcano plots
    for (name in names(results_list_volcano)) {
      
      # Extract individual data frame
      df <- results_list_volcano[[name]]
      
      # Create keyvals.colour specific to this data frame
      keyvals.colour <- dplyr::case_when(
        df$logFC < -1 ~ "darkgoldenrod1",  # Adapt the cutoffs as needed
        df$logFC > 1 ~ "red3",
        is.na(df$logFC) ~ "lightblue3",
        TRUE ~ "lightblue3"
      ) %>% stats::setNames(dplyr::case_when(
        . == "red3" ~ "Up",
        . == "lightblue3" ~ "Not DE",
        . == "darkgoldenrod1" ~ "Down"
      ))
      
      # Create a volcano plot
      suppressMessages(
        plot_volcano <- EnhancedVolcano::EnhancedVolcano(toptable = df, x = "logFC", y = "adj.P.Val", 
                                                         lab = df$gene, 
                                                         pCutoff = 0.0001,
                                                         FCcutoff=1, 
                                                         gridlines.major = T, 
                                                         gridlines.minor = F ,
                                                         axisLabSize=8,
                                                         labCol="pink4", 
                                                         labSize=3,
                                                         labFace=4,
                                                         boxedLabels = F, 
                                                         drawConnectors= F, 
                                                         widthConnectors = 0.3,
                                                         colConnectors = 'lightblue2',
                                                         colAlpha = 0.95,
                                                         xlab = "",
                                                         ylab = "",
                                                         title= "",
                                                         subtitle="",
                                                         legendPosition = 'none',
                                                         legendLabels = "",
                                                         legendIconSize = 3,
                                                         legendLabSize= 7,
                                                         colCustom = keyvals.colour,
                                                         shape = 20,
                                                         pointSize = 2.2,  
                                                         vlineWidth = 1) +  
          theme_pubclean() +
          theme(legend.position = "top", 
                legend.direction = "horizontal",
                legend.box = "horizontal",
                legend.justification = "center",
                legend.key = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_text(size= 9, face = 2, color = "navy"),
                axis.text.y = element_text(vjust = 0.3),
                axis.line.x = element_line(color= "lightskyblue4", linewidth = 1.5 ),
                panel.grid.major.y = element_line(color = "snow3", linewidth = 0.3, linetype = 1)) +
          scale_x_continuous(limits = c(ceiling(abs(min(df$logFC)))*-1,
                                        ceiling(max(df$logFC) ))) +
          # scale_y_continuous(breaks = c(50,100),
          #                    limits = c(0, ceiling(max(abs(log10(df$adj.P.Val)))) )) +
          labs(color="", ylab = "") + 
          coord_equal(ratio = 0.02)
      )
      
      # Store the ggplot object in the list
      volcano_plot_list[[name]] <- plot_volcano
    }
    
    # Determine the number of columns for grid
    num_plots <- length(volcano_plot_list)
    num_cols <- ceiling(sqrt(num_plots))
    
    # Generate labels based on the names of the data frames
    labels <- LETTERS[1:length(names(volcano_plot_list))]
    
    volcano_plot_grid <- plot_grid(
      plotlist = volcano_plot_list,
      ncol = num_cols,
      labels = labels,
      align = "hv"
    )
    
    volcano_plot_grid <- 
      annotate_figure(volcano_plot_grid,
                      bottom = text_grob(bquote(~italic(Log[2])~ "Fold Change"), 
                                         color = "black", hjust = 1, face = 3, size = 10),
                      left = text_grob("adj. p-val", color = "black", 
                                       rot = 90, face = 3, size = 10))
    # Show the final plot grid
    # print(volcano_plot_grid)
  }
  
  # Return the list of results
  return(list("DEGs" = results_list, "decidetest" = sum_fit2, "volcano_plots" = volcano_plot_grid))
  
}
