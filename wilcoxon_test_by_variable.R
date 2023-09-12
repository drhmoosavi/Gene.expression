wilcoxon_test_by_variable <- function(df,
                                      type = NULL, 
                                      p_value = NULL) {
  
  type <- factor(type)
  
  df <- data.frame(t(df),
                   type = type ) %>% 
    gather(key = "variable", value = "value", -c(type))
  
  print(head(df))
  
  # Check input
  if (! all(c("variable", "value", "type") %in% colnames(df)))
    stop(crayon::red("Required columns not found in data frame"))
  
  
  # Split the data by the variable factor and apply the Wilcoxon test to each subset
  cat(crayon::green("Apply Wilcoxon test by the variable `type` as factor...\n\n"))
  
  results <- df %>%
    group_by_at("variable") %>%
    do({
      if (length(unique(.$type)) != 2) 
        stop(crayon::red("Two groups are required for each variable subset"))
      
      group_means <- tapply(.$value, .$type, mean)
      fold_change <- group_means[2] / group_means[1]
      higher_median_group <- names(which.max(tapply(.$value, .$type, median)))
      
      # Wilcoxon test
      test_result <- wilcox.test(as.formula(paste("value", "~", "type")), data = .)
      
      data.frame(variable = unique(.$variable), 
                 p_value = test_result$p.value, 
                 higher_median_group = higher_median_group, 
                 fold_change = fold_change)
    })
  
  cat(crayon::green("The test is done!\n\n"))
  
  cat(crayon::green("Correction for multiple testing using False Discovery Rate (FDR)...\n\n"))
  # Correct the p-values using the False Discovery Rate (FDR) method
  results$p_value <- stats::p.adjust(results$p_value, method = "fdr")
  
  cat(crayon::green("Filter results according to the cut of the ") %+% crayon::red(italic("p-values!\n\n")) )
  results <- results %>% 
    filter(p_value <= p_value ) %>%
    arrange(p_value)
  
  return(results)
}
