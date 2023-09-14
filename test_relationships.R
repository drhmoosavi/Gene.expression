# Define a function to perform and evaluate a series of regression analyses
test_relationships <- function(df, target_prefix = "^PC", ignore.case = F) {
  # Initialization and setup
  all_results <- data.frame()
  target_columns <- grep(target_prefix, names(df), value = TRUE)
  
  # Loop through target and predictor columns
  for (target in target_columns) {
    for (predictor in setdiff(names(df), target)) {
      
      # Model Fitting
      models <- list(
        lm(as.formula(paste(target, "~", predictor)), data = df),
        lm(as.formula(paste(target, "~", predictor, "+ I(", predictor, "^2)")), data = df),
        lm(as.formula(paste(target, "~", predictor, "+ I(", predictor, "^2) + I(", predictor, "^3)")), data = df)
      )
      
      # Function for extracting necessary summary statistics from a model
      extract_summary <- function(model, coef_index) {
        summary_model <- summary(model)
        coef_values <- summary_model$coefficients[coef_index, ]
        return(data.frame(Estimate = coef_values[1],
                          p_value = coef_values[4],
                          R2 = summary_model$r.squared))
      }
      
      # Extract statistics using `extract_summary` function
      extracted_results <- lapply(models, extract_summary, coef_index = 2)
      
      # Combine results and append
      row_results <- data.frame(Target = target, Predictor = predictor, 
                                do.call(cbind, extracted_results), 
                                stringsAsFactors = FALSE)
      all_results <- rbind(all_results, row_results)
    }
  }
  return(all_results)
}


# Call the function
# results <- data.frame(cbind(t(gsva.res_2), dim.r2$pca[,1:3]), check.names = T) %>% 
#   test_relationships(.)
# 
# results %>% 
#   group_by(Target) %>% 
#   slice_max(n = 5, order_by = Linear_R2 , with_ties = F) %>% 
#   print(n= 25) %>% generate_kable
