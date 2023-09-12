everyhting.plot <- function(df, 
                            grp_col, 
                            base_name, 
                            plot_cols, 
                            color_palette, 
                            additional_layers=NULL, 
                            plot.type = "pairs",
                            col.panel = "gray10",
                            shape = 1,size = 2) {
  
  # Check if required libraries are installed
  required_libraries <- c("dplyr", "GGally", "ggplot2", "crayon","plot3D")
  new_libraries <- required_libraries[!(required_libraries %in% installed.packages()[,"Package"])]
  if(length(new_libraries)) stop(cat(crayon::red("Missing packages: "), crayon::blue(toString(new_libraries)), ". Please install them."))
  
  # Filter numeric columns and update user if any column is dropped
  numeric_cols <- select_if(df, is.numeric)
  dropped_cols <- setdiff(names(df), names(numeric_cols))
  if (length(dropped_cols) > 0) {
    cat(crayon::yellow("Dropped non-numeric columns: "), crayon::blue(toString(dropped_cols)), "\n")
  }
  
  # Select only the specified number of columns for plotting
  if (ncol(numeric_cols) < plot_cols) {
    stop(cat(crayon::red("Error: "), crayon::blue("The number of numeric columns in the dataframe is fewer than 'plot_cols'.")))
  }
  selected_cols <- numeric_cols[, 1:plot_cols]
  
  # Set column names
  col_names <- paste0(base_name, seq_along(1:plot_cols))
  names(selected_cols) <- col_names
  
  print(colnames(selected_cols))
  
  
  if (plot.type == "pairs") {
    # Create the GGally plot
    plot <- selected_cols %>%
      dplyr::mutate(type = grp_col) %>%
      GGally::ggpairs(., columns = 1:plot_cols, progress = FALSE,
                      ggplot2::aes(color = type, fill = type),
                      lower = list(continuous = wrap("points", alpha = 0.95, size = 2, pch=21)),
                      diag = list(continuous = wrap("densityDiag", alpha = 0.85, color= NA)),
                      upper = list(continuous = "blank"), legend = 1 ) +
      scale_color_manual(values= color_palette) +
      scale_fill_manual(values= color_palette) +
      theme_test() +
      theme(legend.position = "top")
    
    # Check if grp_col exists in data frame
    cat(crayon::blue("The grouping factor has the following levels:"), crayon::red(paste0(unique(grp_col), collapse = " | ")),"\n")
    
    # Add additional ggplot layers if provided
    if (!is.null(additional_layers)) {
      plot <- plot + additional_layers
    }
    return(plot)
  }
  
  if (plot.type == "3d") {
    plot2 <- plot3D::scatter3D(x = selected_cols[,1], z = selected_cols[,2], y = selected_cols[,3],
                               col = color_palette, colvar = as.numeric(factor(grp_col)), colkey = list(plot = F),
                               pch = shape ,
                               cex = 0.85,
                               theta = 105, phi = 48, lwd.axis= 0.8, bty = "u",
                               expand = 0.2, 
                               r=6, d = 0.2, 
                               col.panel = col.panel, 
                               col.grid = "gray40", 
                               ticktype = "detailed", 
                               nticks = 4,
                               lwd.grid = 0.1, 
                               lwd.panel = 1.5, 
                               col.axis ="white", lwd = 1.5,
                               xlab = paste0(base_name, seq.int(plot_cols)[1]), 
                               zlab = paste0(base_name, seq.int(plot_cols)[2]), 
                               ylab = paste0(base_name, seq.int(plot_cols)[3]), 
                               main = "")
    
    # Check if grp_col exists in data frame
    cat(crayon::blue("The grouping factor has the following levels:"), crayon::red(paste0(unique(grp_col), collapse = " | ")),"\n")
  }
  
  if (plot.type == "simple") {
    
    # Calculate the minimum and maximum x-values from the data
    min_x <- round(min(selected_cols[, 1]))
    max_x <- round(max(selected_cols[, 1]))
    
    # Add or subtract 10 from min and max, keeping the direction same as original value
    min_x_adj <- if (min_x >= 0) min_x + 10 else min_x - 10
    max_x_adj <- if (max_x >= 0) max_x + 10 else max_x - 10
    
    # Generate breaks and labels
    breaks <- seq(min_x_adj, max_x_adj, by = 10)
    limits <- ifelse(breaks >= 0, breaks + 5, breaks - 5)
    
    # Create the plot
    plot3 <- selected_cols %>%
      dplyr::mutate(type = grp_col) %>%
      ggpubr::ggscatter(
        x = paste0(base_name, seq.int(plot_cols)[1]), 
        y = paste0(base_name, seq.int(plot_cols)[2]),
        color = "type", 
        palette = col.mat_2,
        size = size, 
        shape = shape
      ) + 
      theme_classic2() + 
      labs(color = paste0("n = ", nrow(selected_cols))) +
      scale_x_continuous(
        limits = c(min_x_adj, max_x_adj), 
        breaks = breaks
      )
    
    # Check if grp_col exists in data frame
    cat(crayon::blue("The grouping factor has the following levels:"), 
        crayon::red(paste0(unique(grp_col), collapse = " | ")),"\n")
    # Add additional ggplot layers if provided
    if (!is.null(additional_layers)) {
      plot3 <- plot3 + additional_layers
    }
    
    return(plot3)
  }
  
}
