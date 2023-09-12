GSEA.plot <- function(mylist,
                      n.top = 5,
                      x_var = "name",
                      xlab = "",
                      y_var = "logFDR",
                      ylab = "",
                      fill_var = "logFDR",
                      sort_val = "desc",
                      sort_by_groups = TRUE,
                      ncol_grid = 2,
                      ggtheme = theme_cowplot(),
                      lgd.pos = "top",
                      palette.col = "Blue-Red 3",
                      additional_ggplot_layers = NULL,
                      rescale = FALSE,
                      plot.type = "bar",
                      common.legend = TRUE) {
  if (plot.type == "heat") {
    all_df <-
      bind_rows(GSEA.table(mylist, n = n.top), .id = 'DataFrame')
    
    searchdf <-
      bind_rows(GSEA.table(mylist, n = Inf, fdr_threshold = 1), .id = 'DataFrame')
    
    all_df <-
      searchdf[searchdf$name %in% unique(as.character(all_df$name)), ]
    all_df$name <- factor(all_df$name, levels = unique(all_df$name))
    
    # Assuming you're using the first data frame in mylist
    df <- all_df
    
    if (rescale) {
      df <- df %>%
        group_by(name) %>%
        mutate(logFDR = scales::rescale(logFDR, to = c(-1, 1))) %>%
        ungroup()
    }
    
    heatmap_plot <-
      ggplot(df, aes(x = DataFrame, y = name, fill = logFDR)) +
      geom_tile(linewidth = 0.6, color = "white") +
      scale_fill_continuous_diverging(palette = palette.col) +
      theme_bw() +
      labs(
        title = "Gene set enrichment (Limma::camera)",
        x = "",
        y = "",
        fill = "-logFDR.sclaed"
      ) +
      theme(
        legend.position = "left",
        legend.direction = "horizontal",
        panel.border = element_rect(
          color = "white",
          fill = NA,
          linetype = 1,
          linewidth = 0.2
        ),
        axis.text = element_text(face = 4, size = 7),
        axis.title = element_text(face = 2, size = 7),
        axis.text.x = element_text(
          face = 2,
          size = 7,
          color = "royalblue3",
          vjust = 1,
          angle = 30
        ),
        axis.text.y.right = element_text(
          size = 8,
          color = "gray50",
          face = 2,
          hjust = 0
        ),
        axis.ticks = element_blank(),
        legend.margin = margin(
          t = 2,
          r = 4,
          b = 2,
          l = 2,
          unit = "pt"
        ),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(2.5, "mm"),
        legend.title = element_text(
          face = 4,
          size = 8,
          color = "navy",
          vjust = 1.7,
          hjust = 0
        ),
        plot.title = element_text(
          size = 8,
          color = "gray10",
          face = 4,
          hjust = 0.5,
          vjust = 1
        ),
        plot.background = element_blank()
      ) +
      coord_fixed(ratio = 0.55) +
      guides(y = "none",
             y.sec = guide_axis(),
             fill = guide_colorbar(ticks = FALSE, 
                                   frame.colour = "navy", 
                                   frame.linewidth = 0.4,
                                   title.position = "top",
                                   title.hjust = 0.5,
                                   label.theme = element_text(size= 8, face = 2)))
    
    if (rescale) {
      heatmap_plot <- heatmap_plot +
        colorspace::scale_fill_continuous_diverging(breaks = c(-1, 0, 1))
    }
    
    if (!is.null(additional_ggplot_layers)) {
      heatmap_plot <- heatmap_plot + additional_ggplot_layers
    }
    
    return(heatmap_plot)
    
  }
  
  if (plot.type == "bar") {
    plot_list <- list()
    
    mylist <- GSEA.table(mylist, n = n.top)
    
    for (i in seq_along(mylist)) {
      df <- mylist[[i]]
      
      if (rescale) {
        min_val <- min(df[[y_var]], na.rm = TRUE)
        max_val <- max(df[[y_var]], na.rm = TRUE)
        df[[y_var]] <-
          2 * ((df[[y_var]] - min_val) / (max_val - min_val)) - 1
      }
      
      p <- ggbarplot(
        df,
        x = x_var,
        y = y_var,
        fill = fill_var,
        sort.val = sort_val,
        sort.by.groups = sort_by_groups,
        x.text.angle = 0,
        xlab = xlab,
        ylab = ylab,
        legend = FALSE,
        # Turn off individual legends
        rotate = TRUE,
        width = 1,
        ggtheme = ggtheme,
        color = "white",
        size = 0.5
      ) +
        colorspace::scale_fill_continuous_diverging(
          palette = palette.col,
          breaks = c(-1, 0, 1),
          position = "right"
        ) +
        labs(title = names(mylist)[i], tag = LETTERS[i]) + # Add title from data frame name
        theme(
          legend.key.width = unit(7, "mm"),
          legend.key.height = unit(2.5, "mm"),
          plot.tag = element_text(
            size = 11,
            color = "gray20",
            face = 2
          ),
          plot.tag.position = "topleft",
          # legend.background = element_rect(color = "gray80"),
          legend.margin = margin(
            t = 12,
            r = 0,
            b = 4,
            l = 0,
            unit = "pt"
          ),
          legend.title = element_text(
            face = 4,
            size = 8,
            color = "gray30",
            vjust = 1.85,
            hjust = 0.5
          ),
          axis.text.x = element_text(
            size = 8,
            color = "black",
            face = 2
          ),
          axis.text.y.right = element_text(
            size = 7,
            color = "gray30",
            face = 2,
            hjust = 0
          ),
          legend.text = element_text(
            size = 9,
            face = 2,
            color = "navy"
          ),
          legend.justification = "center",
          axis.ticks = element_blank(),
          plot.title = element_text(
            size = 9,
            color = "gray30",
            face = 4
          ),
          axis.line.y = element_line(color = "white")
        ) +
        guides(y = "none",
               y.sec = guide_axis(),
               fill = guide_colorbar(ticks = FALSE, 
                                     frame.colour = "navy", 
                                     frame.linewidth = 0.4,
                                     title.position = "top",
                                     label.theme = element_text(size= 8, face = 2)))+
        geom_hline(aes(yintercept = 0),
                   color = "gray30",
                   linewidth = 1.1)
      
      if (!is.null(additional_ggplot_layers)) {
        p <- p + additional_ggplot_layers
      }
      
      plot_list[[i]] <- p
    }
    
    combined_plot <- ggarrange(
      plotlist = plot_list,
      ncol = ncol_grid,
      nrow = ceiling(length(plot_list) / 2),
      common.legend = common.legend,
      align = "hv"
    )
    
    combined_plot <- annotate_figure(
      combined_plot,
      top = text_grob(
        "Gene set enrichment (Limma::camera)",
        color = "gray20",
        face = 4,
        size = 12,
        hjust = 0.1,
        x = 0.03
      )
    )
    
    combined_plot <-
      combined_plot + background_grid(major = "xy", minor = "none")
    return(combined_plot)
  }
  
  if (plot.type == "none") {
    cat(yellow$bold(
      paste0(
        "\nPrinting the top ",
        n.top,
        " pathways for each conidtion",
        "\n"
      )
    ))
    
    GSEA.table(mylist, n = n.top)
    
  } else {
    cat(yellow$bold("\nInvalid plot.type argument.\n"))
  }
  
}
