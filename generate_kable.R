generate_kable <- function(df, col_name = NULL) {
  
  if (is.list(df)) {
    df <- dplyr::bind_rows(df, .id = col_name)
    col_counts <- table(df[[col_name]])  # Count unique factors in the column
  }
  
  table_list <- kableExtra::kbl(df) %>%
    kableExtra::kable_classic(full_width = F, html_font = "Cambria") %>%
    kableExtra::column_spec(1, bold = T, border_right = T)
  
  if (!missing(col_name)) {
    table_list <- table_list %>% kableExtra::group_rows(index = col_counts)
  }
  
  return(table_list)
}
