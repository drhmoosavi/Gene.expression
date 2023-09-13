convert_geneids <- function(expression_data) {
  if (!requireNamespace("org.Hs.eg.db", quietly=TRUE)) {
    stop("The package org.Hs.eg.db is required for this function to work. Please install it. OR
            set `convert.ids = FALSE`")
  }
  cat(bold(crayon::yellow("Convert Gene IDs to SYMBOLs using `org.Hs.eg.db` package!\n\n")))
  
  rownames(expression_data) <- unname(annotate::getSYMBOL(row.names(expression_data), "org.Hs.eg.db"))
  return(expression_data)
  
}
# gene expression with row ENTREZ IDs 
# convert_geneids(exp_liv.pri)
