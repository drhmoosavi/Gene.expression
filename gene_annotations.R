gene_annotations <- function() {
  # Check and load required packages
  required_packages <- c("org.Hs.eg.db", "biomaRt", "dplyr", "AnnotationDbi")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg , "is not installed. Please install it to use this function."))
    }
    library(pkg, character.only = TRUE)
  }
  
  cat(bold(green("Fetching gene information from org.Hs.eg.db...\n")))
  genes <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                                 columns = c("ENTREZID", "SYMBOL", "ENSEMBL", "MAP"),
                                 keytype = "ENTREZID")
  
  cat(bold("Connecting to Ensembl...\n"))
  #  head(listDatasets(ensembl))
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", verbose = T)
  
  dbname <- listDatasets(ensembl) %>% 
    filter(if_any(everything(), ~ str_detect(., "hsapiens")))
  cat(crayon::green("BioMart database name:\n"), crayon::italic(paste(dbname)) ,"\n\n")
  
  cat(crayon::green("Ensembl version:"), "\n")
  print(listEnsembl()[1,])
  
  
  cat(crayon::green("\nGetting gene biotypes from Ensembl...\n" %+% bold(yellow("It takes a minute...\n\n"))))
  genes_with_biotypes <- getBM(attributes = c("entrezgene_id", "gene_biotype"),
                               filters = "entrezgene_id",
                               values = genes$ENTREZID,
                               mart = ensembl,
                               uniqueRows = TRUE,
                               useCache = FALSE)
  
  genes_with_biotypes$entrezgene_id <- as.character(genes_with_biotypes$entrezgene_id)
  
  cat("Joining gene information with biotypes...\n")
  final_genes <- dplyr::inner_join(genes, genes_with_biotypes,
                                   by = c("ENTREZID" = "entrezgene_id"), 
                                   relationship = "many-to-many")
  
  # Check for NA values
  if (any(is.na(final_genes))) {
    cat(crayon::yellow("Note: There are NA values in the final dataset.\n\n"))
  }
  
  final_genes <- final_genes %>%
    dplyr::mutate(gene_biotype = replace_na(gene_biotype,"NA") ) 
  
  cat("Gene biotype distribution:\n\n")
  
  final_genes %>%
    group_by(gene_biotype) %>%
    summarise(count = n_distinct(ENTREZID)) %>% 
    arrange(-count) %>% 
    setNames(c("Gene_type","Number")) %>% 
    print(n=length(unique(final_genes$gene_biotype)))
  
  return(final_genes)
}
