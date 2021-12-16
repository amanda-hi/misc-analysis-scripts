# Amanda Leonti
# Oct 26, 2018
# Purpose: Filter processed gene expression counts matrix to remove duplicated gene entries.

#####################################################################

filter_dupGenes <- function(counts, gene_id_column = "geneSymbol") {
  # Counts = data frame containing expression data
  # gene_id_column = name of column containing gene ID info
  
  library(dplyr)
  library(stringr)
  library(matrixStats) # for the rowVars function
  library(metaMA) # for the rowVars function
  
  # Getting a list of genes that appear multiple times in the dataset
  duplicatedGenes <- counts[[gene_id_column]][which(duplicated(counts[[gene_id_column]]))]
  
  print(paste("There are", length(unique(duplicatedGenes)), "duplicated genes, for a total of", length(duplicatedGenes), "duplicated entries to remove."))
  
  length(unique(duplicatedGenes))
  length(duplicatedGenes) 
  
  # Filtering out duplicate entries that have low variance (relative to the other duplicates) and only retaining the
  # entry with the highest patient-to-patient variance
  idx <- grep("gene", colnames(counts)) # Getting the index of the gene ID column so it can be excluded from this step
  counts$var <- metaMA::rowVars(counts[,-idx]) # Calculates the rowwise variance, excluding the gene id column
  
  filtered_counts <- counts %>%
    group_by(!!sym(gene_id_column)) %>% 
    filter(!(var < max(var))) %>%
    
    # The above filtering step doesn't remove duplicate entries where the counts are 0 across the board.  
    filter(!(duplicated(!!sym(gene_id_column)))) %>% # This removes any remaining duplicate entries, regardless of variance (since the var is equivalent, i.e. 0)
    ungroup()
  
  filtered_counts <- as.data.frame(filtered_counts) # Casting to a data frame to support rownames (currently it's a tibble)
  rownames(filtered_counts) <- filtered_counts[[gene_id_column]]
  filtered_counts$var <- NULL
  
  if(any(duplicated(rownames(filtered_counts)))){
    print("Duplicated genes remain in the counts matrix - check for errors")
  }else{
    print("No duplicated genes remain!")
  }
  
  return(filtered_counts)
}
  
