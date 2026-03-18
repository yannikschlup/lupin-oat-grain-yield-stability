# Author: Lukas Graz             Date: 2024-11-20
# -----------------------------------------------
# given an AEH object with p-values, this function adjusts the p-values 
# with hyp$H_direction != "equal" using the holm method
# -----------------------------------------------



################################
# ADJUST P-VALUES
adj_nested_pvals <- function(pval_list, method) {
  # Flatten the list into one vector while keeping track of lengths
  lengths <- sapply(pval_list, length)
  names_list <- names(pval_list)
  
  # Combine all p-values into one vector
  all_pvals <- unlist(pval_list)
  
  # Adjust p-values
  adjusted_pvals <- p.adjust(all_pvals, method = method)
  
  # Split back into list structure
  result <- vector("list", length(pval_list))
  names(result) <- names_list
  
  # Restore the original structure
  start_idx <- 1
  for (i in seq_along(lengths)) {
    end_idx <- start_idx + lengths[i] - 1
    result[[i]] <- adjusted_pvals[start_idx:end_idx]
    start_idx <- end_idx + 1
  }
  
  return(result)
}

# Example usage:
if(FALSE){
  example_list <- list(
    group1 = c(a=0.01, 0.02, 0.03),
    group2 = c(0.001, c=0.004),
    group3 = c(0.05, 0.06, 0.07, 0.08)
  )

  # Adjust p-values using Benjamini-Hochberg method
  adjusted_list <- adj_nested_pvals(example_list, method = "none")
}

adj_nested_pvals_AEH <- function(AEH){
  PVALS <- lapply(AEH, \(x) if(x$hyp$H_direction != "equal") x$model$p_vals)
  PVALS <- PVALS[!sapply(PVALS, is.null)]
  names(PVALS) <- paste0(names(PVALS),"XseperatorX")
  PVALS <- adj_nested_pvals(PVALS, method = "holm") 
  names(PVALS) <- gsub("XseperatorX", "", names(PVALS))
  PVALS <- lapply(PVALS, \(x) {names(x) <- gsub("^.*XseperatorX.", "", names(x)); x})
  for (nam in names(PVALS)) {
    AEH[[nam]]$model$p_vals <- PVALS[[nam]]
  }

  attr(AEH, "n.adjpvals") <- sum(sapply(PVALS, length))
  AEH
}
