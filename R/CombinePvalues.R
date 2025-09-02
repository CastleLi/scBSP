#' @title Combine p-values across multiple samples from scBSP
#' @description
#' Given the results from multiple samples with gene names and p-values, this 
#' function merges them by gene and computes a combined p-value for each gene. 
#' Fisher's method or Stouffer's method can be used.
#'
#' @param list_of_pvalues A list of data.frames, each with columns:
#'   \code{GeneNames}, \code{P_values}.
#' @param method Combination method. One of \code{"fisher"} (default) or
#'   \code{"stouffer"}.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item GeneNames
#'     \item Number_Samples: number of datasets contributing to this gene
#'     \item Calibrated_P_values: the combined p-value
#'   }
#' @importFrom stats pchisq pnorm qnorm
#' @examples
#' df1 <- data.frame(GeneNames = c("A","B","C"),
#'                   P_values = c(0.01, 0.20, 0.03))
#' df2 <- data.frame(GeneNames = c("A","C","D"),
#'                   P_values = c(0.04, 0.10, 0.50))
#' df3 <- data.frame(GeneNames = c("B","C","E"),
#'                   P_values = c(0.05, 0.02, 0.80))
#'
#' CombinePvalues(list(df1, df2, df3), method = "fisher")
#'
#' @export
CombinePvalues <- function(list_of_pvalues, method = c("fisher", "stouffer")) {
  method <- match.arg(method)
  
  merged <- Reduce(function(x, y) merge(x, y, by = "GeneNames", all = TRUE),
                   list_of_pvalues)
  
  pval_cols <- paste0("P_value_", seq_along(list_of_pvalues))
  names(merged)[-1] <- pval_cols
  
  results <- apply(merged[-1], 1, function(pvals) {
    valid_p <- as.numeric(pvals[!is.na(pvals)])
    k <- length(valid_p)
    
    if (k == 0) return(c(Number_Samples = 0, Calibrated_P_values = NA))
    
    if (method == "fisher") {
      stat <- -2 * sum(log(valid_p))
      pval <- 1 - pchisq(stat, 2 * k)
    } else if (method == "stouffer") {
      z <- qnorm(1 - valid_p/2) * sign(0.5 - valid_p)
      z_comb <- sum(z) / sqrt(k)
      pval <- 2 * (1 - pnorm(abs(z_comb)))
    }
    
    c(Number_Samples = k, Calibrated_P_values = pval)
  })
  
  results <- as.data.frame(t(results))
  final <- data.frame(GeneNames = merged$GeneNames, results)
  rownames(final) <- NULL
  final$Number_Samples <- as.integer(final$Number_Samples)
  
  return(final)
}
