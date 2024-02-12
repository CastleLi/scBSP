#' @title A function for filtering low expressed genes
#' @usage SpFilter(ExpMat_Sp, Threshold = 5)
#' @param ExpMat_Sp  A sparse, N x M expression matrix in dgCMatrix class with N genes and M spots
#' @param Threshold A threshold set to filter out genes with a total read count below this specified value
#' @return A sparse expression matrix in dgCMatrix class
#' @examples
#' # create a sparse expression matrix
#' Raw_ExpMat <- Matrix::rsparsematrix(nrow = 10000, ncol = 2000, 
#' density = 0.01, rand.x = function(n) rpois(n, 15))
#' Filtered_ExpMat <- SpFilter(Raw_ExpMat)
#' @export

# Define 
SpFilter <- function(ExpMat_Sp, Threshold = 5){
  ExpMat_Out <- ExpMat_Sp[which(Matrix::rowSums(ExpMat_Sp) >= Threshold),]
  return(ExpMat_Out)
}