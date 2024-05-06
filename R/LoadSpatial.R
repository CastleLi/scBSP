#' @title Loading data from a Seurat object or a data frame.
#' @description  A function to load and filter data from a Seurat object or a data frame.
#' @usage LoadSpatial(InputData, Dimension = 2)
#' @param InputData A Seurat spatial object or a M x (D + N) data matrix representing the D-dimensional coordinates and expressions of N genes on M spots. The coordinates should be placed at the first D columns
#' @param Dimension The dimension of coordinates
#' @return A list of two data frame:
#' \item{Coords}{A M x D matrix representing D-dimensional coordinates for M spots}
#' \item{ExpMatrix}{A sparse, N x M expression matrix in dgCMatrix class with N genes and M spots}
#' @export

LoadSpatial <- function(InputData, Dimension = 2){
  # Load in SeuratData
  extract_seurat_sp <- function(InputSeuratData, Slice_Num){
    if("row"%in%colnames(InputSeuratData@images[[Slice_Num]]@coordinates)){
      spmat_seurat <- as.data.frame(cbind(InputSeuratData@images[[Slice_Num]]@coordinates[, c("col", "row")]))
    }
    else{
      spmat_seurat <- as.data.frame(cbind(InputSeuratData@images[[Slice_Num]]@coordinates[, c("x", "y")]))
    }
    colnames(spmat_seurat) <- c("x", "y")
    return(spmat_seurat)
  }
  
  # Load in data
  if(inherits(InputData, "Seurat")){
    N_Slice <- length(InputData@images)
    if(N_Slice == 1){
      Dimension <- 2
      spatialMatrix <- extract_seurat_sp(InputData, 1)
      expMatrix <- InputData@assays$Spatial@data
      geneIndex <- rownames(expMatrix)
    }
    else{
      Dimension <- 3
      spatialMatrix <- lapply(1:N_Slice, function(Slice_Num){
        return(extract_seurat_sp(InputData, Slice_Num))
      })
      spatialMatrix <- do.call("rbind", spatialMatrix)
      spatialMatrix <- cbind(spatialMatrix, InputData$slice)
      colnames(spatialMatrix) <- c("x", "y", "z")
      expMatrix <- InputData@assays$Spatial@data
    }
  }
  else{
    spatialMatrix <- InputData[,1:Dimension]
    expMatrix_ds <- InputData[,-c(1:Dimension)]
    expMatrix <- Matrix::t(Matrix::Matrix(as.matrix(expMatrix_ds), sparse = TRUE))
  }
  Extracated_Data <- list(Coords = spatialMatrix,
                          ExpMatrix = expMatrix)
  return(Extracated_Data)
}




