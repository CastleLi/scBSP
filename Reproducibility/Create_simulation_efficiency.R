
dir.create(file.path("./efficiency"), showWarnings = FALSE)
dir.create(file.path("./efficiency/small_sample"), showWarnings = FALSE)
dir.create(file.path("./efficiency/large_sample"), showWarnings = FALSE)



# small sample ------------------------------------------------------------

## Load in raw gene counts

library(Seurat)
library(SeuratData)
library(SeuratObject)
brain1 <- LoadData("stxBrain", type = "anterior1")
dim(brain1@assays$Spatial@counts)
ExpData <- as.numeric(as.matrix(brain1@assays$Spatial@counts))


## Scenario1: varying number of spots
dir.create(file.path("./efficiency/small_sample/Spots"), showWarnings = FALSE)
dir.create(file.path("./efficiency/small_sample/Spots/Exp"), showWarnings = FALSE)
dir.create(file.path("./efficiency/small_sample/Spots/Loc"), showWarnings = FALSE)

### Number of spots ranging from to 500 to 10000
NSpots <- c(500, seq(1000, 10000, 1000))
NGene <- 2*10^4
for(NSpot in NSpots){
  set.seed(1)
  NGene <- 2*10^4
  Coords <- expand.grid(1:ceiling(sqrt(NSpot)), 1:ceiling(sqrt(NSpot)))
  Coords <- Coords[1:NSpot,]
  Output_Exp <- sample(ExpData, NGene * NSpot, replace = TRUE)
  Output_Exp <- matrix(Output_Exp, nrow = NGene, ncol = NSpot)
  # exclude empty cells
  Output_Exp <- Output_Exp[, which(colSums(Output_Exp)!=0)]
  Coords <- Coords[which(colSums(Output_Exp)!=0),]
  # exclude low expressed genes
  Output_Exp <- scBSP::SpFilter(Output_Exp)
  rownames(Output_Exp) <- paste0("Gene_", 1:nrow(Output_Exp))
  Filtered_ExpMat <- Matrix::Matrix(Output_Exp, sparse = TRUE)

  #system.time({P_values <- scBSP::scBSP(Coords, Filtered_ExpMat)})
  #system.time({P_values <- SPARK::sparkx(Filtered_ExpMat, Coords)})
  
  Matrix::writeMM(Filtered_ExpMat, 
                  paste0("./efficiency/small_sample/Spots/Exp/", 
                         format(NGene, scientific = FALSE), "_", 
                         format(NSpot, scientific = FALSE),".mtx"))
  write.csv(Coords, row.names = FALSE,
            paste0("./efficiency/small_sample/Spots/Loc/", 
                   format(NGene, scientific = FALSE), "_", 
                   format(NSpot, scientific = FALSE),".csv"))
}

rm(list=setdiff(ls(), "ExpData"))



## Scenario2: varying number of genes
dir.create(file.path("./efficiency/small_sample/Genes"), showWarnings = FALSE)
dir.create(file.path("./efficiency/small_sample/Genes/Exp"), showWarnings = FALSE)
dir.create(file.path("./efficiency/small_sample/Genes/Loc"), showWarnings = FALSE)

### Number of genes ranging from to 500 to 10000
NGenes <- c(seq(2000, 10000, 2000), seq(15000, 40000, 5000))
for(NGene in NGenes){
  set.seed(1)
  NSpot <- 3000
  Coords <- expand.grid(1:ceiling(sqrt(NSpot)), 1:ceiling(sqrt(NSpot)))
  Coords <- Coords[1:NSpot,]
  Output_Exp <- sample(ExpData, NGene * NSpot, replace = TRUE)
  Output_Exp <- matrix(Output_Exp, nrow = NGene, ncol = NSpot)
  # exclude empty cells
  Output_Exp <- Output_Exp[, which(colSums(Output_Exp)!=0)]
  Coords <- Coords[which(colSums(Output_Exp)!=0),]
  # exclude low expressed genes
  Output_Exp <- scBSP::SpFilter(Output_Exp)
  rownames(Output_Exp) <- paste0("Gene_", 1:nrow(Output_Exp))
  Filtered_ExpMat <- Matrix::Matrix(Output_Exp, sparse = TRUE)
  
  #system.time({P_values <- scBSP::scBSP(Coords, Filtered_ExpMat)})
  #system.time({P_values <- SPARK::sparkx(Filtered_ExpMat, Coords)})
  
  Matrix::writeMM(Filtered_ExpMat, 
                  paste0("./efficiency/small_sample/Genes/Exp/", 
                         format(NGene, scientific = FALSE), "_", 
                         format(NSpot, scientific = FALSE),".mtx"))
  write.csv(Coords, row.names = FALSE,
            paste0("./efficiency/small_sample/Genes/Loc/", 
                   format(NGene, scientific = FALSE), "_", 
                   format(NSpot, scientific = FALSE),".csv"))
}

rm(list=setdiff(ls(), "ExpData"))













# large sample ------------------------------------------------------------
dir.create(file.path("./efficiency/large_sample/Spots"), showWarnings = FALSE)
dir.create(file.path("./efficiency/large_sample/Spots/Exp"), showWarnings = FALSE)
dir.create(file.path("./efficiency/large_sample/Spots/Loc"), showWarnings = FALSE)

# Scenario1: varying number of spots
NSpots <- seq(50000, 400000, 50000)
NGene <- 2*10^4
Sparsity <- .0005
for(NSpot in NSpots){
  set.seed(1)
  Coords <- expand.grid(1:ceiling(sqrt(NSpot)), 1:ceiling(sqrt(NSpot)))
  Coords <- Coords[1:NSpot,]
  RandFunc <- function(n) floor(1 + 10 * stats::rbeta(n, 1, 5))
  Raw_Exp <- Matrix::rsparsematrix(nrow = NGene, ncol = nrow(Coords), density = Sparsity, rand.x = RandFunc)
  Raw_Exp <- Raw_Exp[, which(Matrix::colSums(Raw_Exp)!=0)]
  Filtered_ExpMat <- scBSP::SpFilter(Raw_Exp)
  Coords <- Coords[which(Matrix::colSums(Filtered_ExpMat)!=0),]
  rownames(Filtered_ExpMat) <- paste0("Gene_", 1:nrow(Filtered_ExpMat))
  Matrix::writeMM(Filtered_ExpMat, 
                  paste0("./efficiency/large_sample/Spots/Exp/", 
                         format(NGene, scientific = FALSE), "_", 
                         format(NSpot, scientific = FALSE), "_", 
                         format(Sparsity, scientific = FALSE), ".mtx"))
  write.csv(Coords, row.names = FALSE,
            paste0("./efficiency/large_sample/Spots/Loc/", 
                   format(NGene, scientific = FALSE), "_", 
                   format(NSpot, scientific = FALSE), "_", 
                   format(Sparsity, scientific = FALSE), ".csv"))
}


rm(list=ls())






# Scenario2: varying number of genes
dir.create(file.path("./efficiency/large_sample/Genes"), showWarnings = FALSE)
dir.create(file.path("./efficiency/large_sample/Genes/Exp"), showWarnings = FALSE)
dir.create(file.path("./efficiency/large_sample/Genes/Loc"), showWarnings = FALSE)
NGenes <- c(seq(2000, 10000, 2000), seq(15000, 40000, 5000))
NSpot <- 100000
Sparsity <- .0005
for(NGene in NGenes){
  set.seed(1)
  Coords <- expand.grid(1:ceiling(sqrt(NSpot)), 1:ceiling(sqrt(NSpot)))
  Coords <- Coords[1:NSpot,]
  RandFunc <- function(n) floor(1 + 10 * stats::rbeta(n, 1, 5))
  Raw_Exp <- Matrix::rsparsematrix(nrow = NGene, ncol = floor(nrow(Coords) * 2), density = Sparsity, rand.x = RandFunc)
  Raw_Exp <- Raw_Exp[, which(Matrix::colSums(Raw_Exp)!=0)]
  Raw_Exp <- Raw_Exp[, 1:NSpot]
  Filtered_ExpMat <- scBSP::SpFilter(Raw_Exp)
  Coords <- Coords[which(Matrix::colSums(Filtered_ExpMat)!=0),]
  rownames(Filtered_ExpMat) <- paste0("Gene_", 1:nrow(Filtered_ExpMat))
  Matrix::writeMM(Filtered_ExpMat, 
                  paste0("./efficiency/large_sample/Genes/Exp/", 
                         format(NGene, scientific = FALSE), "_", 
                         format(NSpot, scientific = FALSE), "_", 
                         format(Sparsity, scientific = FALSE) ,".mtx"))
  write.csv(Coords, row.names = FALSE,
            paste0("./efficiency/large_sample/Genes/Loc/", 
                   format(NGene, scientific = FALSE), "_", 
                   format(NSpot, scientific = FALSE), "_", 
                   format(Sparsity, scientific = FALSE), ".csv"))
}

rm(list=ls())



# Scenario3: varying sparsity
dir.create(file.path("./efficiency/large_sample/Sparsity"), showWarnings = FALSE)
dir.create(file.path("./efficiency/large_sample/Sparsity/Exp"), showWarnings = FALSE)
dir.create(file.path("./efficiency/large_sample/Sparsity/Loc"), showWarnings = FALSE)
NGene <- 20000
NSpot <- 100000
Sparsities <- seq(0.0001, 0.001, 0.0001)
for(Sparsity in Sparsities){
  set.seed(1)
  Coords <- expand.grid(1:ceiling(sqrt(NSpot)), 1:ceiling(sqrt(NSpot)))
  Coords <- Coords[1:NSpot,]
  RandFunc <- function(n) floor(1 + 10 * stats::rbeta(n, 1, 5))
  Raw_Exp <- Matrix::rsparsematrix(nrow = NGene, ncol = floor(nrow(Coords) * 1.5), density = Sparsity, rand.x = RandFunc)
  Raw_Exp <- Raw_Exp[, which(Matrix::colSums(Raw_Exp)!=0)]
  Raw_Exp <- Raw_Exp[, 1:NSpot]
  Filtered_ExpMat <- scBSP::SpFilter(Raw_Exp)
  Coords <- Coords[which(Matrix::colSums(Filtered_ExpMat)!=0),]
  rownames(Filtered_ExpMat) <- paste0("Gene_", 1:nrow(Filtered_ExpMat))
  Matrix::writeMM(Filtered_ExpMat, 
                  paste0("./efficiency/large_sample/Sparsity/Exp/", 
                         format(NGene, scientific = FALSE), "_", 
                         format(NSpot, scientific = FALSE), "_", 
                         format(Sparsity, scientific = FALSE) ,".mtx"))
  write.csv(Coords, row.names = FALSE,
            paste0("./efficiency/large_sample/Sparsity/Loc/", 
                   format(NGene, scientific = FALSE), "_", 
                   format(NSpot, scientific = FALSE), "_", 
                   format(Sparsity, scientific = FALSE), ".csv"))
}

