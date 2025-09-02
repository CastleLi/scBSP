# scBSP - A Fast Tool for Single-Cell Spatially Variable Genes Identifications on Large-Scale Spatially Resolved Transcriptomics Data

This package utilizes a granularity-based dimension-agnostic tool, single-cell big-small patch (scBSP), implementing sparse matrix operation and KD-tree/balltree method for distance calculation, for the identification of spatially variable genes on
large-scale data. A corresponding Python library is available at [https://pypi.org/project/scbsp](https://pypi.org/project/scbsp/).

# Installation
This package can be installed on R CRAN
```
install.packages("scBSP")
```

# Usage

```
# Creating coords and expression matrix
Coords <- expand.grid(1:100,1:100, 1:3)
RandFunc <- function(n) floor(10 * stats::rbeta(n, 1, 5))
Raw_Exp <- Matrix::rsparsematrix(nrow = 10^4, ncol = 3*10^4, density = 0.0001, rand.x = RandFunc)

# Excluding low expressed genes
Filtered_ExpMat <- SpFilter(Raw_Exp)
rownames(Filtered_ExpMat) <- paste0("Gene_", 1:nrow(Filtered_ExpMat))

# Computing p-values
P_values <- scBSP(Coords, Filtered_ExpMat)

```

# Reference
Wang, J., Li, J., Kramer, S.T. et al. Dimension-agnostic and granularity-based spatially variable gene identification using BSP. Nat Commun 14, 7367 (2023). https://doi.org/10.1038/s41467-023-43256-5


