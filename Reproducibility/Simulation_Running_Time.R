#remotes::install_github("lmweber/nnSVG")

dir.create(file.path("./Outputs"), showWarnings = FALSE)
dir.create(file.path("./Outputs/small_sample"), showWarnings = FALSE)
dir.create(file.path("./Outputs/large_sample"), showWarnings = FALSE)


library(peakRAM)

Exp_Names <- c("small_sample", "large_sample")
for(Exp_Name in Exp_Names){
  
  FolderNames <- list.files(paste0("./efficiency/", Exp_Name))
  
  for(FolderName in FolderNames){
    OutputData <- NULL
    FileNames <- list.files(paste0("./efficiency/", Exp_Name, "/", FolderName, "/Exp/"))
    for(FileName in FileNames){
      OutputData_Unit <- FileName
      # Load in Data
      Parameters <- strsplit(FileName, ".mtx")[[1]]
      Parameters <- strsplit(Parameters, "_")[[1]]
      InputExp <- Matrix::readMM(paste0("./efficiency/", Exp_Name, "/", FolderName, "/Exp/",  FileName))
      Coords <- read.csv(paste0("./efficiency/", Exp_Name, "/", FolderName, "/Loc/",  gsub("mtx", "csv",FileName)))
      InputExp <- as(InputExp, "CsparseMatrix")
      rownames(InputExp) <- paste0("Gene_", 1:nrow(InputExp))
      if(Exp_Name == "small_sample"){
        OutputData_Unit <- c(OutputData_Unit, nrow(InputExp), ncol(InputExp))
      }
      else{
        OutputData_Unit <- c(OutputData_Unit, nrow(InputExp), ncol(InputExp), Matrix::nnzero(InputExp)/nrow(InputExp)/ncol(InputExp))
      }
      
      # initialization
      mem <- peakRAM({
        scBSP_time <- system.time({P_values <- scBSP::scBSP(Coords, InputExp)})
      })
      
      # scBSP
      mem <- peakRAM({
        scBSP_time <- system.time({P_values <- scBSP::scBSP(Coords, InputExp)})
      })
      OutputData_Unit <- c(OutputData_Unit, as.numeric(scBSP_time[1]),mem$Peak_RAM_Used_MiB)
      
      
      # SPARK-X
      mem <- peakRAM({
        SPARKX_time <- system.time({P_values <- SPARK::sparkx(InputExp, Coords)})
      })
      OutputData_Unit <- c(OutputData_Unit, as.numeric(SPARKX_time[1]), mem$Peak_RAM_Used_MiB)
      
      
      # nnSVG
      mem <- peakRAM({
        nnSVG_time <- tryCatch({
          R.utils::withTimeout({
            system.time({P_values_nnSVG <- nnSVG::nnSVG(input = InputExp, 
                                                        spatial_coords = Coords)})
          }, timeout = 3600)
        }, TimeoutException = function(ex) {
          NULL  
        })
      })
      if(is.null(nnSVG_time)){
        OutputData_Unit <- c(OutputData_Unit, NA, NA)
      }
      else{
        OutputData_Unit <- c(OutputData_Unit, as.numeric(nnSVG_time[1]), mem$Peak_RAM_Used_MiB)
      }
      
      
      # SPARK
      mem <- peakRAM({
        SPARK_time <- tryCatch({
          R.utils::withTimeout({
            system.time({
              colnames(InputExp) <- rownames(Coords)
              spark <- SPARK::CreateSPARKObject(counts=InputExp, 
                                                location=Coords,
                                                percentage = 0.1, 
                                                min_total_counts = 10)
              spark@lib_size <- apply(spark@counts, 2, sum)
              ## Estimating Parameter Under Null
              spark <- SPARK::spark.vc(spark, 
                                       covariates = NULL, 
                                       lib_size = spark@lib_size, 
                                       num_core = 1,
                                       verbose = F)
              ## Calculating pval
              spark <- SPARK::spark.test(spark, 
                                         check_positive = T, 
                                         verbose = F)
            })
          }, timeout = 3600)
        }, TimeoutException = function(ex) {
          NULL  
        })
      })
      if(is.null(SPARK_time)){
        OutputData_Unit <- c(OutputData_Unit, NA, NA)
      }
      else{
        OutputData_Unit <- c(OutputData_Unit, as.numeric(SPARK_time[1]), mem$Peak_RAM_Used_MiB)
      }
      
      
      # Moran's I
      mem <- peakRAM({
        MoranI_time <- tryCatch({
          R.utils::withTimeout({
            system.time({
              dists.inv <- 1 / as.matrix(dist(Coords))
              diag(dists.inv) <- 0
              P_values <- sapply(1:nrow(InputExp), function(Gene_Index){
                return(ape::Moran.I(InputExp[Gene_Index, ], dists.inv)$p.value)
              })
            })
          }, timeout = 3600)
        }, TimeoutException = function(ex) {
          NULL  
        })
      })
      if(is.null(MoranI_time)){
        OutputData_Unit <- c(OutputData_Unit, NA, NA)
      }
      else{
        OutputData_Unit <- c(OutputData_Unit, as.numeric(MoranI_time[1]), mem$Peak_RAM_Used_MiB)
      }
      
      
      
      OutputData <- rbind(OutputData, OutputData_Unit)
      
    }
    
    OutputData <- as.data.frame(OutputData)
    if(Exp_Name == "small_sample"){
      colnames(OutputData) <- c("FileName", "N_Genes", "N_Spots", "scBSP_time", "scBSP_mem", "SPARKX_time", "SPARKX_mem", "nnSVG_time", "nnSVG_mem", "SPARK_time", "SPARK_mem", "MoranI_time", "MoranI_mem")
    }
    else{
      colnames(OutputData) <- c("FileName", "N_Genes", "N_Spots", "Sparsity", "scBSP_time", "scBSP_mem", "SPARKX_time", "SPARKX_mem", "nnSVG_time", "nnSVG_mem", "SPARK_time", "SPARK_mem", "MoranI_time", "MoranI_mem")
    }
    
    write.csv(OutputData, row.names = FALSE,
              paste0("./Outputs/", Exp_Name, "/", FolderName, ".csv"))
  }
}


