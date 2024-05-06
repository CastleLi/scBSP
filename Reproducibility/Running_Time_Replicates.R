#remotes::install_github("lmweber/nnSVG")

dir.create(file.path("./Outputs_Rep"), showWarnings = FALSE)
dir.create(file.path("./Outputs_Rep/large_sample"), showWarnings = FALSE)


library(peakRAM)

Exp_Names <- c("large_sample")
for(Exp_Name in Exp_Names){
  
  FolderNames <- list.files(paste0("./efficiency/", Exp_Name))
  
  for(FolderName in FolderNames[2:3]){
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
      for(i in 1:10){
        mem <- peakRAM({
          scBSP_time <- system.time({P_values <- scBSP::scBSP(Coords, InputExp)})
        })
        OutputData_Unit <- c(OutputData_Unit, as.numeric(scBSP_time[1]),mem$Peak_RAM_Used_MiB)
        
      }
      
      OutputData <- rbind(OutputData, OutputData_Unit)
    }
    
    OutputData <- as.data.frame(OutputData)
    colnames(OutputData) <- c("FileName", "N_Genes", "N_Spots", "Sparsity",as.character(unlist(outer(c("scBSP_time", "scBSP_mem"), paste0("_",1:10), paste0))))
    
    write.csv(OutputData, row.names = FALSE,
              paste0("./Outputs_Rep/", Exp_Name, "/", FolderName, ".csv"))
  }
}


