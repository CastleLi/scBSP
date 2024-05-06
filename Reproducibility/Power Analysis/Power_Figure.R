
# function ----------------------------------------------------------------


# prepare power matrix
powermatrix <- function(InputData){
  SVG_index <- 1:1000
  NULL_index <- 1001:nrow(InputData)
  Labels <- colnames(InputData)[-1]
  
  Power_FDR <- lapply(Labels, function(Label){
    if(length(sort(InputData[NULL_index, Label])) == 0){
      LeastData <- c(InputData[SVG_index, Label],
                     InputData[NULL_index, Label])
    }
    else{
      LeastData <- c(InputData[SVG_index, Label],
                     sort(InputData[NULL_index, Label]))
    }
    Power_FDR_matrix <- t(sapply(LeastData, function(Input_pvalue){
      FDR <- sum(InputData[NULL_index, Label] < Input_pvalue) / max((sum(InputData[NULL_index, Label] < Input_pvalue) + 
                                                                       sum(InputData[SVG_index, Label] < Input_pvalue)),
                                                                    1)
      Power <- (length(which(InputData[SVG_index, Label] < Input_pvalue)) / length(SVG_index) )
      return(c(FDR, Power))
    }))
    Power_FDR_matrix <- Power_FDR_matrix[order(Power_FDR_matrix[,1], Power_FDR_matrix[,2]),]
    Power_FDR_matrix <- as.data.frame(Power_FDR_matrix)
    return(Power_FDR_matrix)
  })
  Power_FDR <- do.call("rbind", Power_FDR)
  
  OutputData <- data.frame(Power_FDR)
  OutputData <- cbind(Methods = rep(Labels, rep(nrow(InputData), length(Labels))),
                      OutputData)
  colnames(OutputData) <- c("Methods", "FDR", "Power")
  OutputData$Methods <- as.factor(OutputData$Methods)
  OutputData$FDR <- round(OutputData$FDR/0.002)*0.002
  OutputData <- OutputData[!duplicated(OutputData),]
  OutputData2 <- aggregate(OutputData$Power, 
                           list(OutputData$Methods, OutputData$FDR),
                           max)
  colnames(OutputData2) <- c("Methods", "FDR", "Power")
  OutputData2_comp <- NULL
  for(Method in levels(OutputData2$Methods)){
    Method_Power_Step <- 0
    for(FDR_step in seq(0, 0.2,0.002)){
      Method_Power_check <- which(OutputData2$Methods == Method & OutputData2$FDR == FDR_step)
      if(length(Method_Power_check)==0){
        OutputData2_comp <- rbind(OutputData2_comp, 
                                  c(Method, FDR_step, Method_Power_Step))
      }
      else{
        Method_Power_Step <- OutputData2$Power[Method_Power_check]
        OutputData2_comp <- rbind(OutputData2_comp, 
                                  c(Method, FDR_step, Method_Power_Step))
      }
    }
  }
  OutputData2_comp <- as.data.frame(OutputData2_comp)
  for(ind_col in 2:3){
    OutputData2_comp[,ind_col] <- as.numeric(OutputData2_comp[,ind_col])
  }
  colnames(OutputData2_comp) <- c("Methods", "FDR", "Power")
  OutputData2 <- rbind(OutputData2, OutputData2_comp)
  for(Feature in c("FDR", "Power")){
    OutputData2[,Feature] <- as.numeric(OutputData2[,Feature])
  }
  OutputData2 <- OutputData2[!duplicated(OutputData2),]
  return(OutputData2)
}

# calculate power
powermatrixall <- function(InputPath){
  InputFiles <- list.files(InputPath)
  InputFiles <- InputFiles[which(grepl("1_mergedPvalue", InputFiles))]
  
  PowerMatrixAll_List <- list()
  for(i in 1:length(InputFiles)){
    InputData <- read.csv(paste0(InputPath, InputFiles[i]))
    PowerMatrixAll <- powermatrix(InputData)
    colnames(PowerMatrixAll)[3] <- "Power1"
    for(j in 2:10){
      InputData <- read.csv(paste0(InputPath, 
                                   gsub("1_mergedPvalue",
                                        paste0(j, "_mergedPvalue"),
                                        InputFiles[i])))
      PowerMatrix <- powermatrix(InputData)
      colnames(PowerMatrix)[3] <- paste0("Power",j)
      PowerMatrixAll <- merge(PowerMatrixAll, PowerMatrix, 
                              by = c("Methods", "FDR"),
                              all = TRUE)
    }
    PowerMatrixAll$Power <- rowMeans(PowerMatrixAll[,3:12], na.rm = TRUE)
    PowerMatrixAll <- PowerMatrixAll[!duplicated(PowerMatrixAll),]
    PowerMatrixAll_List[[gsub("1_mergedPvalue.csv",
                              "",
                              InputFiles[i])]] <- PowerMatrixAll
  }
  return(PowerMatrixAll_List)
}



powerplot <- function(InputPath, Method_Labels, Legend_title, Ncol, Nrow, Col_plate = NULL){
  ggplot_matrix <- powermatrixall(InputPath)
  Sub_Figs <- names(ggplot_matrix)
  
  ggplot_list <- list()
  for(Sub_Fig in Sub_Figs){
    PowerMatrixAll <- ggplot_matrix[[Sub_Fig]]
    for(i in 1:nrow(Method_Labels)){
      
    }
    PowerMatrixAll$Methods <- Method_Labels[match(PowerMatrixAll$Methods, Method_Labels[,1]),2]
    PowerMatrixAll$Methods <- factor(PowerMatrixAll$Methods, levels = Method_Labels[,2])
    
    if(is.null(Col_plate)){
      Col_plate <- c("#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#00A087FF",  
                     "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#DC0000FF")
    }
    
    ggplot_list[[Sub_Fig]] <- ggplot(PowerMatrixAll, aes(x=FDR, y=Power, group = Methods, color = Methods)) +
      geom_line(size = 2.0) + 
      scale_color_manual(values=Col_plate[1:length(unique(PowerMatrixAll[,"Methods"]))], name = Legend_title) + 
      scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
      scale_x_continuous(breaks = seq(0, 0.1, 0.05), limits = c(0, 0.1))+
      theme(axis.text = element_text(size=15),
            axis.title = element_text(size=20,face="bold"),
            legend.text = element_text(size=15))
  }
  
  gggraph <- do.call(ggarrange, c(ggplot_list, ncol= Ncol, nrow= Nrow, common.legend = TRUE, legend = "bottom"))
  print(gggraph)
}