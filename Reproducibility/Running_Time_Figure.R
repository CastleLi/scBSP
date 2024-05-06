

library(ggplot2)
library(ggpubr)


# Time cost ---------------------------------------------------------------
Col_plate <- c("#DC0000FF", "#4DBBD5FF", "#7E6148FF", "#E64B35FF",
               "#3C5488FF", "#00A087FF", "#8491B4FF","#F39B7FFF")


InputPath <- "./Data/Outputs/large_sample/"
OutputPath <- "./Output/Figures/Efficiency/"

InputData_Gene <- read.csv(paste0(InputPath, "Genes.csv"))
InputData_Spot <- read.csv(paste0(InputPath, "Spots.csv"))
InputData_Spot$N_Spots <- round(InputData_Spot$N_Spots / 1000)
InputData_Sparsity <- read.csv(paste0(InputPath, "Sparsity.csv"))
InputData_Sparsity$Sparsity <- round(InputData_Sparsity$Sparsity * 1000, 4)
InputData <- list(Gene = InputData_Gene,
                  Spot = InputData_Spot,
                  Sparsity = InputData_Sparsity)
ggplot_list <- list()
for(File in c("Gene", "Spot", "Sparsity")){
  InputData_Long <- reshape2::melt(InputData[[File]], id = c("N_Genes", "N_Spots", "Sparsity"), 
                                   measure.vars = c("scBSP_time", "SPARKX_time", "SOMDE_time", "BSP_time"))
  InputData_Long$variable <- factor(InputData_Long$variable, levels = c("scBSP_time",
                                                                        "SPARKX_time",
                                                                        "SOMDE_time",
                                                                        "BSP_time"))
  levels(InputData_Long$variable) <- c("scBSP",
                                       "SPARK-X",
                                       "SOMDE", 
                                       "BSP")
  colnames(InputData_Long) <- c("Gene", "Spot", "Sparsity", "Methods", "Value")
  InputData_Long$IndVar <- InputData_Long[, File]
  ggplot_list[[File]] <- ggplot(InputData_Long, aes(x=IndVar, y=Value, group = Methods, color = Methods)) +
    geom_line(size = 1.0) +
    ylab("Computational Time (sec)") +
    scale_y_log10(limits = c(1, 3600))+
    scale_color_manual(values=Col_plate, drop = FALSE)
  
  if(File == "Gene"){
    ggplot_list[[File]] <- ggplot_list[[File]]
      xlab("Number of Genes")
  }  
  if(File == "Spot"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Spots (k)")+
      scale_x_continuous(breaks=seq(0, 400, 50)) 
    
  }
  
  if(File == "Sparsity"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Density (‰)")+
      scale_x_continuous(breaks=seq(0.1, 1, 0.1)) 
  }
}


gggraph <- do.call(ggarrange, c(ggplot_list, ncol= 3, nrow= 1, common.legend = TRUE, legend = "bottom"))



png(file=paste0(OutputPath, "Time_cost.png"),
    width=9, height=3, units = "in", res = 600)
print(gggraph)
dev.off()




rm(list=ls())


# Memory Usage ---------------------------------------------------------------
Col_plate <- c("#DC0000FF", "#4DBBD5FF", "#7E6148FF", "#E64B35FF",
               "#3C5488FF", "#00A087FF", "#8491B4FF","#F39B7FFF")

InputPath <- "./Data/Outputs/large_sample/"
OutputPath <- "./Output/Figures/Efficiency/"

InputData_Gene <- read.csv(paste0(InputPath, "Genes.csv"))
InputData_Spot <- read.csv(paste0(InputPath, "Spots.csv"))
InputData_Spot$N_Spots <- round(InputData_Spot$N_Spots / 1000)
InputData_Sparsity <- read.csv(paste0(InputPath, "Sparsity.csv"))
InputData_Sparsity$Sparsity <- round(InputData_Sparsity$Sparsity * 1000, 4)
InputData <- list(Gene = InputData_Gene,
                  Spot = InputData_Spot,
                  Sparsity = InputData_Sparsity)
ggplot_list <- list()
for(File in c("Gene", "Spot", "Sparsity")){
  InputData_Long <- reshape2::melt(InputData[[File]], id = c("N_Genes", "N_Spots", "Sparsity"),
                                   measure.vars = c("scBSP_mem", "SPARKX_mem", "SOMDE_mem", "BSP_mem"))
  InputData_Long$variable <- factor(InputData_Long$variable, levels = c("scBSP_mem",
                                                                        "SPARKX_mem",
                                                                        "SOMDE_mem",
                                                                        "BSP_mem"))
  levels(InputData_Long$variable) <- c("scBSP",
                                       "SPARK-X",
                                       "SOMDE",
                                       "BSP")
  colnames(InputData_Long) <- c("Gene", "Spot", "Sparsity", "Methods", "Value")
  InputData_Long$Value <- InputData_Long$Value/1024
  InputData_Long$IndVar <- InputData_Long[, File]
  ggplot_list[[File]] <- ggplot(InputData_Long, aes(x=IndVar, y=Value, group = Methods, color = Methods)) +
    geom_line(size = 1.0) +
    ylab("Memory Usage (GB)") +
    scale_y_log10()+
    scale_color_manual(values=Col_plate, drop = FALSE)
  
  if(File == "Gene"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Genes")
  }  
  if(File == "Spot"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Spots (k)")+
      scale_x_continuous(breaks=seq(0, 400, 50)) 
    
  }
  
  if(File == "Sparsity"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Density (‰)")+
      scale_x_continuous(breaks=seq(0.1, 1, 0.1)) 
  }
}


gggraph <- do.call(ggarrange, c(ggplot_list, ncol= 3, nrow= 1, common.legend = TRUE, legend = "bottom"))



png(file=paste0(OutputPath, "Memory_usage.png"),
    width=9, height=3, units = "in", res = 600)
print(gggraph)
dev.off()




rm(list=ls())



# Time cost small ---------------------------------------------------------------
Col_plate = c("#4DBBD5FF", "#3C5488FF", "#00A087FF", "#E64B35FF", 
              "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#DC0000FF")

InputPath <- "./Data/Outputs/small_sample/"
OutputPath <- "./Output/Figures/Efficiency/"

InputData_Gene <- read.csv(paste0(InputPath, "Genes.csv"))
InputData_Spot <- read.csv(paste0(InputPath, "Spots.csv"))
InputData_Spot$N_Spots <- round(InputData_Spot$N_Spots / 1000)
InputData <- list(Gene = InputData_Gene,
                  Spot = InputData_Spot)
ggplot_list <- list()
for(File in c("Gene", "Spot")){
  InputData_Long <- reshape2::melt(InputData[[File]], id = c("N_Genes", "N_Spots"), 
                                   measure.vars = c("scBSP_time",
                                                    "SPARKX_time",
                                                    "SOMDE_time",
                                                    "BSP_time",
                                                    "nnSVG_time",
                                                    "MoranI_time",
                                                    "SPARK_time",
                                                    "spatialDE_time"))
  InputData_Long$variable <- factor(InputData_Long$variable, levels = c("SPARKX_time",
                                                                        "spatialDE_time",
                                                                        "MoranI_time",
                                                                        "SOMDE_time",
                                                                        "SPARK_time",
                                                                        "BSP_time",
                                                                        "nnSVG_time",
                                                                        "scBSP_time"))
  levels(InputData_Long$variable) <- c("SPARK-X",
                                       "spatialDE",
                                       "Moran's I",
                                       "SOMDE",
                                       "SPARK",
                                       "BSP",
                                       "nnSVG",
                                       "scBSP")
  colnames(InputData_Long) <- c("Gene", "Spot", "Methods", "Value")
  InputData_Long$IndVar <- InputData_Long[, File]
  ggplot_list[[File]] <- ggplot(InputData_Long, aes(x=IndVar, y=Value, group = Methods, color = Methods)) +
    geom_line(size = 1.0) +
    ylab("Computational Time (sec)") +
    scale_y_log10(limits = c(1, 3600))+
    scale_color_manual(values=Col_plate, drop = FALSE)
  
  if(File == "Gene"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Genes")
  }  
  if(File == "Spot"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Spots (k)")+
      scale_x_continuous(breaks=seq(0, 10, 2))
    
  }
  
}


gggraph <- do.call(ggarrange, c(ggplot_list, ncol= 1, nrow= 2, common.legend = TRUE, legend = "bottom"))



png(file=paste0(OutputPath, "Time_cost_small.png"),
    width=4, height=8, units = "in", res = 600)
print(gggraph)
dev.off()




rm(list=ls())




# Memory Usage small ---------------------------------------------------------------
Col_plate = c("#4DBBD5FF", "#3C5488FF", "#00A087FF", "#E64B35FF", 
              "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#DC0000FF")

InputPath <- "./Data/Outputs/small_sample/"
OutputPath <- "./Output/Figures/Efficiency/"

InputData_Gene <- read.csv(paste0(InputPath, "Genes.csv"))
InputData_Spot <- read.csv(paste0(InputPath, "Spots.csv"))
InputData_Spot$N_Spots <- round(InputData_Spot$N_Spots / 1000)
InputData <- list(Gene = InputData_Gene,
                  Spot = InputData_Spot)
ggplot_list <- list()

for(File in c("Gene", "Spot")){
  InputData_Long <- reshape2::melt(InputData[[File]], id = c("N_Genes", "N_Spots"),
                                   measure.vars = c("SPARKX_mem",
                                                    "spatialDE_mem",
                                                    "MoranI_mem",
                                                    "SOMDE_mem",
                                                    "SPARK_mem",
                                                    "BSP_mem",
                                                    "nnSVG_mem",
                                                    "scBSP_mem"))
  InputData_Long$variable <- factor(InputData_Long$variable, levels = c("SPARKX_mem",
                                                                        "spatialDE_mem",
                                                                        "MoranI_mem",
                                                                        "SOMDE_mem",
                                                                        "SPARK_mem",
                                                                        "BSP_mem",
                                                                        "nnSVG_mem",
                                                                        "scBSP_mem"))
  levels(InputData_Long$variable) <- c("SPARK-X",
                                       "spatialDE",
                                       "Moran's I",
                                       "SOMDE",
                                       "SPARK",
                                       "BSP",
                                       "nnSVG",
                                       "scBSP")
  colnames(InputData_Long) <- c("Gene", "Spot", "Methods", "Value")
  InputData_Long$Value <- InputData_Long$Value/1024
  InputData_Long$IndVar <- InputData_Long[, File]
  ggplot_list[[File]] <- ggplot(InputData_Long, aes(x=IndVar, y=Value, group = Methods, color = Methods)) +
    geom_line(size = 1.0) +
    ylab("Memory Usage (GB)") +
    scale_y_log10()+
    scale_color_manual(values=Col_plate, drop = FALSE)
  
  if(File == "Gene"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Genes")
  }  
  if(File == "Spot"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Spots (k)")+
      scale_x_continuous(breaks=seq(0, 10, 2)) 
    
  }
}


gggraph <- do.call(ggarrange, c(ggplot_list, ncol= 1, nrow= 2, common.legend = TRUE, legend = "bottom"))



png(file=paste0(OutputPath, "Memory_usage_small.png"),
    width=4, height=8, units = "in", res = 600)
print(gggraph)
dev.off()




rm(list=ls())

















# Repeated Time ---------------------------------------------------------------
Col_plate <- c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF", 
               "#8491B4FF","#7E6148FF", "#F39B7FFF",  "#DC0000FF")


InputPath <- "./Data/Outputs_Rep/large_sample/"
OutputPath <- "./Output/Figures/Efficiency/"
InputData_Gene <- read.csv(paste0(InputPath, "Genes.csv"))
InputData_Spot <- read.csv(paste0(InputPath, "Spots.csv"))
InputData_Spot$N_Spots <- round(InputData_Spot$N_Spots / 1000)
InputData_Sparsity <- read.csv(paste0(InputPath, "Sparsity.csv"))
InputData_Sparsity$Sparsity <- round(InputData_Sparsity$Sparsity * 1000, 4)
InputData <- list(Gene = InputData_Gene,
                  Spot = InputData_Spot,
                  Sparsity = InputData_Sparsity)

ggplot_list <- list()
for(File in c("Gene", "Spot", "Sparsity")){
  InputData_Long <- reshape2::melt(InputData[[File]], id = c("N_Genes", "N_Spots", "Sparsity"), measure.vars = paste0("scBSP_time_", 1:10))
  colnames(InputData_Long) <- c("Gene", "Spot", "Sparsity", "Methods", "Value")
  InputData_Long$IndVar <- InputData_Long[, File]
  ggplot_list[[File]] <- ggplot(InputData_Long, aes(x=IndVar, y=Value)) +
    geom_point(size = 1.0) +
    geom_smooth(method = "lm")+
    ylab("Computational Time (sec)") +
    scale_color_manual(values=Col_plate[1:3], drop = FALSE)
  if(File == "Gene"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Genes")+
      annotate("text", x = 20000, y = 5.5, 
               label = paste0("Pearson correlation = ", round(cor(InputData_Long$Value, InputData_Long$Gene), 4),
                              "\n95% CI = [", format(round(cor.test(InputData_Long$Value, InputData_Long$Gene, method = "pearson", conf.level = 0.95)$conf.int[1], 4), nsmall = 4),
                              ", ",format(round(cor.test(InputData_Long$Value, InputData_Long$Gene, method = "pearson", conf.level = 0.95)$conf.int[2], 4), nsmall = 4), "]"), 
               color = "black", size = 3)
  }  
  if(File == "Spot"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Spots (k)")+
      scale_x_continuous(breaks=seq(0, 400, 50)) +
      annotate("text", x = 200, y = 13, 
               label = paste0("Pearson correlation = ", round(cor(InputData_Long$Value, InputData_Long$Spot), 4),
                              "\n95% CI = [", format(round(cor.test(InputData_Long$Value, InputData_Long$Spot, method = "pearson", conf.level = 0.95)$conf.int[1], 4), nsmall = 4),
                              ", ",format(round(cor.test(InputData_Long$Value, InputData_Long$Spot, method = "pearson", conf.level = 0.95)$conf.int[2], 4), nsmall = 4), "]"), 
               color = "black", size = 3)
    
  }
  
  if(File == "Sparsity"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Density (‰)")+
      scale_x_continuous(breaks=seq(0.1, 1, 0.1)) +
      annotate("text", x = 0.5, y = 6, 
               label = paste0("Pearson correlation = ", round(cor(InputData_Long$Value, InputData_Long$Sparsity), 4),
                              "\n95% CI = [", format(round(cor.test(InputData_Long$Value, InputData_Long$Sparsity, method = "pearson", conf.level = 0.95)$conf.int[1], 4), nsmall = 4),
                              ", ",format(round(cor.test(InputData_Long$Value, InputData_Long$Sparsity, method = "pearson", conf.level = 0.95)$conf.int[2], 4), nsmall = 4), "]"), 
               color = "black", size = 3)
  }
}


gggraph <- do.call(ggarrange, c(ggplot_list, ncol= 3, nrow= 1, common.legend = TRUE, legend = "bottom"))



png(file=paste0(OutputPath, "Time_cost_repeated.png"),
    width=9, height=3, units = "in", res = 600)
print(gggraph)
dev.off()






# Repeated Memory ---------------------------------------------------------------
Col_plate <- c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF", 
               "#8491B4FF","#7E6148FF", "#F39B7FFF",  "#DC0000FF")


InputPath <- "./Data/Outputs_Rep/large_sample/"
OutputPath <- "./Output/Figures/Efficiency/"
InputData_Gene <- read.csv(paste0(InputPath, "Genes.csv"))
InputData_Spot <- read.csv(paste0(InputPath, "Spots.csv"))
InputData_Spot$N_Spots <- round(InputData_Spot$N_Spots / 1000)
InputData_Sparsity <- read.csv(paste0(InputPath, "Sparsity.csv"))
InputData_Sparsity$Sparsity <- round(InputData_Sparsity$Sparsity*1000, 4)
InputData <- list(Gene = InputData_Gene,
                  Spot = InputData_Spot,
                  Sparsity = InputData_Sparsity)

ggplot_list <- list()
for(File in c("Gene", "Spot", "Sparsity")){
  InputData_Long <- reshape2::melt(InputData[[File]], id = c("N_Genes", "N_Spots", "Sparsity"), measure.vars = paste0("scBSP_mem_", 1:10))
  colnames(InputData_Long) <- c("Gene", "Spot", "Sparsity", "Methods", "Value")
  InputData_Long$IndVar <- InputData_Long[, File]
  ggplot_list[[File]] <- ggplot(InputData_Long, aes(x=IndVar, y=Value)) +
    geom_point(size = 1.0) +
    geom_smooth(method = "lm")+
    ylab("Memory Usage (MB)") +
    scale_color_manual(values=Col_plate[1:3], drop = FALSE)
  
  if(File == "Gene"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Genes")+
      annotate("text", x = 20000, y = 790, 
               label = paste0("Pearson correlation = ", round(cor(InputData_Long$Value, InputData_Long$Gene), 4),
                              "\n95% CI = [", format(round(cor.test(InputData_Long$Value, InputData_Long$Gene, method = "pearson", conf.level = 0.95)$conf.int[1], 4), nsmall = 4),
                              ", ",format(round(cor.test(InputData_Long$Value, InputData_Long$Gene, method = "pearson", conf.level = 0.95)$conf.int[2], 4), nsmall = 4), "]"), 
               color = "black", size = 3)
  }  
  if(File == "Spot"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Number of Spots (k)")+
      scale_x_continuous(breaks=seq(0, 400, 50)) +
      annotate("text", x = 200, y = 1750, 
               label = paste0("Pearson correlation = ", round(cor(InputData_Long$Value, InputData_Long$Spot), 4),
                              "\n95% CI = [", format(round(cor.test(InputData_Long$Value, InputData_Long$Spot, method = "pearson", conf.level = 0.95)$conf.int[1], 4), nsmall = 4),
                              ", ",format(round(cor.test(InputData_Long$Value, InputData_Long$Spot, method = "pearson", conf.level = 0.95)$conf.int[2], 4), nsmall = 4), "]"), 
               color = "black", size = 3)
    
  }
  
  if(File == "Sparsity"){
    ggplot_list[[File]] <- ggplot_list[[File]] + 
      xlab("Density (‰)")+
      scale_x_continuous(breaks=seq(0.1, 1, 0.1)) +
      annotate("text", x = 0.5, y = 770, 
               label = paste0("Pearson correlation = ", round(cor(InputData_Long$Value, InputData_Long$Sparsity), 4),
                              "\n95% CI = [", format(round(cor.test(InputData_Long$Value, InputData_Long$Sparsity, method = "pearson", conf.level = 0.95)$conf.int[1], 4), nsmall = 4),
                              ", ",format(round(cor.test(InputData_Long$Value, InputData_Long$Sparsity, method = "pearson", conf.level = 0.95)$conf.int[2], 4), nsmall = 4), "]"), 
               color = "black", size = 3)
  }
}


gggraph <- do.call(ggarrange, c(ggplot_list, ncol= 3, nrow= 1, common.legend = TRUE, legend = "bottom"))



png(file=paste0(OutputPath, "Memory_usage_repeated.png"),
    width=9, height=3, units = "in", res = 600)
print(gggraph)
dev.off()

