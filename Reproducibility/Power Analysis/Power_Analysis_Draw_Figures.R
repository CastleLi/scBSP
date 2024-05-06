

source("./Power_Figure.R")

# main --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
OutputPath <- "../../Output/Power_Analysis"





# 2D ----------------------------------------------------------------------


InputPath_All <- c("../../Output/Simulations/Merged_p_values/2D/2D_Sim/")
SubFolders <- list.files(InputPath_All)

for(SubFolder in SubFolders){
  InputPath <- paste0(InputPath_All, SubFolder, "/")
  Method_Labels <- data.frame(Raw = c("sparkx_pvalue",
                                      "spark_pvalue",
                                      "spatialde_pvalue",
                                      "bsp_pvalue",
                                      "MoranI_pvalue",
                                      "nnsvg_pvalue",
                                      "somde_pvalue",
                                      "scbsp_pvalue"),
                              Labels = c("SPARK-X", "SPARK", "SpatialDE", "BSP",
                                         "Moran's I", "nnSVG", "SOMDE", "scBSP"))
  png(file=paste0(OutputPath, "/2D_Sim_", SubFolder ,".png"),
      width=12, height=12, units = "in", res = 600)
  powerplot(InputPath = InputPath, Method_Labels = Method_Labels, Legend_title = "Method", Ncol= 3, Nrow = 3, 
            Col_plate = c("#4DBBD5FF", "#3C5488FF", "#00A087FF", "#E64B35FF", 
                          "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#DC0000FF"))
  dev.off()
}



InputPath <- paste0("../../Output/Simulations/Merged_p_values/2D/2D_Sim_Manuscript/")
Method_Labels <- data.frame(Raw = c("sparkx_pvalue",
                                    "spark_pvalue",
                                    "spatialde_pvalue",
                                    "bsp_pvalue",
                                    "MoranI_pvalue",
                                    "nnsvg_pvalue",
                                    "somde_pvalue",
                                    "scbsp_pvalue"),
                            Labels = c("SPARK-X", "SPARK", "SpatialDE", "BSP",
                                       "Moran's I", "nnSVG", "SOMDE", "scBSP"))
png(file=paste0(OutputPath, "/2D_Sim_Manuscript_", SubFolder ,".png"),
    width=12, height=4, units = "in", res = 600)
powerplot(InputPath = InputPath, Method_Labels = Method_Labels, Legend_title = "Method", Ncol= 3, Nrow = 1,
          Col_plate = c("#4DBBD5FF", "#3C5488FF", "#00A087FF", "#E64B35FF", 
                        "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#DC0000FF"))
dev.off()




# 3D ----------------------------------------------------------------------



# 3D Supp

# rename

InputPath_All <- c("../../Output/Simulations/Merged_p_values/3D/3D_Sim_Output/")

#for(FileName in list.files(paste0(InputPath_All, "Scenario2"))){
#  file.rename(paste0(InputPath_All, "Scenario2/", FileName), 
#              paste0(InputPath_All, "Scenario2/", gsub("0-8_", "0-80_",FileName)))
#}
#for(FileName in list.files(paste0(InputPath_All, "Scenario1"))){
#  for(i in 1:3){
#    file.rename(paste0(InputPath_All, "Scenario1/", FileName), 
#                paste0(InputPath_All, "Scenario1/", gsub(paste0("RW", i, "_3_"), 
#                                                         paste0("RW", i, "_3-0_"),FileName)))
#  }
#}

SubFolders <- list.files(InputPath_All)

for(SubFolder in SubFolders){
  InputPath <- paste0(InputPath_All, SubFolder, "/")
  Method_Labels <- data.frame(Raw = c("sparkx_pvalue",
                                      "bsp_pvalue",
                                      "scbsp_pvalue"),
                              Labels = c("SPARK-X", "BSP", "scBSP"))
  png(file=paste0(OutputPath, "/", "3D_Sim_Supp_", SubFolder ,".png"),
      width=12, height=12, units = "in", res = 600)
  powerplot(InputPath = InputPath, Method_Labels = Method_Labels, Legend_title = "Method", Ncol= 3, Nrow = 3,
            Col_plate = c("#4DBBD5FF", "#E64B35FF", "#DC0000FF", "#3C5488FF",
                          "#00A087FF", "#8491B4FF", "#F39B7FFF", "#7E6148FF"))
  dev.off()
}



# 3D manuscript
InputPath <- c("../../Output/Simulations/Merged_p_values/3D/3D_Standard_Sim_Output/")
#for(FileName in list.files(InputPath)){
#  if(grepl("Discrete", FileName)){
#    file.rename(paste0(InputPath, "/", FileName), 
#                paste0(InputPath, "/", gsub("Discrete", "scenario4_Discrete",FileName)))
#  }
#}


Method_Labels <- data.frame(Raw = c("sparkx_pvalue",
                                    "bsp_pvalue",
                                    "scbsp_pvalue"),
                            Labels = c("SPARK-X", "BSP", "scBSP"))
png(file=paste0(OutputPath, "/", "3D_Sim_Manuscript.png"),
    width=12, height=3, units = "in", res = 600)
powerplot(InputPath = InputPath, Method_Labels = Method_Labels, Legend_title = "Method", Ncol= 4, Nrow = 1,
          Col_plate = c("#4DBBD5FF", "#E64B35FF", "#DC0000FF", "#3C5488FF",
                        "#00A087FF", "#8491B4FF", "#F39B7FFF", "#7E6148FF"))
dev.off()



# 3D Additional Experiments
# discrete
InputPath_All <- c("../../Output/Simulations/Merged_p_values/3D/Adl_Sim/discrete/")
SubFolders <- list.files(InputPath_All)
for(SubFolder in SubFolders){
  InputPath <- paste0(InputPath_All, SubFolder, "/")
  Method_Labels <- data.frame(Raw = c("sparkx_pvalue",
                                      "bsp_pvalue",
                                      "scbsp_pvalue"),
                              Labels = c("SPARK-X", "BSP", "scBSP"))
  png(file=paste0(OutputPath, "/", "3D_Discrete_", SubFolder ,".png"),
      width=9, height=4, units = "in", res = 600)
  powerplot(InputPath = InputPath, Method_Labels = Method_Labels, Legend_title = "Method", Ncol= 2, Nrow = 1,
            Col_plate = c("#4DBBD5FF", "#E64B35FF", "#DC0000FF", "#3C5488FF",
                          "#00A087FF", "#8491B4FF", "#F39B7FFF", "#7E6148FF"))
  dev.off()
}

# dropout
InputPath <- c("../../Output/Simulations/Merged_p_values/3D/Adl_Sim/dropout/")
Method_Labels <- data.frame(Raw = c("sparkx_pvalue",
                                    "bsp_pvalue",
                                    "scbsp_pvalue"),
                            Labels = c("SPARK-X", "BSP", "scBSP"))
png(file=paste0(OutputPath, "/", "3D_Dropout.png"),
    width=12, height=12, units = "in", res = 600)
powerplot(InputPath = InputPath, Method_Labels = Method_Labels, Legend_title = "Method", Ncol= 3, Nrow = 3,
          Col_plate = c("#4DBBD5FF", "#E64B35FF", "#DC0000FF", "#3C5488FF",
                        "#00A087FF", "#8491B4FF", "#F39B7FFF", "#7E6148FF"))
dev.off()


# zaxis
InputPath_All <- c("../../Output/Simulations/Merged_p_values/3D/Adl_Sim/zaxis_res/")
SubFolders <- list.files(InputPath_All)
for(SubFolder in SubFolders){
  InputPath <- paste0(InputPath_All, SubFolder, "/")
  Method_Labels <- data.frame(Raw = c("sparkx_pvalue",
                                      "bsp_pvalue",
                                      "scbsp_pvalue"),
                              Labels = c("SPARK-X", "BSP", "scBSP"))
  png(file=paste0(OutputPath, "/", "3D_Zaxis_", SubFolder ,".png"),
      width=9, height=4, units = "in", res = 600)
  powerplot(InputPath = InputPath, Method_Labels = Method_Labels, Legend_title = "Method", Ncol= 2, Nrow = 1,
            Col_plate = c("#4DBBD5FF", "#E64B35FF", "#DC0000FF", "#3C5488FF",
                          "#00A087FF", "#8491B4FF", "#F39B7FFF", "#7E6148FF"))
  dev.off()
}



