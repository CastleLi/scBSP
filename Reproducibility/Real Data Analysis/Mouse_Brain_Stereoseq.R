scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_brain/Mouse_brain_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_brain/Mouse_brain_P_values_SPARK_R.csv")
Merged_Pvalues <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                             scBSP = scBSP_Pvalues$P_values,
                             SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)

rm(list = setdiff(ls(), c("InputExp", "InputCoords", "Merged_Pvalues")))

# GO analysis -------------------------------------------------------------
library(vidger)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
source("./Real_Data_Analysis.R")
GO_Overlapping <- GO_Comparison(Merged_Pvalues, "../../Output/Real Data/Figure/Mouse_brain_Enrich_Plot.png")



