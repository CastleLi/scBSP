
scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/Mouse_olfa_s1_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/Mouse_olfa_s1_P_values_SPARK_R.csv")
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
GO_Overlapping <- GO_Comparison(Merged_Pvalues, "../../Output/Real Data/Figure/Mouse_olfa_s1_Enrich_Plot.png")
#  7.600000e-01   9.000000e-01   5.591919e-01   1.673467e-16 

# Proportion -------------------------------------------------------------
Enriched_SVG_Prop <- Enriched_Prop(Merged_Pvalues)
write.csv(Enriched_SVG_Prop$Results_Threshold, row.names = FALSE,
          "../../Data/Real Data/Mouse_olfa_s1_Threshold.csv")
write.csv(Enriched_SVG_Prop$Results_Rank, row.names = FALSE,
          "../../Data/Real Data/Mouse_olfa_s1_Rank.csv")

