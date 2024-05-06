

scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/Mouse_olfa_s2_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/Mouse_olfa_s2_P_values_SPARK_R.csv")
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
GO_Overlapping <- GO_Comparison(Merged_Pvalues, "../../Output/Real Data/Figure/Mouse_olfa_s2_Enrich_Plot.png")
#5.200000e-01   7.600000e-01   5.074747e-01   7.374269e-14 

# Proportion -------------------------------------------------------------
Enriched_SVG_Prop <- Enriched_Prop(Merged_Pvalues)

write.csv(Enriched_SVG_Prop$Results_Threshold, row.names = FALSE,
          "../../Data/Real Data/Mouse_olfa_s2_Threshold.csv")

write.csv(Enriched_SVG_Prop$Results_Rank, row.names = FALSE,
          "../../Data/Real Data/Mouse_olfa_s2_Rank.csv")
