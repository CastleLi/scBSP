scBSP_Pvalues <- read.csv("../../Data/Real Data/HDST/HDST/HDST_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/HDST/HDST/HDST_P_values_SPARK_R.csv")
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
GO_Overlapping <- GO_Comparison(Merged_Pvalues, "../../Output/Real Data/Figure/HDST_Enrich_Plot.png")
#4.200000e-01   6.800000e-01   3.333333e-01   1.500609e-06

# Proportion -------------------------------------------------------------
Enriched_SVG_Prop <- Enriched_Prop(Merged_Pvalues)

write.csv(Enriched_SVG_Prop$Results_Threshold, row.names = FALSE,
          "../../Data/Real Data/HDST_Threshold.csv")

write.csv(Enriched_SVG_Prop$Results_Rank, row.names = FALSE,
          "../../Data/Real Data/HDST_Rank.csv")
