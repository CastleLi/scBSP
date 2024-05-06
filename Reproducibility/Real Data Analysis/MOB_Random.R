

# stereoseq s1 -------------------------------------------------------------
scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/Mouse_olfa_s1_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/Mouse_olfa_s1_P_values_SPARK_R.csv")
Merged_Pvalues <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                             scBSP = scBSP_Pvalues$P_values,
                             SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)

source("./Real_Data_Analysis.R")
Enriched_SVG_Prop_Threshold <- NULL
Enriched_SVG_Prop_Rank <- NULL
for(i in 1:10){
  set.seed(i)
  Enriched_SVG_Prop <- Enriched_Prop_Base(Merged_Pvalues)
  Enriched_SVG_Prop_Threshold <- rbind(Enriched_SVG_Prop_Threshold, Enriched_SVG_Prop$Results_Threshold)
  Enriched_SVG_Prop_Rank <- rbind(Enriched_SVG_Prop_Rank, Enriched_SVG_Prop$Results_Rank)
}
write.csv(Enriched_SVG_Prop_Threshold, row.names = FALSE,
          "../../Data/Real Data/Mouse_olfa_s1_Random_Threshold.csv")
write.csv(Enriched_SVG_Prop_Rank, row.names = FALSE,
          "../../Data/Real Data/Mouse_olfa_s1_Random_Rank.csv")

rm(list=ls())

# stereoseq s2 -------------------------------------------------------------
scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/Mouse_olfa_s2_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/Mouse_olfa_s2_P_values_SPARK_R.csv")
Merged_Pvalues <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                             scBSP = scBSP_Pvalues$P_values,
                             SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)

source("./Real_Data_Analysis.R")
Enriched_SVG_Prop_Threshold <- NULL
Enriched_SVG_Prop_Rank <- NULL
for(i in 1:10){
  set.seed(i)
  Enriched_SVG_Prop <- Enriched_Prop_Base(Merged_Pvalues)
  Enriched_SVG_Prop_Threshold <- rbind(Enriched_SVG_Prop_Threshold, Enriched_SVG_Prop$Results_Threshold)
  Enriched_SVG_Prop_Rank <- rbind(Enriched_SVG_Prop_Rank, Enriched_SVG_Prop$Results_Rank)
}
write.csv(Enriched_SVG_Prop_Threshold, row.names = FALSE,
          "../../Data/Real Data/Mouse_olfa_s2_Random_Threshold.csv")
write.csv(Enriched_SVG_Prop_Rank, row.names = FALSE,
          "../../Data/Real Data/Mouse_olfa_s2_Random_Rank.csv")


rm(list=ls())

# HDST -------------------------------------------------------------

scBSP_Pvalues <- read.csv("../../Data/Real Data/HDST/HDST/HDST_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/HDST/HDST/HDST_P_values_SPARK_R.csv")
Merged_Pvalues <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                             scBSP = scBSP_Pvalues$P_values,
                             SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)
source("./Real_Data_Analysis.R")
Enriched_SVG_Prop_Threshold <- NULL
Enriched_SVG_Prop_Rank <- NULL
for(i in 1:10){
  set.seed(i)
  Enriched_SVG_Prop <- Enriched_Prop_Base(Merged_Pvalues)
  Enriched_SVG_Prop_Threshold <- rbind(Enriched_SVG_Prop_Threshold, Enriched_SVG_Prop$Results_Threshold)
  Enriched_SVG_Prop_Rank <- rbind(Enriched_SVG_Prop_Rank, Enriched_SVG_Prop$Results_Rank)
}

write.csv(Enriched_SVG_Prop_Threshold, row.names = FALSE,
          "../../Data/Real Data/HDST_Random_Threshold.csv")

write.csv(Enriched_SVG_Prop_Rank, row.names = FALSE,
          "../../Data/Real Data/HDST_Random_Rank.csv")

