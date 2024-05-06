
scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_brain/Mouse_brain_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_brain/Mouse_brain_P_values_SPARK_R.csv")

library(Seurat)
library(SeuratData)
library(SeuratObject)
brain_anterior1 <- LoadData("stxBrain", type = "anterior1")
brain_anterior2 <- LoadData("stxBrain", type = "anterior2")
brain_posterior1 <- LoadData("stxBrain", type = "posterior1")
brain_posterior2 <- LoadData("stxBrain", type = "posterior2")

brain_anterior1 <- SCTransform(brain_anterior1, assay = "Spatial", verbose = FALSE)
brain_anterior2 <- SCTransform(brain_anterior2, assay = "Spatial", verbose = FALSE)
brain_posterior1 <- SCTransform(brain_posterior1, assay = "Spatial", verbose = FALSE)
brain_posterior2 <- SCTransform(brain_posterior2, assay = "Spatial", verbose = FALSE)

brain_anterior1_extracted <- scBSP::LoadSpatial(brain_anterior1)
brain_anterior2_extracted <- scBSP::LoadSpatial(brain_anterior2)
brain_posterior1_extracted <- scBSP::LoadSpatial(brain_posterior1)
brain_posterior2_extracted <- scBSP::LoadSpatial(brain_posterior2)

stxbrain_time <- data.frame(Gene_Count = rep(0, 4),
                            Cell_Count = rep(0 ,4),
                            scBSP = rep(0, 4),
                            SPARKK = rep(0, 4))
ExpMatrix_Filtered <- scBSP::SpFilter(brain_anterior1_extracted$ExpMatrix, Threshold = 1)
stxbrain_time$Gene_Count[1] <- dim(ExpMatrix_Filtered)[1]
stxbrain_time$Cell_Count[1] <- dim(ExpMatrix_Filtered)[2]
stxbrain_time$scBSP[1] <- system.time(brain_anterior1_pvalue_scBSP <- scBSP::scBSP(brain_anterior1_extracted$Coords, ExpMatrix_Filtered))[1]
stxbrain_time$SPARKK[1] <- system.time(brain_anterior1_pvalue_SPARKX <- SPARK::sparkx(ExpMatrix_Filtered, brain_anterior1_extracted$Coords))[1]

ExpMatrix_Filtered <- scBSP::SpFilter(brain_anterior2_extracted$ExpMatrix, Threshold = 1)
stxbrain_time$Gene_Count[2] <- dim(ExpMatrix_Filtered)[1]
stxbrain_time$Cell_Count[2] <- dim(ExpMatrix_Filtered)[2]
stxbrain_time$scBSP[2] <- system.time(brain_anterior2_pvalue_scBSP <- scBSP::scBSP(brain_anterior2_extracted$Coords, ExpMatrix_Filtered))[1]
stxbrain_time$SPARKK[2] <- system.time(brain_anterior2_pvalue_SPARKX <- SPARK::sparkx(ExpMatrix_Filtered, brain_anterior2_extracted$Coords))[1]

ExpMatrix_Filtered <- scBSP::SpFilter(brain_posterior1_extracted$ExpMatrix, Threshold = 1)
stxbrain_time$Gene_Count[3] <- dim(ExpMatrix_Filtered)[1]
stxbrain_time$Cell_Count[3] <- dim(ExpMatrix_Filtered)[2]
stxbrain_time$scBSP[3] <- system.time(brain_posterior1_pvalue_scBSP <- scBSP::scBSP(brain_posterior1_extracted$Coords, ExpMatrix_Filtered))[1]
stxbrain_time$SPARKK[3] <- system.time(brain_posterior1_pvalue_SPARKX <- SPARK::sparkx(ExpMatrix_Filtered, brain_posterior1_extracted$Coords))[1]

ExpMatrix_Filtered <- scBSP::SpFilter(brain_posterior2_extracted$ExpMatrix, Threshold = 1)
stxbrain_time$Gene_Count[4] <- dim(ExpMatrix_Filtered)[1]
stxbrain_time$Cell_Count[4] <- dim(ExpMatrix_Filtered)[2]
stxbrain_time$scBSP[4] <- system.time(brain_posterior2_pvalue_scBSP <- scBSP::scBSP(brain_posterior2_extracted$Coords, ExpMatrix_Filtered))[1]
stxbrain_time$SPARKK[4] <- system.time(brain_posterior2_pvalue_SPARKX <- SPARK::sparkx(ExpMatrix_Filtered, brain_posterior2_extracted$Coords))[1]


stxBrain_pvalues <- list(anterior1 = data.frame(Genes = brain_anterior1_pvalue_scBSP$GeneNames,
                                                scBSP = brain_anterior1_pvalue_scBSP$P_values,
                                                SPARKX = brain_anterior1_pvalue_SPARKX$res_mtest$adjustedPval),
                         anterior2 = data.frame(Genes = brain_anterior2_pvalue_scBSP$GeneNames,
                                                scBSP = brain_anterior2_pvalue_scBSP$P_values,
                                                SPARKX = brain_anterior2_pvalue_SPARKX$res_mtest$adjustedPval),
                         posterior1 = data.frame(Genes = brain_posterior1_pvalue_scBSP$GeneNames,
                                                 scBSP = brain_posterior1_pvalue_scBSP$P_values,
                                                 SPARKX = brain_posterior1_pvalue_SPARKX$res_mtest$adjustedPval),
                         posterior2 = data.frame(Genes = brain_posterior2_pvalue_scBSP$GeneNames,
                                                 scBSP = brain_posterior2_pvalue_scBSP$P_values,
                                                 SPARKX = brain_posterior2_pvalue_SPARKX$res_mtest$adjustedPval),
                         Stereoseq = data.frame(Genes = scBSP_Pvalues$GeneNames,
                                                scBSP = scBSP_Pvalues$P_values,
                                                SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval))

source("./Real_Data_Analysis.R")

for(SeqData in 1:5){
  Enriched_SVG_Prop_Threshold <- NULL
  Enriched_SVG_Prop_Rank <- NULL
  for(i in 1:10){
    set.seed(i)
    Enriched_SVG_Prop <- Enriched_Prop_Base(stxBrain_pvalues[[SeqData]])
    Enriched_SVG_Prop_Threshold <- rbind(Enriched_SVG_Prop_Threshold, Enriched_SVG_Prop$Results_Threshold)
    Enriched_SVG_Prop_Rank <- rbind(Enriched_SVG_Prop_Rank, Enriched_SVG_Prop$Results_Rank)
  }
  write.csv(Enriched_SVG_Prop_Threshold, row.names = FALSE,
            paste0("../../Data/Real Data/Mouse_Brain_", names(stxBrain_pvalues)[SeqData], "_Random_Threshold.csv"))
  write.csv(Enriched_SVG_Prop_Rank, row.names = FALSE,
            paste0("../../Data/Real Data/Mouse_Brain_", names(stxBrain_pvalues)[SeqData], "_Random_Rank.csv"))
  
}

rm(list=ls())

