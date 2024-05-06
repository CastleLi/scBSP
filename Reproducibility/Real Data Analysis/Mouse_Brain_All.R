
scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_brain/Mouse_brain_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_brain/Mouse_brain_P_values_SPARK_R.csv")

library(Seurat)
library(SeuratData)
library(SeuratObject)
library(peakRAM)
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
                            scBSP_mem = rep(0, 4),
                            SPARKK = rep(0, 4),
                            SPARKK_mem = rep(0, 4))
ExpMatrix_Filtered <- scBSP::SpFilter(brain_anterior1_extracted$ExpMatrix, Threshold = 1)
stxbrain_time$Gene_Count[1] <- dim(ExpMatrix_Filtered)[1]
stxbrain_time$Cell_Count[1] <- dim(ExpMatrix_Filtered)[2]
mem <- peakRAM({})
stxbrain_time$scBSP_mem[1] <- peakRAM({
  stxbrain_time$scBSP[1] <- system.time(brain_anterior1_pvalue_scBSP <- scBSP::scBSP(brain_anterior1_extracted$Coords, ExpMatrix_Filtered))[1]
})$Peak_RAM_Used_MiB
stxbrain_time$SPARKK_mem[1] <- peakRAM({
  stxbrain_time$SPARKK[1] <- system.time(brain_anterior1_pvalue_SPARKX <- SPARK::sparkx(ExpMatrix_Filtered, brain_anterior1_extracted$Coords))[1]
})$Peak_RAM_Used_MiB

ExpMatrix_Filtered <- scBSP::SpFilter(brain_anterior2_extracted$ExpMatrix, Threshold = 1)
stxbrain_time$Gene_Count[2] <- dim(ExpMatrix_Filtered)[1]
stxbrain_time$Cell_Count[2] <- dim(ExpMatrix_Filtered)[2]
stxbrain_time$scBSP_mem[2] <- peakRAM({
  stxbrain_time$scBSP[2] <- system.time(brain_anterior2_pvalue_scBSP <- scBSP::scBSP(brain_anterior2_extracted$Coords, ExpMatrix_Filtered))[1]
})$Peak_RAM_Used_MiB
stxbrain_time$SPARKK_mem[2] <- peakRAM({
  stxbrain_time$SPARKK[2] <- system.time(brain_anterior2_pvalue_SPARKX <- SPARK::sparkx(ExpMatrix_Filtered, brain_anterior2_extracted$Coords))[1]
})$Peak_RAM_Used_MiB

ExpMatrix_Filtered <- scBSP::SpFilter(brain_posterior1_extracted$ExpMatrix, Threshold = 1)
stxbrain_time$Gene_Count[3] <- dim(ExpMatrix_Filtered)[1]
stxbrain_time$Cell_Count[3] <- dim(ExpMatrix_Filtered)[2]
stxbrain_time$scBSP_mem[3] <- peakRAM({
  stxbrain_time$scBSP[3] <- system.time(brain_posterior1_pvalue_scBSP <- scBSP::scBSP(brain_posterior1_extracted$Coords, ExpMatrix_Filtered))[1]
})$Peak_RAM_Used_MiB
stxbrain_time$SPARKK_mem[3] <- peakRAM({
  stxbrain_time$SPARKK[3] <- system.time(brain_posterior1_pvalue_SPARKX <- SPARK::sparkx(ExpMatrix_Filtered, brain_posterior1_extracted$Coords))[1]
})$Peak_RAM_Used_MiB
ExpMatrix_Filtered <- scBSP::SpFilter(brain_posterior2_extracted$ExpMatrix, Threshold = 1)
stxbrain_time$Gene_Count[4] <- dim(ExpMatrix_Filtered)[1]
stxbrain_time$Cell_Count[4] <- dim(ExpMatrix_Filtered)[2]
stxbrain_time$scBSP_mem[4] <- peakRAM({
  stxbrain_time$scBSP[4] <- system.time(brain_posterior2_pvalue_scBSP <- scBSP::scBSP(brain_posterior2_extracted$Coords, ExpMatrix_Filtered))[1]
})$Peak_RAM_Used_MiB
stxbrain_time$SPARKK_mem[4] <- peakRAM({
  stxbrain_time$SPARKK[4] <- system.time(brain_posterior2_pvalue_SPARKX <- SPARK::sparkx(ExpMatrix_Filtered, brain_posterior2_extracted$Coords))[1]
})$Peak_RAM_Used_MiB

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
rownames(stxbrain_time) <- names(stxBrain_pvalues)[1:4]
write.csv(stxbrain_time, 
          "../../Data/Real Data/Mouse_Brain_All_time.csv")


jaccard_similarity <- function(list1, list2) {
  List_intersection <- length(intersect(list1, list2))
  List_union <- length(union(list1, list2))
  return(as.numeric(List_intersection / List_union))
}
jaccard_similarity_matrix <- function(list_all, Method){
  Output_Mat <- matrix(0, nrow = length(list_all), ncol = length(list_all))
  for(i in 1:length(list_all)){
    for(j in 1:length(list_all)){
      Output_Mat[i,j] <- jaccard_similarity(list_all[[i]][which(list_all[[i]][,Method]<0.05),"Genes"], 
                                            list_all[[j]][which(list_all[[j]][,Method]<0.05),"Genes"])
    }
  }
  return(Output_Mat)
}
scBSP_MouseBrain <- jaccard_similarity_matrix(stxBrain_pvalues, "scBSP")
SPARKX_MouseBrain <- jaccard_similarity_matrix(stxBrain_pvalues, "SPARKX")
rownames(scBSP_MouseBrain) <- colnames(scBSP_MouseBrain) <- names(stxBrain_pvalues)
rownames(SPARKX_MouseBrain) <- colnames(SPARKX_MouseBrain) <- names(stxBrain_pvalues)


write.csv(scBSP_MouseBrain, 
          "../../Output/Real Data/Tables/Mouse_Brain_scBSP_similarity.csv")
write.csv(SPARKX_MouseBrain, 
          "../../Output/Real Data/Tables/Mouse_Brain_SPARKX_similarity.csv")




library(UpSetR)
library(ComplexHeatmap)

scBSP_SVGs <- list(anterior1 = stxBrain_pvalues[[1]]$Genes[which(stxBrain_pvalues[[1]]$scBSP<0.05)],
                   anterior2 = stxBrain_pvalues[[2]]$Genes[which(stxBrain_pvalues[[2]]$scBSP<0.05)],
                   posterior1 = stxBrain_pvalues[[3]]$Genes[which(stxBrain_pvalues[[3]]$scBSP<0.05)],
                   posterior2 = stxBrain_pvalues[[4]]$Genes[which(stxBrain_pvalues[[4]]$scBSP<0.05)],
                   Stereoseq = stxBrain_pvalues[[5]]$Genes[which(stxBrain_pvalues[[5]]$scBSP<0.05)])

scBSP_SVGs_Comb <- make_comb_mat(scBSP_SVGs)
scBSP_SVGs_Comb <- normalize_comb_mat(scBSP_SVGs_Comb, full_comb_sets = TRUE)
set_name(scBSP_SVGs_Comb) <- c("Anterior 1", "Anterior 2", "Posterior 1", "Posterior 2", "Stereo-seq")

SPARKX_SVGs <- list(anterior1 = stxBrain_pvalues[[1]]$Genes[which(stxBrain_pvalues[[1]]$SPARKX<0.05)],
                    anterior2 = stxBrain_pvalues[[2]]$Genes[which(stxBrain_pvalues[[2]]$SPARKX<0.05)],
                    posterior1 = stxBrain_pvalues[[3]]$Genes[which(stxBrain_pvalues[[3]]$SPARKX<0.05)],
                    posterior2 = stxBrain_pvalues[[4]]$Genes[which(stxBrain_pvalues[[4]]$SPARKX<0.05)],
                    Stereoseq = stxBrain_pvalues[[5]]$Genes[which(stxBrain_pvalues[[5]]$SPARKX<0.05)])

SPARKX_SVGs_Comb <- make_comb_mat(SPARKX_SVGs)
SPARKX_SVGs_Comb <- normalize_comb_mat(SPARKX_SVGs_Comb, full_comb_sets = TRUE)
set_name(SPARKX_SVGs_Comb) <- c("Anterior 1", "Anterior 2", "Posterior 1", "Posterior 2", "Stereo-seq")


Col_plate <- c("#4DBBD5FF", "#3C5488FF", "#DC0000FF",  "#00A087FF",
               "#7E6148FF", "#8491B4FF", "#E64B35FF","#F39B7FFF")

png(file="../../Output/Real Data/Figure/Mouse_brain_Upset.png",
    width=6, height=6, units = "in", res = 600)

UpSet(scBSP_SVGs_Comb, set_order = c("Anterior 1", "Anterior 2", "Posterior 1", "Posterior 2", "Stereo-seq"),
      top_annotation = upset_top_annotation(scBSP_SVGs_Comb, add_numbers = TRUE),
      comb_col = Col_plate[comb_degree(scBSP_SVGs_Comb)]) %v%
  UpSet(SPARKX_SVGs_Comb, set_order = c("Anterior 1", "Anterior 2", "Posterior 1", "Posterior 2", "Stereo-seq"),
        top_annotation = upset_top_annotation(SPARKX_SVGs_Comb, add_numbers = TRUE),
        comb_col = Col_plate[comb_degree(scBSP_SVGs_Comb)])

dev.off()



# GO analysis -------------------------------------------------------------

library(vidger)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
source("./Real_Data_Analysis.R")

for(i in 1:4){
  GO_Overlapping <- GO_Comparison(stxBrain_pvalues[[i]], 
                                  paste0("../../Output/Real Data/Figure/Mouse_Brain_10x_",
                                         names(stxBrain_pvalues)[i],
                                         
                                         "_Enrich_Plot.png"))
}


#  7.600000e-01   9.000000e-01   5.591919e-01   1.673467e-16 






# Proportion -------------------------------------------------------------
library(vidger)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
source("./Real_Data_Analysis.R")
stxBrain_prop <- lapply(stxBrain_pvalues, Enriched_Prop)

for(i in 1:length(stxBrain_prop)){
  for(j in 1:length(stxBrain_prop[[i]])){
    write.csv(stxBrain_prop[[i]][[j]], row.names = FALSE,
              paste0("../../Data/Real Data/Mouse_Brain_", names(stxBrain_pvalues)[i], 
                     "_", names(stxBrain_prop[[i]])[j], ".csv"))
  }
}

rm(list=ls())


stxBrain_prop <- lapply(c("anterior1", "anterior2", "posterior1", "posterior2", "Stereoseq"), function(FileName){
  Prop_Threshold <- read.csv(paste0("../../Data/Real Data/Mouse_Brain_", FileName, "_Results_Threshold.csv"))
  Prop_Threshold_Rand <- read.csv(paste0("../../Data/Real Data/Mouse_Brain_", FileName, "_Random_Threshold.csv"))
  
  Prop_Threshold_Rand[is.na(Prop_Threshold_Rand)] <- 0
  Prop_Threshold$Baseline_scBSP <- aggregate(Prop_Threshold_Rand$scBSP, list(Prop_Threshold_Rand$Threshold), mean)[,2]
  Prop_Threshold$Baseline_SPARKX <- aggregate(Prop_Threshold_Rand$SPARKX, list(Prop_Threshold_Rand$Threshold), mean)[,2]
  
  return(Prop_Threshold)
})



Col_plate <- c("#4DBBD5FF", "#DC0000FF", "#E64B35FF", "#3C5488FF",  
               "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#00A087FF")
Max_Prop <- 0.15


Plot_anterior1_Threshold <- ggplot(tidyr::pivot_longer(stxBrain_prop[[1]], cols = -Threshold, names_to = "Method"), 
       aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c( "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))


Plot_anterior2_Threshold <-ggplot(tidyr::pivot_longer(stxBrain_prop[[2]], cols = -Threshold, names_to = "Method"), 
                                  aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c( "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))

Plot_posterior1_Threshold <- ggplot(tidyr::pivot_longer(stxBrain_prop[[3]], cols = -Threshold, names_to = "Method"), 
                                    aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c( "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))


Plot_posterior2_Threshold <- ggplot(tidyr::pivot_longer(stxBrain_prop[[4]], cols = -Threshold, names_to = "Method"), 
                                    aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c( "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))


Plot_stereo_Threshold <- ggplot(tidyr::pivot_longer(stxBrain_prop[[5]], cols = -Threshold, names_to = "Method"), 
                                aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c( "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))

png(file="../../Output/Real Data/Figure/Mouse_Brain_All_Enriched_Prop.png",
    width=15, height=3, units = "in", res = 600)

ggpubr::ggarrange(Plot_anterior1_Threshold, 
                  Plot_anterior2_Threshold, 
                  Plot_posterior1_Threshold,
                  Plot_posterior2_Threshold,
                  Plot_stereo_Threshold, 
                  ncol = 5, nrow = 1, 
                  common.legend = TRUE, legend = "bottom")

dev.off()






