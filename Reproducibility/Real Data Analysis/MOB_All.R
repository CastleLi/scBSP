
# Computational efficiency All ------------------------------------------------


Comp_Eff <- read.csv("../../Data/Real Data/Comp_Eff_Summary.csv")
Comp_Eff$Sequencing <- sapply(Comp_Eff$Dataset, function(i){
  return(strsplit(i, '_')[[1]][1])
})
Comp_Eff_scBSP_Mean <- aggregate(Comp_Eff$scBSP_time, list(Comp_Eff$Sequencing), mean)
Comp_Eff_scBSP_Var <- aggregate(Comp_Eff$scBSP_time, list(Comp_Eff$Sequencing), sd)
Comp_Eff_SPARKK_Mean <- aggregate(Comp_Eff$SPARK_time, list(Comp_Eff$Sequencing), mean)
Comp_Eff_SPARKX_Var <- aggregate(Comp_Eff$SPARK_time, list(Comp_Eff$Sequencing), sd)
Comp_Eff_Gene_Mean <- aggregate(Comp_Eff$Permuted_Gene_Count, list(Comp_Eff$Sequencing), mean)
Comp_Eff_Cell_Mean <- aggregate(Comp_Eff$Cell_Count, list(Comp_Eff$Sequencing), mean)
Comp_Eff_All <- data.frame(Method = Comp_Eff_Gene_Mean[,1],
                           Gene_Count = Comp_Eff_Gene_Mean[,2],
                           Cell_Count = Comp_Eff_Cell_Mean[,2],
                           scBSP_time = Comp_Eff_scBSP_Mean[,2],
                           SPARKX_time = Comp_Eff_SPARKK_Mean[,2])

Comp_Eff_time <- reshape2::melt(Comp_Eff_All, 
                                id = c("Method", "Gene_Count", "Cell_Count"), 
                                measure.vars = c("scBSP_time", "SPARKX_time"))
colnames(Comp_Eff_time) <- c("Sequencing", "GeneCount", "CellCount", "Method", "Time")
Comp_Eff_time$Data <- factor(Comp_Eff_time$Sequencing, levels = c("Visium", "Steroseq", 
                                                            "CosMx", "Xenium", 
                                                            "HDST"))
Comp_Eff_time$Method <- as.factor(Comp_Eff_time$Method)
levels(Comp_Eff_time$Method) <- c("scBSP", "SPARK-X")


library(ggplot2)
Col_plate <- c("#4DBBD5FF", "#3C5488FF", "#DC0000FF", "#E64B35FF", 
               "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#00A087FF")

png(file="../../Output/Real Data/Figure/Real_Data_Overall_RunningTime.png",
    width=4, height=4, units = "in", res = 600)
ggplot(Comp_Eff_time, aes(x = factor(Sequencing, levels = c("Visium", "Steroseq", 
                                                            "CosMx", "Xenium", 
                                                            "HDST")), y = Time)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = Method)) +
  geom_point(data = Comp_Eff_time, aes(x = factor(Sequencing, levels = c("Visium", "Steroseq", 
                                                                         "CosMx", "Xenium", 
                                                                         "HDST")), y = CellCount/1000, group = 1),
             color="#00A087FF",size=5)+
  geom_line(data = Comp_Eff_time, aes(x = factor(Sequencing, levels = c("Visium", "Steroseq", 
                                                                        "CosMx", "Xenium", 
                                                                        "HDST")), y = CellCount/1000, group = 1),
            color="#00A087FF",size=2)+ 
  scale_fill_manual(values = c("#DC0000FF", "#4DBBD5FF")) +
  labs(x = "Sequencing Method",
       y = "Computational time (sec)") + 
  scale_y_continuous(sec.axis=sec_axis(~.*1000,name="Number of spots"))  + 
  theme(legend.position = c(0.1, 0.9))
dev.off()


# Computational efficiency ------------------------------------------------


Comp_Eff <- read.csv("../../Data/Real Data/Real_Data_Computational_Efficiency.csv")
Comp_Eff_time <- reshape2::melt(Comp_Eff, 
                                id = c("Dataset", "Gene_Count", "Cell_Count"), 
                                measure.vars = c("scBSP_time", "SPARKX_time"))
colnames(Comp_Eff_time) <- c("Data", "GeneCount", "CellCount", "Method", "Time")
Comp_Eff_time$Data <- factor(Comp_Eff_time$Data, levels = c("anterior1", "anterior2", 
                                                            "posterior1", "posterior2", 
                                                            "Mouse_brain", "Mouse_olfa_s1", "Mouse_olfa_s2",
                                                            "HDST"))
levels(Comp_Eff_time$Data) <- c("10X Visium \nmouse brain \nanterior 1", "10X Visium \nmouse brain \nanterior 2", 
                                "10X Visium \nmouse brain \nposterior 1", "10X Visium \nmouse brain \nposterior 2", 
                                "Stereo-seq \nmouse brain", 
                                "Stereo-seq \nolfactory bulb \nsection 1", "Stereo-seq \nolfactory bulb \nsection 2",
                                "HDST \nolfactory bulb")
Comp_Eff_time$Method <- as.factor(Comp_Eff_time$Method)
levels(Comp_Eff_time$Method) <- c("scBSP", "SPARK-X")


library(ggplot2)
Col_plate <- c("#4DBBD5FF", "#3C5488FF", "#DC0000FF", "#E64B35FF", 
               "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#00A087FF")

png(file="../../Output/Real Data/Figure/Real_Data_RunningTime.png",
    width=8, height=6, units = "in", res = 600)
ggplot(Comp_Eff_time, aes(x = Data, y = Time)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = Method)) +
  geom_point(data = Comp_Eff_time, aes(x = Data, y = CellCount/3000, group = 1),
            color="#00A087FF",size=5)+
  geom_line(data = Comp_Eff_time, aes(x = Data, y = CellCount/3000, group = 1),
            color="#00A087FF",size=2)+ 
  scale_fill_manual(values = c("#DC0000FF", "#4DBBD5FF")) +
  labs(x = "Data",
       y = "Computational time (sec)") + 
  scale_y_continuous(sec.axis=sec_axis(~.*3000,name="Number of spots"))  + 
  theme(legend.position = c(0.1, 0.9))
dev.off()


# SVG proportion ------------------------------------------------


Stereo_seq_s1_Thres <- read.csv("../../Data/Real Data/Mouse_olfa_s1_Threshold.csv")
Stereo_seq_s1_Rank <- read.csv("../../Data/Real Data/Mouse_olfa_s1_Rank.csv")
Stereo_seq_s2_Thres <- read.csv("../../Data/Real Data/Mouse_olfa_s2_Threshold.csv")
Stereo_seq_s2_Rank <- read.csv("../../Data/Real Data/Mouse_olfa_s2_Rank.csv")
HDST_Thres <- read.csv("../../Data/Real Data/HDST_Threshold.csv")
HDST_Rank <- read.csv("../../Data/Real Data/HDST_Rank.csv")

Stereo_seq_s1_Random_Thres <- read.csv("../../Data/Real Data/Mouse_olfa_s1_Random_Threshold.csv")
Stereo_seq_s1_Random_Rank <- read.csv("../../Data/Real Data/Mouse_olfa_s1_Random_Rank.csv")
Stereo_seq_s2_Random_Thres <- read.csv("../../Data/Real Data/Mouse_olfa_s2_Random_Threshold.csv")
Stereo_seq_s2_Random_Rank <- read.csv("../../Data/Real Data/Mouse_olfa_s2_Random_Rank.csv")
HDST_Random_Thres <- read.csv("../../Data/Real Data/HDST_Random_Threshold.csv")
HDST_Random_Rank <- read.csv("../../Data/Real Data/HDST_Random_Rank.csv")

# merge results from scBSP/SPARK-X and baselines
Stereo_seq_s1_Random_Rank[is.na(Stereo_seq_s1_Random_Rank)] <- 0
Stereo_seq_s1_Rank$Baseline <- aggregate(Stereo_seq_s1_Random_Rank$scBSP, list(Stereo_seq_s1_Random_Rank$Threshold), mean)[,2]
Stereo_seq_s1_Random_Thres[is.na(Stereo_seq_s1_Random_Thres)] <- 0
Stereo_seq_s1_Thres$Baseline_scBSP <- aggregate(Stereo_seq_s1_Random_Thres$scBSP, list(Stereo_seq_s1_Random_Thres$Threshold), mean)[,2]
Stereo_seq_s1_Thres$Baseline_SPARKX <- aggregate(Stereo_seq_s1_Random_Thres$SPARKX, list(Stereo_seq_s1_Random_Thres$Threshold), mean)[,2]

Stereo_seq_s2_Random_Rank[is.na(Stereo_seq_s2_Random_Rank)] <- 0
Stereo_seq_s2_Rank$Baseline <- aggregate(Stereo_seq_s2_Random_Rank$scBSP, list(Stereo_seq_s2_Random_Rank$Threshold), mean)[,2]
Stereo_seq_s2_Random_Thres[is.na(Stereo_seq_s2_Random_Thres)] <- 0
Stereo_seq_s2_Thres$Baseline_scBSP <- aggregate(Stereo_seq_s2_Random_Thres$scBSP, list(Stereo_seq_s2_Random_Thres$Threshold), mean)[,2]
Stereo_seq_s2_Thres$Baseline_SPARKX <- aggregate(Stereo_seq_s2_Random_Thres$SPARKX, list(Stereo_seq_s2_Random_Thres$Threshold), mean)[,2]

HDST_Random_Rank[is.na(HDST_Random_Rank)] <- 0
HDST_Rank$Baseline <- aggregate(HDST_Random_Rank$scBSP, list(HDST_Random_Rank$Threshold), mean)[,2]
HDST_Random_Thres[is.na(HDST_Random_Thres)] <- 0
HDST_Thres$Baseline_scBSP <- aggregate(HDST_Random_Thres$scBSP, list(HDST_Random_Thres$Threshold), mean)[,2]
HDST_Thres$Baseline_SPARKX <- aggregate(HDST_Random_Thres$SPARKX, list(HDST_Random_Thres$Threshold), mean)[,2]


Col_plate <- c("#4DBBD5FF", "#DC0000FF", "#E64B35FF", "#3C5488FF",  
               "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#00A087FF")
Max_Prop <- 0.15

Plot_Stereo_seq_s1_Rank <- ggplot(tidyr::pivot_longer(Stereo_seq_s1_Rank, cols = -Threshold, names_to = "Method"), 
                                  aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
                size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline" = "black", 
                                  "Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                      values =  c("Baseline" = "dotted", 
                                  "Baseline_scBSP" = "dotted",
                                  "Baseline_SPARKX" = "dotted",
                                  "scBSP" = "solid",
                                  "SPARKX" = "solid"),
                      drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  guides(colour = guide_legend(nrow = 2))+
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))

Plot_Stereo_seq_s1_Threshold <- ggplot(tidyr::pivot_longer(Stereo_seq_s1_Thres, cols = -Threshold, names_to = "Method"), 
                                       aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline" = "black", 
                                  "Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline" = "dotted", 
                                    "Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  guides(colour = guide_legend(nrow = 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))

Plot_Stereo_seq_s2_Rank <- ggplot(tidyr::pivot_longer(Stereo_seq_s2_Rank, cols = -Threshold, names_to = "Method"), 
                             aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline" = "black", 
                                  "Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline" = "dotted", 
                                    "Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  guides(colour = guide_legend(nrow = 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))

Plot_Stereo_seq_s2_Threshold <- ggplot(tidyr::pivot_longer(Stereo_seq_s2_Thres, cols = -Threshold, names_to = "Method"), 
                                  aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline" = "black", 
                                  "Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline" = "dotted", 
                                    "Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  guides(colour = guide_legend(nrow = 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))

Plot_HDST_Rank <- ggplot(tidyr::pivot_longer(HDST_Rank, cols = -Threshold, names_to = "Method"), 
                             aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline" = "black", 
                                  "Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline" = "dotted", 
                                    "Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  guides(colour = guide_legend(nrow = 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))

Plot_HDST_Threshold <- ggplot(tidyr::pivot_longer(HDST_Thres, cols = -Threshold, names_to = "Method"), 
                                  aes(x = Threshold)) +
  geom_line(aes(x = Threshold, y = value, 
                colour = factor(Method,
                                c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX")), 
                linetype = factor(Method,
                                  c("Baseline", "Baseline_scBSP", "Baseline_SPARKX", "scBSP", "SPARKX"))), 
            size = 1.5) +
  geom_point(aes(x = Threshold, y = value, colour = Method), size = 3) + 
  scale_colour_manual(name = "Method",
                      values =  c("Baseline" = "black", 
                                  "Baseline_scBSP" = Col_plate[2],
                                  "Baseline_SPARKX" = Col_plate[1],
                                  "scBSP" = Col_plate[2],
                                  "SPARKX" = Col_plate[1]),
                      drop = FALSE)+
  scale_linetype_manual(name = "Method",
                        values =  c("Baseline" = "dotted", 
                                    "Baseline_scBSP" = "dotted",
                                    "Baseline_SPARKX" = "dotted",
                                    "scBSP" = "solid",
                                    "SPARKX" = "solid"),
                        drop = FALSE)+
  ylim(0, round(Max_Prop*1.2, 2))+
  guides(colour = guide_legend(nrow = 2))+
  labs(x = "Top K Genes", y = "Enriched SVG Proportion") +
  theme_minimal()+
  theme(legend.key.size = unit(0.5, 'inch'))

png(file="../../Output/Real Data/Figure/Mouse_olfa_Enriched_Prop.png",
    width=6, height=9, units = "in", res = 600)

ggpubr::ggarrange(Plot_Stereo_seq_s1_Rank, Plot_Stereo_seq_s1_Threshold, 
                  Plot_Stereo_seq_s2_Rank, Plot_Stereo_seq_s2_Threshold, 
                  Plot_HDST_Rank, Plot_HDST_Threshold,
                  ncol = 2, nrow = 3, 
                  common.legend = TRUE, legend = "bottom")

dev.off()








# SVG Inference ------------------------------------------------




scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/Mouse_olfa_s1_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/Mouse_olfa_s1_P_values_SPARK_R.csv")
Merged_Pvalues1 <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                             scBSP = scBSP_Pvalues$P_values,
                             SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)

rm(list = setdiff(ls(), c("Merged_Pvalues1")))

scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/Mouse_olfa_s2_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/Mouse_olfa_s2_P_values_SPARK_R.csv")
Merged_Pvalues2 <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                              scBSP = scBSP_Pvalues$P_values,
                              SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)

rm(list = setdiff(ls(), c("Merged_Pvalues1", "Merged_Pvalues2")))

scBSP_Pvalues <- read.csv("../../Data/Real Data/HDST/HDST/HDST_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/HDST/HDST/HDST_P_values_SPARK_R.csv")
Merged_Pvalues3 <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                             scBSP = scBSP_Pvalues$P_values,
                             SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)

Merged_Pvalues <- list(Merged_Pvalues_s1 = Merged_Pvalues1,
                       Merged_Pvalues_s2 = Merged_Pvalues2,
                       Merged_Pvalues_HDST = Merged_Pvalues3)

rm(list = setdiff(ls(), c("Merged_Pvalues")))

# GO analysis -------------------------------------------------------------
library(vidger)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)



SVG_scBSP_s1 <- Merged_Pvalues$Merged_Pvalues_s1$Genes[which(Merged_Pvalues$Merged_Pvalues_s1$scBSP<0.05)]
# 2139
SVG_scBSP_s2 <- Merged_Pvalues$Merged_Pvalues_s2$Genes[which(Merged_Pvalues$Merged_Pvalues_s2$scBSP<0.05)]
# 2057
SVG_scBSP_HDST <- Merged_Pvalues$Merged_Pvalues_HDST$Genes[which(Merged_Pvalues$Merged_Pvalues_HDST$scBSP<0.05)]
# 2292
SVG_SPARKX_s1 <- Merged_Pvalues$Merged_Pvalues_s1$Genes[which(Merged_Pvalues$Merged_Pvalues_s1$SPARKX<0.05)]
# 5972
SVG_SPARKX_s2 <- Merged_Pvalues$Merged_Pvalues_s2$Genes[which(Merged_Pvalues$Merged_Pvalues_s2$SPARKX<0.05)]
# 7715
SVG_SPARKX_HDST <- Merged_Pvalues$Merged_Pvalues_HDST$Genes[which(Merged_Pvalues$Merged_Pvalues_HDST$SPARKX<0.05)]
# 133



# Upset plot --------------------------------------------------------------




library(UpSetR)
library(ComplexHeatmap)

scBSP_SVGs <- list(Stereo_Seq_s1 = SVG_scBSP_s1,
                   Stereo_Seq_s2 = SVG_scBSP_s2,
                   HDST = SVG_scBSP_HDST)

scBSP_SVGs_Comb <- make_comb_mat(scBSP_SVGs)
scBSP_SVGs_Comb <- normalize_comb_mat(scBSP_SVGs_Comb, full_comb_sets = TRUE)
set_name(scBSP_SVGs_Comb) <- c("HDST", "Stereo-seq Section 1", "Stereo-seq Section 2")

SPARKX_SVGs <- list(Stereo_Seq_s1 = SVG_SPARKX_s1,
                    Stereo_Seq_s2 = SVG_SPARKX_s2,
                    HDST = SVG_SPARKX_HDST)

SPARKX_SVGs_Comb <- make_comb_mat(SPARKX_SVGs)
SPARKX_SVGs_Comb <- normalize_comb_mat(SPARKX_SVGs_Comb, full_comb_sets = TRUE)
set_name(SPARKX_SVGs_Comb) <- c("HDST", "Stereo-seq Section 1", "Stereo-seq Section 2")


Col_plate <- c("#4DBBD5FF", "#3C5488FF", "#DC0000FF", "#E64B35FF", 
               "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#00A087FF")

png(file="../../Output/Real Data/Figure/Olfactory_Upset.png",
    width=6, height=6, units = "in", res = 600)

UpSet(scBSP_SVGs_Comb,set_order = c("Stereo-seq Section 1", "Stereo-seq Section 2", "HDST"),
      top_annotation = upset_top_annotation(scBSP_SVGs_Comb, add_numbers = TRUE),
      comb_col = Col_plate[comb_degree(scBSP_SVGs_Comb)]) %v%
  UpSet(SPARKX_SVGs_Comb, set_order = c("Stereo-seq Section 1", "Stereo-seq Section 2", "HDST"),
        top_annotation = upset_top_annotation(SPARKX_SVGs_Comb, add_numbers = TRUE),
        comb_col = Col_plate[comb_degree(scBSP_SVGs_Comb)])

dev.off()


length(intersect(intersect(scBSP_SVGs[[1]], scBSP_SVGs[[2]]), scBSP_SVGs[[3]]))
length(intersect(intersect(SPARKX_SVGs[[1]], SPARKX_SVGs[[2]]), SPARKX_SVGs[[3]]))
length(unique(c(scBSP_SVGs[[1]], scBSP_SVGs[[2]], scBSP_SVGs[[3]])))
length(unique(c(SPARKX_SVGs[[1]], SPARKX_SVGs[[2]], SPARKX_SVGs[[3]])))
length(intersect(intersect(scBSP_SVGs[[1]], scBSP_SVGs[[2]]), scBSP_SVGs[[3]]))/ length(unique(c(scBSP_SVGs[[1]], scBSP_SVGs[[2]], scBSP_SVGs[[3]])))
length(intersect(intersect(SPARKX_SVGs[[1]], SPARKX_SVGs[[2]]), SPARKX_SVGs[[3]]))/length(unique(c(SPARKX_SVGs[[1]], SPARKX_SVGs[[2]], SPARKX_SVGs[[3]])))


jaccard_similarity <- function(list1, list2) {
  List_intersection <- length(intersect(list1, list2))
  List_union <- length(union(list1, list2))
  return(as.numeric(List_intersection / List_union))
}
jaccard_similarity(scBSP_SVGs[[1]], scBSP_SVGs[[2]])
jaccard_similarity(SPARKX_SVGs[[1]], SPARKX_SVGs[[2]])

jaccard_similarity(scBSP_SVGs[[1]], scBSP_SVGs[[3]])
jaccard_similarity(scBSP_SVGs[[2]], scBSP_SVGs[[3]])
jaccard_similarity(SPARKX_SVGs[[1]], SPARKX_SVGs[[3]])
jaccard_similarity(SPARKX_SVGs[[2]], SPARKX_SVGs[[3]])





# visualization -----------------------------------------------------------
rm(list=ls())

HDST_InputExp <- Matrix::readMM("../../Data/Real Data/HDST/HDST/HDST_CN24-D1_mouseMOB_expression_spa.mtx")
HDST_InputExp <- as(HDST_InputExp, "CsparseMatrix")
HDST_InputExp <- Matrix::t(HDST_InputExp)
HDST_InputCoords <- read.csv("../../Data/Real Data/HDST/HDST/HDST_CN24-D1_mouseMOB_metaData.csv")
InputData_GeneNames <- read.csv("../../Data/Real Data/HDST/HDST/CN24-D1_mouseMOB_geneName.csv", header = FALSE)
rownames(HDST_InputExp) <- as.character(InputData_GeneNames[,1])
HDST_InputCoords <- HDST_InputCoords[, c("pxl_col_in_fullres", "pxl_row_in_fullres")]
scBSP_Pvalues <- read.csv("../../Data/Real Data/HDST/HDST/HDST_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/HDST/HDST/HDST_P_values_SPARK_R.csv")
HDST_Merged_Pvalues <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                                  scBSP = scBSP_Pvalues$P_values,
                                  SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)

Stereoseq_s1_InputExp <- Matrix::readMM("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/Exp_Sparse_Data_counts.mtx")
Stereoseq_s1_InputCoords <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/spa.csv")
scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/Mouse_olfa_s1_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s1/Mouse_olfa_s1_P_values_SPARK_R.csv")
Stereoseq_s1_Merged_Pvalues <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                                          scBSP = scBSP_Pvalues$P_values,
                                          SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)

Stereoseq_s2_InputExp <- Matrix::readMM("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/Exp_Sparse_Data_counts.mtx")
Stereoseq_s2_InputCoords <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/spa.csv")
scBSP_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/Mouse_olfa_s2_P_values_R.csv")
SPARKX_Pvalues <- read.csv("../../Data/Real Data/Steroseq/Mouse_Olfa_s2/Mouse_olfa_s2_P_values_SPARK_R.csv")
Stereoseq_s2_Merged_Pvalues <- data.frame(Genes = scBSP_Pvalues$GeneNames,
                                          scBSP = scBSP_Pvalues$P_values,
                                          SPARKX = SPARKX_Pvalues$res_mtest.adjustedPval)

rm(list = setdiff(ls(), c("HDST_InputExp", "HDST_InputCoords", "HDST_Merged_Pvalues",
                          "Stereoseq_s1_InputExp", "Stereoseq_s1_InputCoords", "Stereoseq_s1_Merged_Pvalues",
                          "Stereoseq_s2_InputExp", "Stereoseq_s2_InputCoords", "Stereoseq_s2_Merged_Pvalues")))

scBSP_SVGs <- intersect(intersect(HDST_Merged_Pvalues$Genes[which(HDST_Merged_Pvalues$scBSP<0.05)],
                                  Stereoseq_s1_Merged_Pvalues$Genes[which(Stereoseq_s1_Merged_Pvalues$scBSP<0.05)]),
                        Stereoseq_s2_Merged_Pvalues$Genes[which(Stereoseq_s2_Merged_Pvalues$scBSP<0.05)])
SPARKX_SVGs <- intersect(intersect(HDST_Merged_Pvalues$Genes[which(HDST_Merged_Pvalues$SPARKX<0.05)],
                                   Stereoseq_s1_Merged_Pvalues$Genes[which(Stereoseq_s1_Merged_Pvalues$SPARKX<0.05)]),
                         Stereoseq_s2_Merged_Pvalues$Genes[which(Stereoseq_s2_Merged_Pvalues$SPARKX<0.05)])

# scBSP only
scBSP_psum <- lapply(scBSP_SVGs, function(GeneName){
  return(c(HDST_Merged_Pvalues$scBSP[which(HDST_Merged_Pvalues$Genes==GeneName)],
           Stereoseq_s1_Merged_Pvalues$scBSP[which(Stereoseq_s1_Merged_Pvalues$Genes==GeneName)],
           Stereoseq_s2_Merged_Pvalues$scBSP[which(Stereoseq_s2_Merged_Pvalues$Genes==GeneName)]))
})
scBSP_psum <- do.call("cbind", scBSP_psum)
scBSP_psum <- as.data.frame(t(scBSP_psum))
rownames(scBSP_psum) <- scBSP_SVGs
scBSP_psum <- rowSums(scBSP_psum)
sort(scBSP_psum[setdiff(scBSP_SVGs, SPARKX_SVGs)])

# scBSP only
scBSP_SVG_Exp_HDST <- data.frame(X = HDST_InputCoords[,1],
                               Y = HDST_InputCoords[,2],
                               Exp = HDST_InputExp[which(rownames(HDST_InputExp) == "Ptprd"),])

scBSP_SVG_Exp_S1 <- data.frame(X = Stereoseq_s1_InputCoords[,1],
                               Y = Stereoseq_s1_InputCoords[,2],
                               Exp = Stereoseq_s1_InputExp[which(Stereoseq_s1_Merged_Pvalues$Genes == "Ptprd"),])

scBSP_SVG_Exp_S2 <- data.frame(X = Stereoseq_s2_InputCoords[,1],
                               Y = Stereoseq_s2_InputCoords[,2],
                               Exp = Stereoseq_s2_InputExp[which(Stereoseq_s2_Merged_Pvalues$Genes == "Ptprd"),])

scBSP_SVG_Exp_HDST$Exp[which(scBSP_SVG_Exp_HDST$Exp>3)] <- 3
png(file="../../Output/Real Data/Figure/scBSP_Only_Ptprd_HDST.png",
    width=18, height=18, units = "in", res = 600)
ggplot(scBSP_SVG_Exp_HDST, aes(x = Y, y = X, colour = Exp))+
  geom_point(size = 8)+
  scale_color_gradientn(colours = c("#FEFEFE05", "#1E865466", "#004533FF", "#004533FF"),
                        breaks = c(0, 1, 2, 3),
                        labels = c("Low", "", "","High")) +
  theme_void()+ 
  theme(legend.position = "None")
dev.off()


png(file="../../Output/Real Data/Figure/scBSP_Only_Ptprd_Stereoseq_s1.png",
    width=18, height=18, units = "in", res = 600)
ggplot(scBSP_SVG_Exp_S1, aes(x = X, y = Y, colour = Exp))+
  geom_point(size = 2)+
  scale_color_gradientn(colours = c("#FEFEE305", "#1E865466", "#004533FF"),
                        breaks = c(0, 1, 2),
                        labels = c("Low", "", "High")) +
  theme_void()+ 
  theme(legend.position = "None")
dev.off()

png(file="../../Output/Real Data/Figure/scBSP_Only_Ptprd_Stereoseq_s2.png",
    width=18, height=18, units = "in", res = 600)
ggplot(scBSP_SVG_Exp_S2, aes(x = Y, y = X, colour = Exp))+
  geom_point(size = 2)+
  scale_color_gradientn(colours = c("#FEFEE305", "#1E865466", "#004533FF"),
                        breaks = c(0, 1, 2),
                        labels = c("Low", "", "High")) +
  theme_void()+ 
  theme(legend.position = "None")
dev.off()

