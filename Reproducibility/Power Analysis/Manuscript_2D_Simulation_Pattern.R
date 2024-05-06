
library(ggplot2)
library(scales)
library(ggpubr)
Input_data <- read.csv(".\\..\\..\\Data\\2D\\Signal_Str\\sim_MOB_pattern2_fc5_tau50_count_power1.csv")


Figure1 <- ggplot(Input_data, aes(x = x, y = y,
                       color = Input_data[,"gene1"])) +
  geom_point(size = 3) + 
  scale_color_gradientn(colours = c("#FEFEE3", "#1E8654", "#004533"),
                        breaks = quantile(Input_data[,"gene1"], c(0, 0.1, 1)),
                        labels = c("Low", "", "High")) +
  theme_void() +
  labs(color = "Expression") 


Input_data2 <- read.csv(".\\..\\..\\Data\\2D\\Signal_Str\\sim_MOB_pattern3_fc5_tau50_count_power1.csv")
Figure2 <- ggplot(Input_data2, aes(x = x, y = y,
                       color = Input_data2[,"gene1"])) +
  geom_point(size = 3) + 
  scale_color_gradientn(colours = c("#FEFEE3", "#1E8654", "#004533"),
                        breaks = quantile(Input_data[,"gene1"], c(0, 0.1, 1)),
                        labels = c("Low", "", "High")) +
  theme_void() +
  labs(color = "Expression") 



Input_data3 <- read.csv(".\\..\\..\\Data\\2D\\Signal_Str\\sim_MOB_pattern4_fc5_tau50_count_power1.csv")
Figure3 <- ggplot(Input_data3, aes(x = x, y = y,
                       color = Input_data3[,"gene1"])) +
  geom_point(size = 3) + 
  scale_color_gradientn(colours = c("#FEFEE3", "#1E8654", "#004533"),
                        breaks = quantile(Input_data[,"gene1"], c(0, 0.1, 1)),
                        labels = c("Low", "", "High")) +
  theme_void() +
  labs(color = "Expression") 


Final_figure <- ggarrange(Figure1, Figure2, Figure3, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")


png(file=".\\..\\..\\Output\\Power_Analysis\\Manuscript_2D_Simulation_Pattern.png",width = 10,height = 3.5, units = "in",
    res = 600)
print(Final_figure)
dev.off()

