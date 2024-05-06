
library(ggplot2)
library(scales)
library(rayshader)

# Load in raw figure ------------------------------------------------------


# raw figure
Figure_Raw <- png::readPNG("./Output/Fig1_Pattern/Workflow_Pattern.png")
Figure_Gray <- 0.2989 * Figure_Raw[,,1] + 0.587 * Figure_Raw[,,2] + 0.114 * Figure_Raw[,,3]

Figure_Gray <- as.data.frame(Figure_Gray)
Figure <- data.frame(x = rep(1:ncol(Figure_Gray), rep(nrow(Figure_Gray), ncol(Figure_Gray))),
                     y = rep(c((nrow(Figure_Gray) + 1) - 1:nrow(Figure_Gray)), ncol(Figure_Gray)), 
                     value = as.numeric(unlist(as.data.frame(Figure_Gray))))

ggplot(Figure, aes(x = x, y = y, color = value)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "black")



# Plot of raw figure ------------------------------------------------------


# Define regions
## backgrounds
InputExp <- Figure
InputExp$x <- floor(InputExp$x / (max(InputExp$x) - min(InputExp$x)) * 30)
InputExp$y <- floor(InputExp$y / (max(InputExp$y) - min(InputExp$y)) * 30)
InputExp <- InputExp[!duplicated(InputExp[,c("x", "y")]),]
set.seed(1)
InputExp$Exp <- sapply(InputExp$value, function(value){
  if(value > 0.9){
    return(NA)
  }
  else if(value > 0.3){
    return(rpois(1, 10))
  }
  else{
    return(rpois(1, 30))
  }
})
InputExp <- na.omit(InputExp)
for(i in 1:4){
  InputExp[,i] <- as.numeric(InputExp[,i])
}

InputExp <- na.omit(InputExp)

# 2D plot
png(file=".\\Output\\Fig1_Pattern\\Figure1_Raw_2d.png",width = 6,height = 6,units = "in",
    res = 600)
ggplot(InputExp, aes(x = x, y = y, color = Exp)) +
  geom_point(size = 5) +
  scale_color_viridis_c() + 
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none")
dev.off()


# 3D plot
ggvolcano <- InputExp %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = Exp)) +
#  geom_contour(aes(x = x, y = y, z = Exp), breaks = c(0,  20, 100)) +
  scale_x_continuous("X", expand = c(0, 0)) +
  scale_y_continuous("Y", expand = c(0, 0)) +
  scale_fill_gradientn("Z", colours = viridis::viridis(10)) +
  coord_fixed() + theme(line = element_blank(),
                        text = element_blank(),
                        title = element_blank(), 
                        panel.background = element_rect(fill = "white"),
                        legend.position = "none")

plot_gg(ggvolcano, multicore = TRUE, raytrace = TRUE, width = 7, height = 4, 
        scale = 300, windowsize = c(1000, 1000), zoom = 0.6, phi = 30, theta = 30)



# Plot of local means -----------------------------------------------------


LocalMeans <- function(DK) {
  KDBinary <- RANN::nn2(InputExp[,c("x","y")], k = K_NN, treetype = treetype, 
                        searchtype = "radius", radius = DK)$nn.idx
  KDBinary <- Matrix::Matrix(KDBinary, sparse = TRUE)
  KDBinary_rowp <- sparseMatrixStats::rowCounts(KDBinary > 
                                                  0)
  KDBinary_rowind <- rep(1:length(KDBinary_rowp), KDBinary_rowp)
  KDBinary <- spam::as.spam.dgCMatrix(KDBinary)
  PatchesCells <- Matrix::sparseMatrix(x = rep(1, length(KDBinary@entries)), 
                                       i = KDBinary_rowind, j = KDBinary@entries)
  PatchesCells_Centroid <- Matrix::Diagonal(x = (Matrix::colSums(PatchesCells) > 
                                                   1))
  PatchesCells <- PatchesCells - PatchesCells_Centroid
  diagMatrix_sparse <- Matrix::Diagonal(x = 1/Matrix::colSums(PatchesCells))
  X_kj_Inter <- PatchesCells %*% diagMatrix_sparse
  X_kj <- InputExp$Exp %*% X_kj_Inter
  return(X_kj)
}


LocalMeans <- function(DK) {
  KDBinary <- RANN::nn2(InputExp[,c("x","y")], k = 100, treetype = "kd", 
                        searchtype = "radius", radius = DK)$nn.idx
  KDBinary <- Matrix::Matrix(KDBinary, sparse = TRUE)
  KDBinary_rowp <- sparseMatrixStats::rowCounts(KDBinary > 
                                                  0)
  KDBinary_rowind <- rep(1:length(KDBinary_rowp), KDBinary_rowp)
  KDBinary <- spam::as.spam.dgCMatrix(KDBinary)
  PatchesCells <- Matrix::sparseMatrix(x = rep(1, length(KDBinary@entries)), 
                                       i = KDBinary_rowind, j = KDBinary@entries)
  PatchesCells_Centroid <- Matrix::Diagonal(x = (Matrix::colSums(PatchesCells) > 
                                                   1))
  PatchesCells <- PatchesCells - PatchesCells_Centroid
  diagMatrix_sparse <- Matrix::Diagonal(x = 1/Matrix::colSums(PatchesCells))
  X_kj_Inter <- PatchesCells %*% diagMatrix_sparse
  X_kj <- InputExp$Exp %*% X_kj_Inter
  return(as.numeric(X_kj))
}

InputExp$Exp_SG1 <- LocalMeans(1.0)
InputExp$Exp_SG3 <- LocalMeans(3.0)



# 2D plot
png(file=".\\Output\\Fig1_Pattern\\Figure1_SG1_2d.png",width = 6,height = 6,units = "in",
    res = 600)
ggplot(InputExp, aes(x = x, y = y, color = Exp_SG1)) +
  geom_point(size = 5) +
  scale_color_viridis_c() + 
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none")
dev.off()

png(file=".\\Output\\Fig1_Pattern\\Figure1_SG3_2d.png",width = 6,height = 6,units = "in",
    res = 600)
ggplot(InputExp, aes(x = x, y = y, color = Exp_SG3)) +
  geom_point(size = 5) +
  scale_color_viridis_c() + 
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none")
dev.off()




# 3D plot
ggvolcano <- InputExp %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = Exp_SG1)) +
  #  geom_contour(aes(x = x, y = y, z = Exp), breaks = c(0,  20, 100)) +
  scale_x_continuous("X", expand = c(0, 0)) +
  scale_y_continuous("Y", expand = c(0, 0)) +
  scale_fill_gradientn("Z", colours = viridis::viridis(10)) +
  coord_fixed() + theme(line = element_blank(),
                        text = element_blank(),
                        title = element_blank(), 
                        panel.background = element_rect(fill = "white"),
                        legend.position = "none")

plot_gg(ggvolcano, multicore = TRUE, raytrace = TRUE, width = 7, height = 4, 
        scale = 300, windowsize = c(1000, 1000), zoom = 0.6, phi = 30, theta = 30)




ggvolcano <- InputExp %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = Exp_SG3)) +
  #  geom_contour(aes(x = x, y = y, z = Exp), breaks = c(0,  20, 100)) +
  scale_x_continuous("X", expand = c(0, 0)) +
  scale_y_continuous("Y", expand = c(0, 0)) +
  scale_fill_gradientn("Z", colours = viridis::viridis(10)) +
  coord_fixed() + theme(line = element_blank(),
                        text = element_blank(),
                        title = element_blank(), 
                        panel.background = element_rect(fill = "white"),
                        legend.position = "none")

plot_gg(ggvolcano, multicore = TRUE, raytrace = TRUE, width = 7, height = 4, 
        scale = 300, windowsize = c(1000, 1000), zoom = 0.6, phi = 30, theta = 30)





set.seed(1)
Lognorm_Sample <- rlnorm(100000)
Lognorm_Sample_density <- density(Lognorm_Sample)

Lognorm_Sample_density <- data.frame(x = Lognorm_Sample_density$x, y = Lognorm_Sample_density$y)
Lognorm_Sample_poly <- Lognorm_Sample_density[which(Lognorm_Sample_density$x>4),]
png(file=".\\Output\\Fig1_Pattern\\LognormDensity.png",width = 6,height = 6,units = "in",
    res = 600)
ggplot(Lognorm_Sample_density, aes(x = x, y = y)) +
  geom_area(data = Lognorm_Sample_poly, aes(x = x, y = y), fill = "darkblue", alpha = 0.5) +
  geom_line() +
  xlim(-1, 10) +
  theme_void()
dev.off()
