
library(ggplot2)
library(scales)
library(ggpubr)
library(plot3D)

library(plotly)


Input_Folder <- ".\\..\\..\\Data\\3D_Pattern\\Simulation_3D_Example\\"
for(Pattern_index in 1:3){
  Input_File <- paste0("Pattern_",Pattern_index,".csv")
  Input_Data <- read.csv(paste0(Input_Folder, Input_File))
  
  SEED = Pattern_index 
  N = 225
  numSlices = 10
  width = 3.5
  SubPattern = Pattern_index
  set.seed(SEED)
  SimData <- data.frame()
  for(i in 1:numSlices){
    Coords <- spatstat.random::rpoispp(N, win = spatstat.geom::owin(c(0, 1), c(0, 1)))
    SimData <- rbind(SimData, data.frame(x = Coords$x * (sqrt(N) - 1), y = Coords$y * (sqrt(N) - 1), z = i))
  }
  coordDF <- SimData
  
  set.seed(SEED)
  if(SubPattern == "1"){
    Phi_const <- pi * round(runif(1, 0, 2)) - pi/2
    Phi <- pi*runif(20) + Phi_const
    Theta_const <- 0.5 * pi * round(runif(1, 0, 3))
    Theta <- 0.5 * pi*runif(20) + Theta_const
  }
  if(SubPattern == "2"){
    Phi_const <- pi * round(runif(1, 0, 2)) - pi/2
    Phi <- pi*runif(20) + Phi_const
    Theta_const <- pi * round(runif(1, 0, 2))
    Theta <- pi * runif(20) + Theta_const
  }                
  if(SubPattern == "3"){
    Phi <- 2 * pi * runif(20)
    Theta <- 2 * pi * runif(20)
  }
  Center_pts <- cbind(cumsum(2 * cos(Phi) * sin(Theta)) + sqrt(N)/2,
                      cumsum(2 * cos(Phi) * cos(Theta)) + sqrt(N)/2,
                      cumsum(2 * sin(Phi)) + numSlices/2)
  
  #plot(cumsum(2*sin(Theta)), cumsum(2*cos(Theta)),type = "b", xlim = c(-10,10), ylim =c(-10,10) )
  
  spikedCoords <- sapply(1:length(Phi), function(Center_pt){
    # get distance from every cell to the center of the coordinate set; assuming hot spot is centered here
    distFromCenter <- sp::spDists(
      x = coordDF %>% dplyr::select(x, y, z) %>% as.matrix(),
      y = matrix(Center_pts[Center_pt, ], ncol = 3),
      longlat = FALSE
    ) %>% as.vector()
    
    # identify coordinates in hot spot
    spikedCoord <- coordDF[distFromCenter <= (width-1),] %>% row.names() %>% as.numeric()
  })
  spikedCoords <- unique(unlist(spikedCoords))
  
  
  
  
  Input_Data_Plot <- Input_Data[,1:4]
  # trim extreme values
  Input_Data_Plot[Input_Data_Plot[,4]>quantile(Input_Data_Plot[,4],0.95),4] <- quantile(Input_Data_Plot[,4],0.95)
  Input_Data_PlotSub <- Input_Data_Plot[spikedCoords,]
  Center_pts <- as.data.frame(Center_pts)
  colnames(Center_pts) <- c("x", "y", "z")
  Center_pts <- Center_pts[which(apply(Center_pts, 1, max)<sqrt(N) & apply(Center_pts, 1, min)>0),]
  my_palette <- c("#FEFEE333",
                  "#F3F9D533",
                  "#004533")

  Fig <- plot_ly() %>% 
    add_trace(data = Input_Data_Plot, 
              x = ~x, y = ~y, z = ~z, 
              marker = list(size = 5),
              mode = "markers", 
              type = "scatter3d", 
              opacity = 0.3,
              color = ~SVG_1, 
              colors = my_palette) %>%
    add_trace(
      data = Input_Data_PlotSub,
      x = ~x,
      y = ~y,
      z = ~z, 
      marker = list(size = 5),
      mode = "markers", 
      type = "scatter3d", 
      opacity = 1.0,
      color = ~SVG_1,
      colors = my_palette,
      showlegend = FALSE)%>%
    add_trace(
      data = Center_pts,
      x = ~x,
      y = ~y,
      z = ~z, 
      type="scatter3d", 
      opacity = 0.3,
      mode="lines",
      line = list(color = "#1E8654", width = 5),
      showlegend = FALSE)%>%
    layout(showlegend = FALSE,
           margin = list(
             l = 20,
             r = 1,
             b = 20,
             t = 1,
             pad = 1
           ),
           scene = list(camera = list(eye = list(x = -1.25, y = 1.25, z = 0.75)),
                        xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
                        yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
                        zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>% 
    hide_colorbar()
  
  
  
  htmlwidgets::saveWidget(as_widget(Fig), 
                          paste0(".\\..\\..\\Output\\Power_Analysis\\Manuscript 3D\\Manuscript_3D_Pattern_",Pattern_index,".html"))
}









rm(list=ls())






InputData <- read.csv(".\\..\\..\\Data\\3D_Pattern\\Simulation_3D_Example\\Discrete_900by10_width2_qt88_Noise0_pw1.csv")

spikedIndex <- read.csv(".\\..\\..\\Data\\3D_Pattern\\Simulation_3D_Example\\Discrete_Index_900by10_width2_qt88_Noise0_pw1.csv")

Input_Data_PlotSub <- InputData[spikedIndex$Indx,]



Createcolorbar <- function(InputValues){
  Cont_part1 <- ramp.col(n = 100, c("#FEFEE333", 
                                    "#004533"))
  Cont_part2 <- ramp.col(n = floor(30 / (max(quantile(InputValues)[4] - min(InputValues), 1)) * (max(InputValues) - min(InputValues))) - 100, 
                         c("#004533", "#004533"), alpha = 1.0)
  Cont_part1_trans <- stringr::str_to_upper(as.hexmode(floor(seq(0.6, 1.0, length.out =100)*256)))
  Cont_part1_rev <- sapply(1:length(Cont_part1), function(i){
    paste0(substr(Cont_part1[i], 1, 7), Cont_part1_trans[i])
  })
  return(c(Cont_part1_rev[-100], Cont_part2))
}

# 3d real data 1
my_palette1 <- Createcolorbar(InputData[,"SVG_1"])
Fig1 <- plot_ly() %>% 
  add_trace(data = InputData, 
            x = ~x, y = ~y, z = ~z, 
            marker = list(size = 3),
            mode = "markers", 
            opacity = 0.3,
            type = "scatter3d", 
            color = ~SVG_1, 
            colors = my_palette1) %>%
  add_trace(
    data = Input_Data_PlotSub,
    x = ~x,
    y = ~y,
    z = ~z, 
    marker = list(size = 3),
    mode = "markers", 
    type = "scatter3d", 
    opacity = 1.0,
    color = ~SVG_1,
    colors = my_palette1,
    showlegend = FALSE) %>%
  layout(showlegend = FALSE,
         margin = list(
           l = 20,
           r = 1,
           b = 20,
           t = 1,
           pad = 1
         ),
         scene = list(camera = list(center = list(x = 0, y = 0, z = 0),
                                    eye = list(x = 0.75, y = -2.0, z = 2.0)),
                      aspectratio = list(x = 1, y = 1, z = 0.5),
                      xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
                      yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
                      zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>% 
  hide_colorbar()

htmlwidgets::saveWidget(as_widget(Fig1), ".\\..\\..\\Output\\Power_Analysis\\Manuscript 3D\\Manuscript_3D_DiscretePattern.html")
