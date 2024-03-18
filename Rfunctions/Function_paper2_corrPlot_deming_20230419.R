library(ggplot2) 
library(deming)

getMaxMin <- function(x_value, mydeming_lines){
  x_val_fit <- ( mydeming_lines$slope * x_value) + mydeming_lines$intercept
  xFit.max <- max(x_val_fit)
  xFit.min <- min(x_val_fit)
  return(c(xFit.min, xFit.max))
}

default_pointsCol <- "#40404040"
orange_pointsCol <- "#D95F0280"
purple_pointsCol <- "#6e02d980"

paper2_corrPlot_deming <- function(myData,myTitle = "corrPlot_deming.pdf", myWidth=1.5, myHeight=1.5, myXlab = "X variable", myYlab = "Y variable", myPlotTitle = NULL, 
                                   negIdentity=FALSE, myYlim = NULL, myXlim=NULL, printRsq = FALSE, pointsCol = default_pointsCol, demingBehind = FALSE){
  print(summary(myData))
  
  if(is.null(myYlim) & is.null(myXlim)){
    myLim <- makeSymmetric_X_Y(myData$x, myData$y)
    myData_min_x <- myLim[1] - 0.2
    myData_max_x <- myLim[2] + 0.2
    myData_min_y <- myLim[1] - 0.2
    myData_max_y <- myLim[2] + 0.2
  } else {
    myData_min_x <- myXlim[1] - 0.2
    myData_max_x <- myXlim[2] + 0.2
    myData_min_y <- myYlim[1] - 0.2
    myData_max_y <- myYlim[2] + 0.2
  }

  my.deming <- deming(formula = myData$y ~ myData$x) 
  deming.slope <- my.deming$coefficients[2]
  deming.int <- my.deming$coefficients[1]
  
  #get the 95% CI for the slope and intercept
  deming.slope_seq <- seq(my.deming$ci[2,1], my.deming$ci[2,2], 0.001)
  deming.int_seq <- seq(my.deming$ci[1,1], my.deming$ci[1,2], 0.0001)
  
  x_range <- seq(myData_min_x - 0.2, myData_max_x + 0.2, 0.05)
  deming.ribbon_table <- data.frame("x" = x_range, "fit" = ((deming.slope * x_range) + deming.int), "y.min" = numeric(length = length(x_range)), "y.max" = numeric(length=length(x_range)))
  deming.95ci.lines <- data.frame("intercept" = rep(deming.int_seq, each = length(deming.slope_seq)), "slope" = deming.slope_seq)
  
  #for each x_range, calculate all of the deming.95ci.line fits and get max and min
  deming.ribbon_table[,c("y.min","y.max")] <- t(sapply(x_range, getMaxMin, mydeming_lines = deming.95ci.lines))
  
  myCor_pear <- cor.test(myData$x, myData$y, method=c("pearson"))
  print(myCor_pear)
  
  if(demingBehind){
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max), fill = default_pointsCol) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0) +
      geom_abline(slope = deming.slope, intercept = deming.int, col="black", size=0.5) + 
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle
      ) 
  } else {
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max)) +
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle
      ) 
  }
  
  
  if(printRsq){
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #bottom left  
          y = myData_min_y + 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 0,
          size=2,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2 ,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2,
          col="black"
        ) 
    }
  } else {
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_max_x - 0.1, #top right
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 1,
          vjust = 1,
          size=2,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2,
          col="black"
        ) 
    }
  }
  

  pdf(file = myTitle, width = myWidth, height = myHeight)
  print(p)
  dev.off()
}

paper2_corrPlot_deming_Large <- function(myData,myTitle = "corrPlot_deming.pdf", myWidth=1.5, myHeight=1.5, myXlab = "X variable", myYlab = "Y variable", myPlotTitle = NULL, 
                                   negIdentity=FALSE, myYlim = NULL, myXlim=NULL, printRsq = FALSE, pointsCol = default_pointsCol, demingBehind = FALSE){
  print(summary(myData))
  
  if(is.null(myYlim) & is.null(myXlim)){
    myLim <- makeSymmetric_X_Y(myData$x, myData$y)
    myData_min_x <- myLim[1] - 0.2
    myData_max_x <- myLim[2] + 0.2
    myData_min_y <- myLim[1] - 0.2
    myData_max_y <- myLim[2] + 0.2
  } else {
    myData_min_x <- myXlim[1] - 0.2
    myData_max_x <- myXlim[2] + 0.2
    myData_min_y <- myYlim[1] - 0.2
    myData_max_y <- myYlim[2] + 0.2
  }
  
  my.deming <- deming(formula = myData$y ~ myData$x) 
  deming.slope <- my.deming$coefficients[2]
  deming.int <- my.deming$coefficients[1]
  
  #get the 95% CI for the slope and intercept
  deming.slope_seq <- seq(my.deming$ci[2,1], my.deming$ci[2,2], 0.01)
  deming.int_seq <- seq(my.deming$ci[1,1], my.deming$ci[1,2], 0.01)
  
  x_range <- seq(myData_min_x - 0.2, myData_max_x + 0.2, 0.1)
  deming.ribbon_table <- data.frame("x" = x_range, "fit" = ((deming.slope * x_range) + deming.int), "y.min" = numeric(length = length(x_range)), "y.max" = numeric(length=length(x_range)))
  deming.95ci.lines <- data.frame("intercept" = rep(deming.int_seq, each = length(deming.slope_seq)), "slope" = deming.slope_seq)
  
  #for each x_range, calculate all of the deming.95ci.line fits and get max and min
  deming.ribbon_table[,c("y.min","y.max")] <- t(sapply(x_range, getMaxMin, mydeming_lines = deming.95ci.lines))
  
  myCor_pear <- cor.test(myData$x, myData$y, method=c("pearson"))
  print(myCor_pear)
  
  if(demingBehind){
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max), fill = default_pointsCol) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0) +
      geom_abline(slope = deming.slope, intercept = deming.int, col="black", size=0.5) + 
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle
      ) 
  } else {
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max)) +
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle
      ) 
  }
  
  
  if(printRsq){
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #bottom left  
          y = myData_min_y + 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 0,
          size=2,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2 ,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2,
          col="black"
        ) 
    }
  } else {
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_max_x - 0.1, #top right
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 1,
          vjust = 1,
          size=2,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2,
          col="black"
        ) 
    }
  }
  
  
  pdf(file = myTitle, width = myWidth, height = myHeight)
  print(p)
  dev.off()
}

paper2_corrPlot_deming_Large_unequal <- function(myData,myTitle = "corrPlot_deming.pdf", myWidth=1.5, myHeight=1.5, myXlab = "X variable", myYlab = "Y variable", myPlotTitle = NULL, 
                                         negIdentity=FALSE, myYlim = NULL, myXlim=NULL, printRsq = FALSE, pointsCol = default_pointsCol, demingBehind = FALSE){
  print(summary(myData))
  
  if(is.null(myYlim) & is.null(myXlim)){
    myLim <- makeSymmetric_X_Y(myData$x, myData$y)
    myData_min_x <- myLim[1] - 0.2
    myData_max_x <- myLim[2] + 0.2
    myData_min_y <- myLim[1] - 0.2
    myData_max_y <- myLim[2] + 0.2
  } else {
    myData_min_x <- myXlim[1] - 0.2
    myData_max_x <- myXlim[2] + 0.2
    myData_min_y <- myYlim[1] - 0.2
    myData_max_y <- myYlim[2] + 0.2
  }
  
  my.deming <- deming(formula = myData$y ~ myData$x) 
  deming.slope <- my.deming$coefficients[2]
  deming.int <- my.deming$coefficients[1]
  
  #get the 95% CI for the slope and intercept
  deming.slope_seq <- seq(my.deming$ci[2,1], my.deming$ci[2,2], 0.01)
  deming.int_seq <- seq(my.deming$ci[1,1], my.deming$ci[1,2], 0.01)
  
  x_range <- seq(myData_min_x - 0.2, myData_max_x + 0.2, 0.1)
  deming.ribbon_table <- data.frame("x" = x_range, "fit" = ((deming.slope * x_range) + deming.int), "y.min" = numeric(length = length(x_range)), "y.max" = numeric(length=length(x_range)))
  deming.95ci.lines <- data.frame("intercept" = rep(deming.int_seq, each = length(deming.slope_seq)), "slope" = deming.slope_seq)
  
  #for each x_range, calculate all of the deming.95ci.line fits and get max and min
  deming.ribbon_table[,c("y.min","y.max")] <- t(sapply(x_range, getMaxMin, mydeming_lines = deming.95ci.lines))
  
  myCor_pear <- cor.test(myData$x, myData$y, method=c("pearson"))
  print(myCor_pear)
  
  if(demingBehind){
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max), fill = default_pointsCol) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0) +
      geom_abline(slope = deming.slope, intercept = deming.int, col="black", size=0.5) + 
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25)
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle
      ) 
  } else {
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max)) +
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle
      ) 
  }
  
  
  if(printRsq){
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #bottom left  
          y = myData_min_y + 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 0,
          size=2,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2 ,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2,
          col="black"
        ) 
    }
  } else {
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_max_x - 0.1, #top right
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 1,
          vjust = 1,
          size=2,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\nP = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2,
          col="black"
        ) 
    }
  }
  
  
  pdf(file = myTitle, width = myWidth, height = myHeight)
  print(p)
  dev.off()
}


paper2_corrPlot_deming_err <- function(myData,myTitle = "corrPlot_deming_err.pdf", myWidth=1.5, myHeight=1.5, myXlab = "X variable", myYlab = "Y variable", myPlotTitle = NULL, 
                                   negIdentity=FALSE, myYlim = NULL, myXlim=NULL, printRsq = FALSE, pointsCol = default_pointsCol, demingBehind = FALSE){
  print(summary(myData))
  
  if(is.null(myYlim) & is.null(myXlim)){
    myLim <- makeSymmetric_X_Y(myData$x, myData$y)
    myData_min_x <- myLim[1] - 0.2
    myData_max_x <- myLim[2] + 0.2
    myData_min_y <- myLim[1] - 0.2
    myData_max_y <- myLim[2] + 0.2
  } else {
    myData_min_x <- myXlim[1] - 0.2
    myData_max_x <- myXlim[2] + 0.2
    myData_min_y <- myYlim[1] - 0.2
    myData_max_y <- myYlim[2] + 0.2
  }
  
  my.deming <- deming(formula = myData$y ~ myData$x, xstd = myData$x_err, ystd = myData$y_err) 
  deming.slope <- my.deming$coefficients[2]
  deming.int <- my.deming$coefficients[1]
  
  #get the 95% CI for the slope and intercept
  deming.slope_seq <- seq(my.deming$ci[2,1], my.deming$ci[2,2], 0.01)
  deming.int_seq <- seq(my.deming$ci[1,1], my.deming$ci[1,2], 0.001)
  
  x_range <- seq(myData_min_x - 0.2, myData_max_x + 0.2, 0.05)
  deming.ribbon_table <- data.frame("x" = x_range, "fit" = ((deming.slope * x_range) + deming.int), "y.min" = numeric(length = length(x_range)), "y.max" = numeric(length=length(x_range)))
  deming.95ci.lines <- data.frame("intercept" = rep(deming.int_seq, each = length(deming.slope_seq)), "slope" = deming.slope_seq)
  
  #for each x_range, calculate all of the deming.95ci.line fits and get max and min
  deming.ribbon_table[,c("y.min","y.max")] <- t(sapply(x_range, getMaxMin, mydeming_lines = deming.95ci.lines))
  
  myCor_pear <- cor.test(myData$x, myData$y, method=c("pearson"))
  print(myCor_pear)
  
  if(demingBehind){
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", linewidth=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", linewidth=0.25) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max), fill = default_pointsCol, alpha = 0.3) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0, size = 2) +
      geom_abline(slope = deming.slope, intercept = deming.int, col="black", size=0.5) + 
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 7) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(linewidth = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle,
        subtitle = paste0("n = ", dim(myData)[1], " genes")
      ) 
  } else {
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0, size = 2) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max), alpha = 0.3) +
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 7) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle,
        subtitle = paste0("n = ", dim(myData)[1], " genes")
      ) 
  }
  
  
  if(printRsq){
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #bottom left  
          y = myData_min_y + 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 0,
          size=2.5,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2 ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    }
  } else {
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_max_x - 0.1, #top right
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 1,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    }
  } 
  
  
  pdf(file = myTitle, width = myWidth, height = myHeight)
  print(p)
  dev.off()
}

paper2_corrPlot_deming_text_err <- function(myData,myTitle = "corrPlot_deming_err.pdf", myWidth=1.5, myHeight=1.5, myXlab = "X variable", myYlab = "Y variable", myPlotTitle = NULL, 
                                       negIdentity=FALSE, myYlim = NULL, myXlim=NULL, printRsq = FALSE, pointsCol = default_pointsCol, demingBehind = FALSE){
  print(summary(myData))
  
  if(is.null(myYlim) & is.null(myXlim)){
    myLim <- makeSymmetric_X_Y(myData$x, myData$y)
    myData_min_x <- myLim[1] - 0.2
    myData_max_x <- myLim[2] + 0.2
    myData_min_y <- myLim[1] - 0.2
    myData_max_y <- myLim[2] + 0.2
  } else {
    myData_min_x <- myXlim[1] - 0.2
    myData_max_x <- myXlim[2] + 0.2
    myData_min_y <- myYlim[1] - 0.2
    myData_max_y <- myYlim[2] + 0.2
  }
  
  my.deming <- deming(formula = myData$y ~ myData$x, xstd = myData$x_err, ystd = myData$y_err) 
  deming.slope <- my.deming$coefficients[2]
  deming.int <- my.deming$coefficients[1]
  
  #get the 95% CI for the slope and intercept
  deming.slope_seq <- seq(my.deming$ci[2,1], my.deming$ci[2,2], 0.01)
  deming.int_seq <- seq(my.deming$ci[1,1], my.deming$ci[1,2], 0.001)
  
  x_range <- seq(myData_min_x - 0.2, myData_max_x + 0.2, 0.05)
  deming.ribbon_table <- data.frame("x" = x_range, "fit" = ((deming.slope * x_range) + deming.int), "y.min" = numeric(length = length(x_range)), "y.max" = numeric(length=length(x_range)))
  deming.95ci.lines <- data.frame("intercept" = rep(deming.int_seq, each = length(deming.slope_seq)), "slope" = deming.slope_seq)
  
  #for each x_range, calculate all of the deming.95ci.line fits and get max and min
  deming.ribbon_table[,c("y.min","y.max")] <- t(sapply(x_range, getMaxMin, mydeming_lines = deming.95ci.lines))
  
  myCor_pear <- cor.test(myData$x, myData$y, method=c("pearson"))
  print(myCor_pear)
  
  if(demingBehind){
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", linewidth=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", linewidth=0.25) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max), fill = default_pointsCol, alpha = 0.3) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0, size = 2) +
      geom_abline(slope = deming.slope, intercept = deming.int, col="black", size=0.5) + 
      geom_text_repel(data = myData, aes(x = x, y=y, label = Gene), size = 1.75, fontface= "italic",max.overlaps = 20) +
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 7) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle,
        subtitle = paste0("n = ", dim(myData)[1], " genes")
      ) 
  } else {
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0, size = 2) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max), alpha = 0.3) +
      geom_text_repel(data = myData, aes(x=x, y = y, label = Gene), size = 1.5, max.overlaps = 45, fontface= "italic") +
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 7) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle,
        subtitle = paste0("n = ", dim(myData)[1], " genes")
      ) 
  }
  
  
  if(printRsq){
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #bottom left  
          y = myData_min_y + 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 0,
          size=2.5,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2 ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    }
  } else {
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_max_x - 0.1, #top right
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 1,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    }
  } 
  
  
  pdf(file = myTitle, width = myWidth, height = myHeight)
  print(p)
  dev.off()
}

paper2_corrPlot_deming_err_noCI <- function(myData,myTitle = "corrPlot_deming_err.pdf", myWidth=1.5, myHeight=1.5, myXlab = "X variable", myYlab = "Y variable", myPlotTitle = NULL, 
                                       negIdentity=FALSE, myYlim = NULL, myXlim=NULL, printRsq = FALSE, pointsCol = default_pointsCol, demingBehind = FALSE){
  print(summary(myData))
  
  if(is.null(myYlim) & is.null(myXlim)){
    myLim <- makeSymmetric_X_Y(myData$x, myData$y)
    myData_min_x <- myLim[1] - 0.2
    myData_max_x <- myLim[2] + 0.2
    myData_min_y <- myLim[1] - 0.2
    myData_max_y <- myLim[2] + 0.2
  } else {
    myData_min_x <- myXlim[1] - 0.2
    myData_max_x <- myXlim[2] + 0.2
    myData_min_y <- myYlim[1] - 0.2
    myData_max_y <- myYlim[2] + 0.2
  }
  
  my.deming <- deming(formula = myData$y ~ myData$x, xstd = myData$x_err, ystd = myData$y_err) 
  deming.slope <- my.deming$coefficients[2]
  deming.int <- my.deming$coefficients[1]
  
  #get the 95% CI for the slope and intercept
  deming.slope_seq <- seq(my.deming$ci[2,1], my.deming$ci[2,2], 0.01)
  deming.int_seq <- seq(my.deming$ci[1,1], my.deming$ci[1,2], 0.001)
  
  myCor_pear <- cor.test(myData$x, myData$y, method=c("pearson"))
  print(myCor_pear)
  
  if(demingBehind){
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0, size = 2) +
      geom_abline(slope = deming.slope, intercept = deming.int, col="black", size=0.5) + 
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 7) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle,
        subtitle = paste0("n = ", dim(myData)[1], " genes")
      ) 
  } else {
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", size=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", size=0.25) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0, size = 2) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max), alpha = 0.3) +
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 7) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(size = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle,
        subtitle = paste0("n = ", dim(myData)[1], " genes")
      ) 
  }
  
  
  if(printRsq){
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #bottom left  
          y = myData_min_y + 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 0,
          size=2.5,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2 ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    }
  } else {
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_max_x - 0.1, #top right
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 1,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    }
  } 
  
  
  pdf(file = myTitle, width = myWidth, height = myHeight)
  print(p)
  dev.off()
}

paper2_corrPlot_deming_text_err_noCI <- function(myData,myTitle = "corrPlot_deming_err.pdf", myWidth=1.5, myHeight=1.5, myXlab = "X variable", myYlab = "Y variable", myPlotTitle = NULL, 
                                            negIdentity=FALSE, myYlim = NULL, myXlim=NULL, printRsq = FALSE, pointsCol = default_pointsCol, demingBehind = FALSE){
  print(summary(myData))
  
  if(is.null(myYlim) & is.null(myXlim)){
    myLim <- makeSymmetric_X_Y(myData$x, myData$y)
    myData_min_x <- myLim[1] - 0.2
    myData_max_x <- myLim[2] + 0.2
    myData_min_y <- myLim[1] - 0.2
    myData_max_y <- myLim[2] + 0.2
  } else {
    myData_min_x <- myXlim[1] - 0.2
    myData_max_x <- myXlim[2] + 0.2
    myData_min_y <- myYlim[1] - 0.2
    myData_max_y <- myYlim[2] + 0.2
  }
  
  my.deming <- deming(formula = myData$y ~ myData$x, xstd = myData$x_err, ystd = myData$y_err) 
  deming.slope <- my.deming$coefficients[2]
  deming.int <- my.deming$coefficients[1]
  
  #get the 95% CI for the slope and intercept
  deming.slope_seq <- seq(my.deming$ci[2,1], my.deming$ci[2,2], 0.01)
  deming.int_seq <- seq(my.deming$ci[1,1], my.deming$ci[1,2], 0.001)
  
  myCor_pear <- cor.test(myData$x, myData$y, method=c("pearson"))
  print(myCor_pear)
  
  if(demingBehind){
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", linewidth=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", linewidth=0.25) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0, size = 2) +
      geom_abline(slope = deming.slope, intercept = deming.int, col="black", size=0.5) + 
      geom_text_repel(data = myData, aes(x = x, y=y, label = Gene), size = 1.75, fontface= "italic",max.overlaps = 20) +
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 7) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(linewidth = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle,
        subtitle = paste0("n = ", dim(myData)[1], " genes")
      ) 
  } else {
    p <- ggplot() +
      geom_hline(yintercept =0,col="darkgrey", linewidth=0.25) +
      geom_vline(xintercept = 0,  col="darkgrey", linewidth=0.25) +
      geom_point(data=myData, mapping=aes(x=x, y=y),col=pointsCol, stroke=0, size = 2) +
      geom_ribbon(data = deming.ribbon_table, mapping = aes(x=x, ymin =y.min, ymax=y.max), alpha = 0.3) +
      geom_text_repel(data = myData, aes(x=x, y = y, label = Gene), size = 1.5, max.overlaps = 45, fontface= "italic") +
      coord_cartesian(xlim = c(myData_min_x, myData_max_x), ylim =c(myData_min_y, myData_max_y) ) +
      theme_classic(base_size = 7) +
      theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.line = element_line(linewidth = 0.25),
        aspect.ratio=1
      )  +  labs(
        x= myXlab,
        y= myYlab,
        title= myPlotTitle,
        subtitle = paste0("n = ", dim(myData)[1], " genes")
      ) 
  }
  
  
  if(printRsq){
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #bottom left  
          y = myData_min_y + 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 0,
          size=2.5,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('R2 = ', formatC(myCor_pear$estimate ^ 2 ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    }
  } else {
    if (negIdentity){
      p <- p + 
        geom_abline(slope=-1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_max_x - 0.1, #top right
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 1,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    } else {
      p <- p + 
        geom_abline(slope=1, intercept=0, col="#02a3d9", size=0.25, lty=2) +
        annotate(
          geom = "text",
          x = myData_min_x + 0.1, #top left
          y = myData_max_y - 0.1,
          label = paste0('r = ', formatC(myCor_pear$estimate ,digits = 2,format = "f"), "\np = ", formatC(myCor_pear$p.value, digits=2,format = "e")),
          hjust = 0,
          vjust = 1,
          size=2.5,
          col="black"
        ) 
    }
  } 
  
  
  pdf(file = myTitle, width = myWidth, height = myHeight)
  print(p)
  dev.off()
}