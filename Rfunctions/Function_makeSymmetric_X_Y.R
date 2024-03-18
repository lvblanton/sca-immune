makeSymmetric_X_Y <- function(x_vector, y_vector){
  myMin.x <- abs(round(min(na.omit(x_vector)), digits = 1))
  myMax.x <- abs(round(max(na.omit(x_vector)),digits=1))
  myMin.y <- abs(round(min(na.omit(y_vector)), digits = 1))
  myMax.y <- abs(round(max(na.omit(y_vector)),digits=1))
  
  #get max of X and Y
  myMax_all <- max(myMin.x, myMax.x, myMin.y, myMax.y)
  myLim <- c(-myMax_all, myMax_all)
  return(myLim)
}