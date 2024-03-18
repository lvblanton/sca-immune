getMaxMin <- function(x_value){
  x_val_fit <- ( deming.95ci.lines$slope * x_value) + deming.95ci.lines$intercept
  xFit.max <- max(x_val_fit)
  xFit.min <- min(x_val_fit)
  return(c(xFit.min, xFit.max))
}
