#Generate random displacement angles (between 0 - 2pi)

#Argument: length: Number of observations.
randomAngle<- function(length = NULL){
  angle<- (stats::runif(length, min = 0, max = 2*pi))
  return(angle)
}
