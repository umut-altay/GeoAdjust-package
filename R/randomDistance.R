# Generate random distances with respect to the selected jittering scale.
# type : a vector of location types : U for urban, R for rural.
# s = scaling factor (1, 3 ,5, 10).
randomDistance<- function(type = NULL, s = NULL){
  distance<- rep(0, length(type))
  for (i in 1:length(type)){
    if (type[[i]]=="U"){
      distance[[i]]=stats::runif(1, min = 0, max = 2*s)}
    else {
      if (stats::runif(1) < 0.01){
        distance[[i]]=stats::runif(1, min = 0, max = 10*s)
      } else{
        distance[[i]]=stats::runif(1, min = 0, max = 5*s)
      }
      list(rand.dist=distance)
    }}
  return(distance)
}
