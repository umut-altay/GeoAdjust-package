# Displace locations w.r.t distance and angle.
relocate <- function(east = NULL, north = NULL, angle = NULL, distance = NULL){
  locx=rep(0, length(north))
  locy=rep(0, length(north))
  for (i in 1:length(north)){
    (locy[i]=north[i]+((distance[i])*sin(angle[i])))&
      (locx[i]=east[i]+((distance[i])*cos(angle[i])))}
  results=data.frame(locx=locx, locy=locy)
  return(results)
}
