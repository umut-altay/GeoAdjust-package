displace = function(scale = NULL, locKM = NULL, urbanRural = NULL, AdminShapeFile = NULL, check1 = NULL, boundary = NULL){
  eastOriginal = locKM[, "east"]
  northOriginal = locKM[,"north"]
  nLoc = length(eastOriginal)
  jitteredCoords = list()
  for (j in 1:length(scale)){
    newLocationSet=data.frame(east = rep(NA, nLoc), north = rep(NA, nLoc))
    for (i in 1:nLoc){
      repeat{
        east = eastOriginal[i]; north = northOriginal[i]; angle = randomAngle(1)
        distance = randomDistance(type = urbanRural[i], s = scale[j])
        newPoint_eastNorth = relocate(east = east, north = north, angle = angle, distance = distance)
        if (boundary == "TRUE"){
          newPoint_spatialPointsObject = sp::SpatialPoints(cbind(newPoint_eastNorth[,1], newPoint_eastNorth[,2]), proj4string = sp::CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"), bbox = NULL)
          newPoint_longLat <- sp::spTransform(newPoint_spatialPointsObject, Admin2ShapeFile@proj4string)
          check2 <- sp::over(newPoint_longLat, Admin2ShapeFile, returnList = FALSE)
          if ((is.na(check2[,"NAME_2"][[1]]) == FALSE) & (check2[,"NAME_2"][[1]] == check1[,"NAME_2"][[i]])){
            break
          }else{next}
        }else{break}
      }
      newLocationSet[[i,1]] = newPoint_eastNorth[[1,1]]
      newLocationSet[[i,2]] = newPoint_eastNorth[[1,2]]
    }
    jitteredCoords[[j]] = newLocationSet
  }
  return(jitteredCoords)
}
