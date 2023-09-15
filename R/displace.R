displace = function(scale = NULL, locKM = NULL, urbanRural = NULL, AdminShapeFile = NULL, check1 = NULL, boundary = NULL){

  crs_km = sf::st_crs(locKM)
  crs_shapefile = sf::st_crs(AdminShapeFile)

  eastOriginal = sf::st_coordinates(locKM)[,1]
  northOriginal = sf::st_coordinates(locKM)[,2]

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
          newPoint = data.frame(east = newPoint_eastNorth[,1], north = newPoint_eastNorth[,2])

          newPointKM_SF = sf::st_as_sf(newPoint, coords=c("east","north"), crs = crs_km)

          newPointDegree_SF = sf::st_transform(newPointKM_SF, crs_shapefile)

          check2 <- sf::st_join(newPointDegree_SF, AdminShapeFile)

          if ((is.na(check2[,"NAME_2"][[1]]) == FALSE) & (check2[,"NAME_2"][[1]] == check1[,"NAME_2"][[i]])){
            break
          }else{next}
        }else{break}
      }
      newLocationSet[[i,1]] = newPoint_eastNorth[[1,1]]
      newLocationSet[[i,2]] = newPoint_eastNorth[[1,2]]
    }
    newLocationSet = sf::st_as_sf(newLocationSet, coords=c("east","north"), crs = crs_km)
    jitteredCoords[[j]] = newLocationSet
  }
  return(jitteredCoords)
}
