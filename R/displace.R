#' Randomly displaces (Jitters) a set of locations
#'
#' @param scale Scaling factor. Use 1 to apply default maximum jittering
#' distances of DHS. Values other than 1 will scale the default values.
#' @param locKM A data frame of the coordinates (in kilometers) of the
#' corresponding locations that are to be jittered. The first column is easting
#'  (in kilometers) and the second column is northing (in kilometers)
#' @param urbanRural A vector containing the urbanicity classification types of
#' the administrative areas that the corresponding locations are located within
#'   : U for urban areas, R for rural areas
#' @param AdminShapeFile A shape file containing the borders of the
#' administrative areas that will be respected (or not) while the locations are
#' being jittered
#' @param check1 A data frame containing the names of the administrative areas
#' that each location is initially located within. It is being used to compare
#' the initial administrative area of each location with the administrative
#' area that they land into after jittering. The data frame can be obtained
#' by sp::over(yourLocations, AdminShapeFile, returnList = FALSE)
#' @param boundary A logical constant (TRUE/FALSE). When TRUE is chosen,
#' administrative area borders are respected and the jittering is repeated until
#'  the jittered location lands into the same administrative area that it was
#'  initially located within. When FALSE is chosen the administrative area
#'  borders are not respected while jittering
#' @return A matrix containing the coordinates of the displaced locations
#' @examples
#' \dontrun{
#' data("clusterData")
#' locKM = cbind(clusterData$east, clusterData$north)
#' scale = 1
#' urbanRural = clusterData$urbanRural
#' # here we use admin2 borders. It can be different for different countries:
#' admin2 = rgdal::readOGR(dsn = "dataFiles/gadm40_NGA_shp",
#' layer = "gadm40_NGA_2")
#' #remove the lake (Water body) :
#' admin2 = admin2[-160,]
#' # add administrative area IDs
#' admin2@data[["OBJECTID"]] =1:774
#' loc = cbind(clusterData$long, clusterData$lat)
#' colnames(loc) = c("long", "lat")
#  loc = sp::SpatialPoints(loc, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"), bbox = NULL)
#' check1 <- sp::over(loc, admin2, returnList = FALSE)
#' locJittered <- displace(scale = scale, locKM = locKM,
#' urbanRural = urbanRural, AdminShapeFile = Admin2ShapeFile, check1 = check1,
#' boundary = TRUE)
#' }
#' @export
#' @import sp
displace = function(scale, locKM, urbanRural, AdminShapeFile, check1, boundary){
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
          newPoint_spatialPointsObject = sp::SpatialPoints(cbind(newPoint_eastNorth[,1], newPoint_eastNorth[,2]), proj4string = CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"), bbox = NULL)
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
