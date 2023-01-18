#' Convert a set of coordinates in degrees into a new set of coordinates in kilometers.
#'
#' @param loc A two column matrix of coordinates (The first column is longitude and the second column is latitude).
#' @return A two column matrix of coordinates (The first column is easting and the second column is northing).
#' @examples
#' convertDegToKM(loc = NULL)
#' @export
#' @import sp
convertDegToKM = function(loc, crs = sp::CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs")){
  locLatLon = sp::SpatialPoints(loc,
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
  locKM = sp::spTransform(locLatLon,
                          crs)
  return(locKM@coords[,c(1,2)])
}
