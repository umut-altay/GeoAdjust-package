#' Converts a set of coordinates in degrees into a new set of coordinates in The Universal Transverse Mercator (UTM) zone:37 coordinate system
#'
#'
#' @param loc A two column matrix of coordinates (The first column is longitude and the second column is latitude).
#' @param crs The corresponding coordinate reference system
#' @return A two column matrix of The Universal Transverse Mercator (UTM) zone:37 coordinate system (https://www.usgs.gov/faqs/what-does-term-utm-mean-utm-better-or-more-accurate-latitudelongitude).
#' @examples
#' data("clusterData")
#' loc = cbind(clusterData$long, clusterData$lat)
#' locKM <- convertDegToKM(loc = loc)
#' head(locKM)
#' @export
#' @import sp
convertDegToKM = function(loc, crs = sp::CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs")){
  locLatLon = sp::SpatialPoints(loc,
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
  locKM = sp::spTransform(locLatLon,
                          crs)
  return(locKM@coords[,c(1,2)])
}
