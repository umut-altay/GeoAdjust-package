#' Converts a set of coordinates in degrees into a new set of coordinates in
#' The Universal Transverse Mercator (UTM) zone:37 coordinate system
#'
#'
#' @param loc A two column matrix of coordinates (The first column is longitude
#' and the second column is latitude).
#' @return A two column matrix of  coordinates in The Universal Transverse
#' Mercator (UTM) zone:37 (https://www.usgs.gov/faqs/what-does-term-utm-mean-
#' utm-better-or-more-accurate-latitudelongitude).
#' @examples
#' path1 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' locDegree <- cbind(surveyData$long, surveyData$lat)
#' locKM <- convertDegToKM(loc = locDegree)
#' head(locKM)
#' @export
convertDegToKM = function(loc){
  crs = sp::CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs")
  locLatLon = sp::SpatialPoints(loc,
                            proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  locKM = sp::spTransform(locLatLon,
                          crs)
  return(locKM@coords[,c(1,2)])
}
