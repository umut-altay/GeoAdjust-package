#' Converts a set of coordinates in degrees into a new set of coordinates in Mollweide coordinate system.
#'
#' @param loc A two column matrix of coordinates (The first column is longitude and the second column is latitude).
#' @return A two column matrix of coordinates in Mollweide (https://pubs.usgs.gov/pp/1395/report.pdf) coordinate system.
#' @examples
#' path1 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' locDegree <- cbind(surveyData$long, surveyData$lat)
#' locMoll <- convertDegToMollweide(loc = locDegree)
#' head(locMoll)
#' @export
convertDegToMollweide = function(loc){
  locLatLon = sp::SpatialPoints(loc,
                            proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  locMollweide = sp::spTransform(locLatLon,
                                 sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  return(locMollweide@coords[,c(1,2)])
}
