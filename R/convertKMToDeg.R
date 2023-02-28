#' Converts a set of coordinates in The Universal Transverse Mercator (UTM) zone:37 coordinate system into a new set of coordinates in degrees.
#'
#' @param loc A two column matrix of coordinates in The Universal Transverse Mercator (UTM) zone:37 coordinate system (https://www.usgs.gov/faqs/what-does-term-utm-mean-utm-better-or-more-accurate-latitudelongitude).
#' @return A two column matrix of coordinates in degrees (the first column is longitude and the second column is latitude).
#' @examples
#' path1 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' locKM <- cbind(surveyData$east, surveyData$north)
#' locDegree <- convertKMToDeg(loc = locKM)
#' head(locDegree)
#' @export
convertKMToDeg = function(loc) {
  locSP = sp::SpatialPoints(loc, proj4string=sp::CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))
  lonLatCoords = sp::spTransform(locSP, sp::CRS("+proj=longlat +datum=WGS84"))
  attr(lonLatCoords, "coords")
}
