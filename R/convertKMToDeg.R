#' Converts a set of coordinates in The Universal Transverse Mercator (UTM) zone:37 coordinate system into a new set of coordinates in degrees.
#'
#' @param loc A two column matrix of coordinates in The Universal Transverse Mercator (UTM) zone:37 coordinate system (https://www.usgs.gov/faqs/what-does-term-utm-mean-utm-better-or-more-accurate-latitudelongitude).
#' @return A two column matrix of coordinates in degrees (the first column is longitude and the second column is latitude).
#' @examples
#' locDegree <- convertKMToDeg(loc = locKM)
#' @export
#' @import sp
convertKMToDeg = function(loc) {
  locSP = SpatialPoints(loc, proj4string=CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))
  lonLatCoords = spTransform(locSP, CRS("+proj=longlat +datum=WGS84"))
  attr(lonLatCoords, "coords")
}
