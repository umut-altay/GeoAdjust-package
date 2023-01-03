#' Convert a set of coordinates in kilometers into a new set of coordinates in degrees.
#'
#' @param loc A two column matrix of coordinates (the first column is easting and the second column is northing).
#' @return A two column matrix of coordinates (the first column is longitude and the second column is latitude).
#' @examples
#' convertKMToDeg(loc = NULL)
#' @export
#' @import sp
convertKMToDeg = function(loc) {
  locSP = SpatialPoints(loc, proj4string=CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))
  lonLatCoords = spTransform(locSP, CRS("+proj=longlat +datum=WGS84"))
  attr(lonLatCoords, "coords")
}
