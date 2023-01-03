#' Convert a set of coordinates in degrees into a new set of coordinates in Mollweide coordinate system.
#'
#' @param loc A two column matrix of coordinates (The first column is longitude and the second column is latitude).
#' @return A two column matrix of coordinates in Mollweide coordinate system.
#' @examples
#' convertDegToMollweide(loc = NULL)
#' @export
#' @import sp
convertDegToMollweide = function(loc){
  locLatLon = SpatialPoints(loc,
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
  locMollweide = spTransform(locLatLon,
                             CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  return(locMollweide@coords[,c(1,2)])
}
