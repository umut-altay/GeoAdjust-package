#' Creates a constrained refined Delaunay triangulation mesh based on the country borders.
#'
#' @param admin0 A SpatialPolygonsDataFrame representing the country borders.
#' @param max.edge A vector of two values. The first and the second elements of
#' the vector represent the largest allowed triangle lengths for the inner and outer mesh, respectively.
#' @param offset offset value representing the extension distance for the inla.mesh.2d object
#' @return A constrained refined Delaunay triangulation mesh created based on the country borders.
#' @export
#' @import INLA
#' @examples
#' \dontrun{
#' mesh.s <- meshCountry(admin0 = admin0, max.edge = c(25, 50), offset = -.08)
#' }
#' @export
meshCountry = function(admin0,max.edge,offset){
  mesh.s <- INLA::inla.mesh.2d(boundary = admin0,
                         offset=-.08,
                         cutoff=4,
                         max.edge = c(25, 50))
  return(mesh.s)
}









