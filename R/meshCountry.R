#' Creates a constrained refined Delaunay triangulation mesh based on the country borders.
#'
#' @param admin0 A SpatialPolygonsDataFrame representing the country borders,
#' in UTM:zone 37 coordinate system.
#' @param max.edge A vector of two values. The first and the second elements of
#' the vector represent the largest allowed triangle lengths for the inner and outer mesh, respectively.
#' @param offset A value representing the extension distance for the inla.mesh.2d object
#' @return A constrained refined Delaunay triangulation mesh created based on the country borders.
#' @export
#' @examples
#' if(requireNamespace("INLA")){
#' path1 <- system.file("extdata", "adm0UTM37.rda", package = "GeoAdjust")
#' load(path1)
#' mesh.s <- meshCountry(admin0 = adm0UTM37, max.edge = c(25, 50), offset = -.08)
#' }
#' @export
meshCountry = function(admin0 = NULL,max.edge = NULL,offset = NULL){

  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }

  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }
  }

  mesh.s <- INLA::inla.mesh.2d(boundary = admin0,
                         offset=-.08,
                         cutoff=4,
                         max.edge = c(25, 50))
  return(mesh.s)
}









