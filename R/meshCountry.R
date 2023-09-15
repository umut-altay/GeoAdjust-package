#' Creates a constrained refined Delaunay triangulation mesh based on the country borders.
#'
#' @param admin0 An sf class multipolygon representing the country borders.
#' @param max.edge A vector of two values. The first and the second elements of
#' the vector represent the largest allowed triangle lengths for the inner and
#' outer mesh, respectively.
#' @param cutoff The minimum allowed distance of the vertices to each other.
#' @param offset A value representing the extension distance for the mesh.
#' @param target_crs A projection string representing the desired coordinate
#' reference system according to which the mesh will be constructed. The
#' measurement unit of the target_crs should be in kilometers.
#' @return A constrained refined Delaunay triangulation mesh created based on
#' the country borders.
#' @export
#' @examples
#' path1 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' crs_KM = "+units=km +proj=utm +zone=37 +ellps=clrk80
#' +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
#' mesh.s <- meshCountry(admin0= adm0, max.edge = c(25, 50), offset = -.08,
#' cutoff=4, target_crs = crs_KM)
#' @export
meshCountry = function(admin0 = NULL,max.edge = NULL,cutoff = NULL, offset = NULL, target_crs = NULL){

  admin0_trnsfrmd = sf::st_transform(admin0, target_crs)

  mesh.s <- fmesher::fm_mesh_2d(boundary = admin0_trnsfrmd,
                         offset=offset,
                         cutoff=cutoff,
                         max.edge = max.edge)
  return(mesh.s)
}









