#' Creates a grid of locations within the bounding box of the national borders of a country of interest.
#'
#' @param admin0 An sf class MULTIPOLYGON representing the country borders.
#' @param res A value representing the resolution in kilometers.
#' @param target_crs A projection string representing the desired coordinate
#' reference system according to which the prediction grid will be constructed.
#' The measurement unit of the target_crs should be in kilometers.
#' @return A list. The first element of the list, predRast, is a SpatRaster object.
#' The second element of the list, loc.pred, is an sf class POINT object containing
#' coordinates of the cell centers (in target_crs) of the prediction raster.
#' @examples
#' path1 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' crs_KM = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
#' grid = gridCountry(admin0 = adm0, res = 5, target_crs = crs_KM)
#' @export
gridCountry = function(admin0 = NULL, res = NULL, target_crs = NULL){

  admin0_trnsfrmd = sf::st_transform(admin0, target_crs)

  xmin = sf::st_bbox(admin0_trnsfrmd)[[1]]
  xmax = sf::st_bbox(admin0_trnsfrmd)[[3]]
  ymin = sf::st_bbox(admin0_trnsfrmd)[[2]]
  ymax = sf::st_bbox(admin0_trnsfrmd)[[4]]

  predRast <- terra::rast(xmin = xmin , xmax = xmax , ymin = ymin , ymax = ymax)
  terra::res(predRast) = res
  terra::crs(predRast) = target_crs

  idx = 1:terra::ncell(predRast)
  loc.pred = terra::xyFromCell(predRast, idx)
  loc.pred = data.frame(east = loc.pred[,1], north = loc.pred[,2])
  loc.pred = sf::st_as_sf(loc.pred, coords = c("east", "north"), crs = target_crs)
  return(list(loc.pred = loc.pred, predRast=predRast))
}

