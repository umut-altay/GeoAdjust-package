#' Creates a grid of locations within the bounding box of the national borders of a country of interest.
#'
#' @param admin0 A SpatialPolygonsDataFrame object representing the national (admin0) level borders of the country.
#' @param res A value representing the resolution in kilometers.
#' @return A list. The first element of the list, predRast, is the prediction raster. The second element of the list, loc.pred, is a data frame containing the grid of coordinates of the cell centers (both in degrees and in kilometers) of the prediction raster.
#' @examples
#' \donttest{
#' path1 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' grid <- gridCountry(admin0 = adm0, res = 5)
#' }
#' @export
gridCountry = function(admin0 = NULL, res = NULL){

  proj = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
  admin0_trnsfrmd = sp::spTransform(admin0, proj)
  xmin = admin0_trnsfrmd@bbox[[1,1]]
  xmax = admin0_trnsfrmd@bbox[[1,2]]
  ymin = admin0_trnsfrmd@bbox[[2,1]]
  ymax = admin0_trnsfrmd@bbox[[2,2]]

  predRast <- raster::raster(xmn = xmin , xmx = xmax , ymn = ymin , ymx = ymax, resolution =res,
                     crs="+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs")

  idx = 1:raster::ncell(predRast)
  loc.pred = raster::xyFromCell(predRast, idx)
  loc.pred = data.frame(east = loc.pred[,1], north = loc.pred[,2])
  loc.pred[,c("long", "lat")] = convertKMToDeg(loc.pred[,c("east", "north")])
  return(list(loc.pred = loc.pred, predRast=predRast))
}

