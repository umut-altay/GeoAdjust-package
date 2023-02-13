#' Creates a grid of locations within the bounding box of the national borders of a country of interest.
#'
#' @param admin0 A SpatialPolygonsDataFrame object representing the national (admin0) level borders of the country
#' @param res An integer representing the resolution in kilometers.
#' @return A list. The first element of the list, "predRast", is the prediction raster. The second element of the list, "loc.pred", is a data frame containing the grid of coordinates of the cell centers (both in degrees and in kilometers) of the prediction raster.
#' @examples
#' \dontrun{
#' res is resolution in kilometers. Example : res = 5 means that the raster cells are created as 5 km by 5 km squares.
#' grid <- gridCountry(admin0 = admin0, res = res)
#' }
#' @export
#' @import spatialEco
gridCountry = function(admin0, res){

  # xx = seq(admin0@bbox[1,1], admin0@bbox[1,2], length.out = m)
  # yy = seq(admin0@bbox[2,1], admin0@bbox[2,2], length.out = n)
  # grid = cbind(rep(xx, each = length(yy)), rep(yy, length(xx)))
  #
  # #Convert prediction grid into a SpatialPointsDataFrame object
  # grid = data.frame(xCoor = grid[,1], yCoor = grid[, 2])
  # #
  # grid <- sp::SpatialPointsDataFrame(grid, data.frame(id=1:length(grid[,1]))) #we have 2500 points in the initial grid
  # grid@proj4string@projargs = admin0@proj4string@projargs
  # #
  # grid=spatialEco::erase.point(grid, admin0, inside = FALSE)
  #
  # loc.pred = cbind(grid@coords[ ,1], grid@coords[ ,2])
  # nPred = length(loc.pred[,1])
  # loc.pred = data.frame(long = loc.pred[,1], lat = loc.pred[,2], east = rep(NA, nPred), north = rep(NA, nPred))
  # loc.pred[,c("east", "north")] = convertDegToKM(loc.pred[,c("long", "lat")])
  # return(loc.pred)

  proj = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
  admin0_trnsfrmd = spTransform(admin0, proj)
  xmin = admin0_trnsfrmd@bbox[[1,1]]
  xmax = admin0_trnsfrmd@bbox[[1,2]]
  ymin = admin0_trnsfrmd@bbox[[2,1]]
  ymax = admin0_trnsfrmd@bbox[[2,2]]

  res=res

  predRast <- raster(xmn = xmin , xmx = xmax , ymn = ymin , ymx = ymax, resolution =res,
                     crs="+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs")

  idx = 1:ncell(predRast)
  loc.pred = xyFromCell(predRast, idx)
  loc.pred = data.frame(east = loc.pred[,1], north = loc.pred[,2])
  loc.pred[,c("long", "lat")] = convertKMToDeg(loc.pred[,c("east", "north")])
  return(list(loc.pred = loc.pred, predRast=predRast))
}

