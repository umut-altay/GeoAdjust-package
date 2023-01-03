#' Create a set of prediction locations, chosen randomly from a prediction grid
#'
#' @param admin0 A SpatialPolygonsdataframe object reresenting the country (admin0) borders
#' @param m A value representing the length of the segments that the longitude range of the country will be divided into
#' @param n A value representing the length of the segments that the latitude range of the country will be divided into
#' @param nPred The number of prediction points which are to be randomly choosen from the created prediction grid
#' @return A data frame containing the coordinates (in degrees and in kilometers) of a set of prediction points
#' @examples
#' gridCountry(admin0 = NULL, m = NULL, n = NULL, nPred = NULL)
#' @export
#' @import spatialEco
gridCountry = function(admin0, m, n, nPred){

  xx = seq(admin0@bbox[1,1], admin0@bbox[1,2], length.out = m)
  yy = seq(admin0@bbox[2,1], admin0@bbox[2,2], length.out = n)
  grid = cbind(rep(xx, each = length(yy)), rep(yy, length(xx)))
  #
  # #Convert prediction grid into a SpatialPointsDataFrame object
  grid = data.frame(xCoor = grid[,1], yCoor = grid[, 2])
  #
  grid <- SpatialPointsDataFrame(grid, data.frame(id=1:2500)) #we have 2500 points in the initial grid
  grid@proj4string@projargs = admin0@proj4string@projargs
  #
  grid=erase.point(grid, admin0, inside = FALSE)

  loc.pred = cbind(grid@coords[ ,1], grid@coords[ ,2])
  idx = sample.int(dim(loc.pred)[1], size = nPred)
  loc.pred = loc.pred[idx,]
  loc.pred = data.frame(long = loc.pred[,1], lat = loc.pred[,2], east = rep(NA, nPred), north = rep(NA, nPred))
  loc.pred[,c("east", "north")] = convertDegToKM(loc.pred[,c("long", "lat")])
  return(loc.pred)
}

