#' Creates a Matern covariance matrix.
#'
#' @param dMat A distance matrix between the locations.
#' @param range Spatial range in kilometers.
#' @param stdDev The marginal variance.
#' @return Matern covariance matrix.
#' @examples
#' data("clusterData")
#' loc = cbind(clusterData$east, clusterData$north)
#' space.range = 114
#' space.sigma = 1
#' covMat <- covMatern(dMat = as.matrix(dist(loc)),
#' range = space.range, stdDev = space.sigma)
#' @export
#' @import INLA
covMatern = function(dMat = NULL, range = NULL, stdDev = NULL){
  Sig = INLA::inla.matern.cov(nu = 1,
                        kappa = sqrt(8*1)/range,
                        x = dMat,
                        corr = TRUE)
  Sig = stdDev^2*Sig
}
