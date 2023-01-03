#' Create a Matern covariance matrix.
#'
#' @param dMat Distance matrix.
#' @param range Spatial range in kilometers.
#' @param stdDev Marginal variance.
#' @return Matern covariance matrix.
#' @examples
#' covMatern(dMat = NULL, range = NULL, stdDev = NULL)
#' @export
#' @import INLA
covMatern = function(dMat, range, stdDev){
  Sig = inla.matern.cov(nu = 1,
                        kappa = sqrt(8*1)/range,
                        x = dMat,
                        corr = TRUE)
  Sig = stdDev^2*Sig
}
