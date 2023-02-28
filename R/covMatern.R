#' Creates a Matern covariance matrix.
#'
#' @param dMat A distance matrix between the locations.
#' @param range Spatial range in kilometers.
#' @param stdDev The marginal variance.
#' @return Matern covariance matrix.
#' @examples
#' if(requireNamespace("INLA")){
#' path1 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' loc <- cbind(surveyData$east, surveyData$north)
#' space.range <- 114
#' space.sigma <- 1
#' covMat <- covMatern(dMat = as.matrix(dist(loc)),
#' range = space.range, stdDev = space.sigma)
#' }
#' @export
covMatern = function(dMat = NULL, range = NULL, stdDev = NULL){

  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }

  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }
  }

  Sig = INLA::inla.matern.cov(nu = 1,
                        kappa = sqrt(8*1)/range,
                        x = dMat,
                        corr = TRUE)
  Sig = stdDev^2*Sig
}
