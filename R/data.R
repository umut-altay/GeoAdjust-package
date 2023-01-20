#' Example data set
#'
#' Simulated data set to mimic the DHS cluster level survey data
#'
#' @format ## `clusterData`
#' A data frame with 1,500 rows and 8 columns:
#' \describe{
#'   \item{clusterID}{Simulated cluster ID}
#'   \item{long}{Simulated longitudes for each cluster center}
#'   \item{lat}{Simulated latitudes for each cluster center}
#'   \item{ys}{Simulated binomial successes for each cluster center}
#'   \item{ns}{Simulated binomial trials for each cluster center}
#'   \item{east}{Simulated UTM zone37 easting coordinates for each cluster center}
#'   \item{north}{Simulated UTM zone37 northing coordinates for each cluster center}
#'   \item{urbanRural}{Simulated urbanicity strata for each cluster center}
#'   ...
#' }
#' @source the data is simulated
"clusterData"
