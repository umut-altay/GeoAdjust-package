#' Predicts model outcomes at new locations.
#'
#' @param mesh.s A mesh created based on the country borders.
#' @param covariateData A list containing the covariate rasters. Each covariate
#' raster should be SpatRast object and contain crs information.
#' @param obj The optimized core model object returned by estimateModel() function.
#' @param draws A matrix containing 10.000 sampled values for each covariate effect
#' size and 10.000 sampled values of random effect coefficients for each mesh node.
#' It is one of the elements of the returning output list of estimateModel() function.
#' @param predCoords An sf class POINT object containing the coordinates of the
#' prediction locations. It should contain crs information.
#' @param flag A value indicating the type of the likelihood that will be used.
#' Pass 0 for Gaussian, 1 for binomial and 2 for Poisson likelihoods.
#' @return A matrix containing the mean, median,standard deviation and the lower
#' and the upper bounds of 95% credible intervals of the predictions.
#' @examples
#' \donttest{
#' path1 <- system.file("extdata", "exampleInputData.rda", package = "GeoAdjust")
#' path2 <- system.file("extdata", "exampleMesh.rda", package = "GeoAdjust")
#' path3 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' load(path2)
#' load(path3)
#' results <- estimateModel(data = exampleInputData,
#' options = list(random = 1, covariates = 1), priors = list(beta = c(0,1),
#' range = 114), USpatial = 1, alphaSpatial = 0.05, UNugget = 1, alphaNug = 0.05, n.sims = 1000)
#' crs_KM = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
#' exampleGrid <- gridCountry(admin0 = adm0, res = 5, target_crs = crs_KM)
#' pred = predRes(obj = results[["obj"]],
#' predCoords = exampleGrid[["loc.pred"]],
#' draws = results[["draws"]],
#' mesh.s = exampleMesh, covariateData = NULL, flag = 1)
#' }
#' @importFrom stats median quantile sd
#' @export
predRes = function(obj = NULL, predCoords  = NULL, draws  = NULL, covariateData  = NULL, mesh.s  = NULL, flag  = NULL){

  t.draws = draws
  prd_coor = sf::st_coordinates(predCoords)
  prd_coor = as.matrix(cbind(prd_coor[,1], prd_coor[,2]))
  A.pred = fmesher::fm_evaluator(mesh = mesh.s, loc = prd_coor)
  A.pred = A.pred[["proj"]][["A"]]
  par <- obj$env$last.par

  idx1 = which(names(par)=="theta")

  mu = par[-c(idx1)]

  parnames <- c(names(mu))
  epsilon_draws  <- t.draws[parnames == 'Epsilon_s',]

  beta_draws<- t.draws[parnames == 'beta',]

    if(!is.null(covariateData)){

      for (i in 1:length(covariateData)){

        crs_CovRaster = sf::st_crs(covariateData[[i]])
        coorVectorPred = terra::vect(sf::st_transform(predCoords, crs_CovRaster[["wkt"]]))

        #Extract covariate values from data rasters at predCoords
        assign(paste0("covariatePred", i), terra::extract(covariateData[[i]], coorVectorPred)[,2])

      }

}
    nLoc_pred = length(sf::st_coordinates(predCoords)[,1])

    covariatesAtPred = rep(1, nLoc_pred)

    if(!is.null(covariateData)){
    for(i in 1:length(covariateData)){
      covariatesAtPred = cbind(covariatesAtPred, get(paste0("covariatePred", i)))
    }
}
    covariatesAtPred[is.nan(covariatesAtPred)] = 0
    covariatesAtPred[is.na(covariatesAtPred)] = 0

    if(!is.null(covariateData)){
    eta.samples = covariatesAtPred%*%beta_draws + as.matrix(A.pred%*%epsilon_draws)
    } else {
      eta.samples = as.matrix(covariatesAtPred)%*%t(as.matrix(beta_draws)) + as.matrix(A.pred%*%epsilon_draws)
    }
    #

    if (flag == 1){
      #   # Convert to probability scale
      eta.samples = SUMMER::expit(eta.samples)
    } else if (flag ==2){
      eta.samples = exp(eta.samples)
    } else{
      eta.samples = eta.samples
    }

    # # Find the median and sd across draws, as well as 95% intervals
    PredictedResponses <- cbind(mean = (apply(eta.samples, 1, mean)),
                                median = (apply(eta.samples, 1, median)),
                                sd     = (apply(eta.samples, 1, sd)),
                                lower = (apply(eta.samples, 1, quantile, .025, na.rm = TRUE)),  #na.rm = TRUE is newly added
                                upper = (apply(eta.samples, 1, quantile, .975, na.rm = TRUE)))



    return(PredictedResponses = PredictedResponses)

  }



