#' Predicts model outcomes at new locations.
#'
#' @param nCov A value showing the number of covariates (including the intercept).
#' @param mesh.s A mesh created based on the country borders.
#' @param covariateData A list containing the covariate rasters.
#' @param obj The optimized core model object returned by estimateModel() function.
#' @param draws A matrix containing 10.000 sampled values for each covariate effect size and 10.000 sampled values of random effect coefficients for each mesh node.
#' It is one of the elements of the returning output list of estimateModel() function.
#' @param predCoords A matrix containing the coordinates of the prediction locations in kilometers (UTM zone:37).
#' @param flag A value indicating the type of the likelihood that will be used. Pass 0 for Gaussian, 1 for binomial and 2 for Poisson likelihoods.
#' @return A matrix containing the mean, median,standard deviation and the lower and the upper bounds of 95% credible intervals of the predictions.
#' @examples
#' \donttest{
#' if(requireNamespace("INLA")){
#' path1 <- system.file("extdata", "exampleInputData.rda", package = "GeoAdjust")
#' path2 <- system.file("extdata", "exampleMesh.rda", package = "GeoAdjust")
#' path3 <- system.file("extdata", "exampleGrid.rda", package = "GeoAdjust")
#' load(path1)
#' load(path2)
#' load(path3)
#' nNodes = exampleMesh[['n']]
#' results <- estimateModel(data = exampleInputData, nNodes = nNodes,
#' options = list(random = 1, covariates = 1), priors = list(beta = c(0,1),
#' range = 114, USpatial = 1, alphaSpatial = 0.05, UNugget = 1, alphaNug = 0.05))
#' pred = predRes(obj = results[["obj"]],
#' predCoords = cbind(exampleGrid[["loc.pred"]]["east"],
#' exampleGrid[["loc.pred"]]["north"]),
#' draws = results[["draws"]], nCov = 1,
#' mesh.s = exampleMesh, covariateData = NULL, flag = 1)
#' }
#' }
#' @importFrom stats median quantile sd
#' @export
predRes = function(obj = NULL, predCoords  = NULL, draws  = NULL, nCov  = NULL, covariateData  = NULL, mesh.s  = NULL, flag  = NULL){

  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }

  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }
  }


  t.draws = draws
  predCoords = as.matrix(cbind(predCoords["east"], predCoords["north"]))
  A.pred = INLA::inla.spde.make.A(mesh = mesh.s, loc = predCoords)

  par <- obj$env$last.par

  idx1 = which(names(par)=="theta")
  # if(flag ==0){
  #   idx2 = which(names(par)=="log_nug_std")
  #   idx3 = which(names(par)=="hp")
  #   mu = par[-c(idx1, idx2, idx3)]
  # }else{
    mu = par[-c(idx1)]
  # }

  parnames <- c(names(mu))
  epsilon_draws  <- t.draws[parnames == 'Epsilon_s',]


  beta_draws<- t.draws[parnames == 'beta',]

    predCoordsDegree = convertKMToDeg(predCoords)
    predCoordsDegree = sp::SpatialPoints(cbind(predCoordsDegree[,1], predCoordsDegree[,2]), proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"), bbox = NULL)

    if(!is.null(covariateData)){
    # Extract the corresponding covariate values at prediction locations
    for (i in 1:length(covariateData)){
      #Extract covariate values from data rasters at predCoords
      assign(paste0("covariatePred", i), raster::extract(covariateData[[1]], predCoordsDegree, ncol=2))
    }
}
    nLoc_pred = length(predCoords[,1])

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



