#' Predicts new outcomes for the response variable
#'
#' @param nCov Number of covariates (including the intercept)
#' @param mesh.s A mesh created based on the country borders
#' @param covariateData A list containing the covariate rasters
#' @param obj The optimized core model object returned by estimateModel() function.
#' @param draws A matrix containing 10.000 sampled values for each covariate effect size and 10.000 sampled values of random effect coefficients for each mesh node.
#' @param predCoords A matrix containing the coordinates of the prediction locations in UTM zone:37
#' @param flag A value indicating the type of the likelihood that will be used. Pass 0 for Gaussian, 1 for binomial and 2 for Poisson likelihoods.
#' @return A list containing a matrix called "PredictedResponses", and another matrix called "eta.samples".  The matrix "PredictedResponses" contains the mean, median,
#' standard deviation and the lower and the upper bounds of 95% credible intervals. The matrix "eta.samples" contains the sampled values
#' @examples
#' \dontrun{
#' predRes(obj = obj, predCoords = predCoords, nCov = nCov, mesh = mesh, covariateData = covariateData)
#' }
#' @export
predRes = function(obj = obj, predCoords = predCoords, draws = draws, nCov = nCov, covariateData = covariateData, mesh.s = mesh.s, flag = flag){

  t.draws = draws
  predCoords = as.matrix(cbind(predCoords["east"], predCoords["north"]))
  A.pred = inla.spde.make.A(mesh = mesh.s, loc = predCoords)

  par <- obj$env$last.par
  mu = par[names(par) != c("log_tau","log_kappa")]
  parnames <- c(names(mu))
  epsilon_draws  <- t.draws[parnames == 'Epsilon_s',]
  beta_draws<- t.draws[parnames == 'beta',]
  # median = sd = quantile = NULL
  # par <- obj$env$last.par
  # Qtest = obj$env$spHess(par, random = TRUE)
  #
  # #Sampling
  # mu <- c(par[-c(nCov+1,nCov+2)]) # nCov+1 refers to log_tau
  #                                 # nCov+2 refers to log_kappa
  # # Simulate draws
  # rmvnorm_prec <- function(mu, chol_prec, n.sims) {
  #   z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  #   L <- chol_prec #Cholesky(prec, super=TRUE)
  #   z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  #   z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  #   z <- as.matrix(z)
  #   mu + z
  #
  # }
    # prec = Qtest
    # L = Matrix::Cholesky(prec, super = T)

    # A.pred = INLA::inla.spde.make.A(mesh = mesh.s, loc = as.matrix(predCoords))
    # #
    # t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = 10000)
    # parnames <- c(names(mu))
    #
    # epsilon_draws  <- t.draws[parnames == 'Epsilon_s',]
    # beta_draws<- t.draws[parnames == 'beta',]

    predCoordsDegree = convertKMToDeg(predCoords)
    predCoordsDegree = sp::SpatialPoints(cbind(predCoordsDegree[,1], predCoordsDegree[,2]), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), bbox = NULL)

    # Extract the corresponding covariate values at prediction locations
    for (i in 1:length(covariateData)){
      #Extract covariate values from data rasters at predCoords
      assign(paste0("covariatePred", i), raster::extract(covariateData[[1]], predCoordsDegree, ncol=2))
    }

    nLoc_pred = length(predCoords[,1])

    covariatesAtPred = rep(1, nLoc_pred)

    for(i in 1:length(covariateData)){
      covariatesAtPred = cbind(covariatesAtPred, get(paste0("covariatePred", i)))
    }

    covariatesAtPred[is.nan(covariatesAtPred)] = 0
    covariatesAtPred[is.na(covariatesAtPred)] = 0

    eta.samples = covariatesAtPred%*%beta_draws + as.matrix(A.pred%*%epsilon_draws)
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



