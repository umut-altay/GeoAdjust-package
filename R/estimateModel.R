#' Estimates model parameters
#'
#' @param data A data input list that is created by prepareInput() function.
#' @param nNodes number of mesh nodes.
#' @param n.sims number of samples to be drawn for each model parameter
#' @param options A list containing two components, namely, random and covariates, representing the spatial random field and covariates.
#' Values of 1 and 0 turn the accounting for jittering in these components on and off.
#' @param priors A list of six components. Beta is a vector of two elements and passes the parameters of the Gaussian prior that will be assigned
#' to the covariates (including the intercept). The first element of it is the mean and the second one is the
#' standard deviation of Gaussian prior. Range is a value representing the median range in kilometers, which will be used for constructing the PC (Penalized-complexity) priors.
#' USpatial and alphaSpatial are the threshold and probability of crossing the threshold for the variance prior.
#' UNugget and alphaNug are the threshold and probability of crossing the threshold for the prior on the nugget standard deviation.
#' UNugget and alphaNug should be included in the priors argument, but they will only be used when the likelihood is Gaussian.
#' @return Model estimation results of class called res. The output consists of four elements:
#' A data frame containing the estimated model parameters and the corresponding 95% credible interval lengths,
#' The optimized core model object from autodifferentiation of TMB,
#' A matrix containing the sampled coefficient effect sizes and the random effect coefficients,
#' A character string indicating the likelihood type in the model.
#' @examples
#' \donttest{
#' path1 <- system.file("extdata", "exampleInputData.rda", package = "GeoAdjust")
#' path2 <- system.file("extdata", "exampleMesh.rda", package = "GeoAdjust")
#' load(path1)
#' load(path2)
#' nNodes = exampleMesh[['n']]
#' results <- estimateModel(data = exampleInputData, nNodes = nNodes,
#' options = list(random = 1, covariates = 1), priors = list(beta = c(0,1),
#' range = 114, USpatial = 1, alphaSpatial = 0.05, UNugget = 1, alphaNug = 0.05), n.sims = 1000)
#' }
#' @importFrom stats optim
#' @export
estimateModel = function(data = NULL, nNodes = NULL, options = NULL, priors = NULL, n.sims = NULL){

  flagRandomField = options[["random"]]
  flagCovariates = options[["covariates"]]

  beta_pri = priors[["beta"]]
  rangeMaternPri = priors[["range"]]
  USpatial  = priors[["USpatial"]]
  alphaSpatial  = priors[["alphaSpatial"]]
  nugStd = priors[["UNugget"]]
  alphaNug = priors[["alphaNug"]]

  matern_pri = c(rangeMaternPri, 0.5, USpatial , alphaSpatial)
  nug_pri = c(nugStd, alphaNug)


  data[["flagRandomField"]] = flagRandomField
  data[["flagCovariates"]] = flagCovariates
  data[["beta_pri"]] = beta_pri
  data[["matern_pri"]] = matern_pri
  data[["nug_pri"]] = nug_pri


  if (data[["flag2"]] ==0){
    likelihood = "Gaussian"
  } else if (data[["flag2"]] ==1){
    likelihood =  "binomial"
  } else {
    likelihood =  "Poisson"
    }


  nCov = length(data[["X_betaUrban"]][1,]) # number of covariates in the model, including the intercept

  log_tau = 1.87 # Log tau (i.e. log spatial precision, Epsilon)
  log_kappa = -3.69 # SPDE parameter related to the range

  if (data[["flag2"]] != 0){
    tmb_params <- list(beta=rep(0, nCov),
                       Epsilon_s = rep(0, nNodes),#,  RE on mesh vertices
                       theta = c(log_kappa, log_tau)
                       )
  }else{
    log_nug_std_initial = -2
    tmb_params <- list(beta=rep(0, nCov),
                       Epsilon_s = rep(0, nNodes),#,  RE on mesh vertices
                       theta = c(log_kappa, log_tau, log_nug_std_initial))
  }

  # random effects
  rand_effs <- c("Epsilon_s", "beta")



  obj <- TMB::MakeADFun(data=data,
                   parameters=tmb_params,
                   random = rand_effs,
                   hessian=TRUE,
                   DLL='GeoAdjust')

  flag1 = 1
  obj <- TMB::normalize(obj, flag="flag1", value = 0)
  if(data[["flag2"]] == 0){
    opt0 = stats::optim(par=obj$par, fn = obj$fn, gr = obj$gr,
                        method = c("BFGS"), hessian = FALSE, control=list(parscale=c(.1, .1, .1)))
  }else{
    opt0 = stats::optim(par=obj$par, fn = obj$fn, gr = obj$gr,
                        method = c("BFGS"), hessian = FALSE, control=list(parscale=c(.1, .1)))
  }

  par <- obj$env$last.par
  Qtest = obj$env$spHess(par, random = TRUE)

  idx1 = which(names(par)=="theta")

  mu = par[-c(idx1)]

  log_kappa = par[[idx1[1]]]
  log_tau   = par[[idx1[2]]]


  if(likelihood =="Gaussian"){
    logNugStd_estimate = par[[idx1[3]]]
  }

  idxCov = which(names(par)=="beta")
  beta_estimate = c(par[idxCov])
  range_estimate = sqrt(8.0)/exp(log_kappa)
  sigma_estimate = 1.0 / sqrt(4.0*3.14159265359*exp(2.0*log_tau)*exp(2.0*log_kappa))

  # Simulate draws
  rmvnorm_prec <- function(mu, chol_prec, n.sims) {
    z <- matrix(stats::rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- chol_prec #Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }
  prec = Qtest
  L = Matrix::Cholesky(prec, super = T)
  t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = n.sims)
  parnames <- c(names(mu))

  beta_draws<- t.draws[parnames == 'beta',]

  if(is.null(dim(beta_draws))){
    intercept = beta_draws
  } else {
  intercept = beta_draws[1,]
}
  if(nCov>1){
  npar = nCov-1
  parNames = rep(0, npar)
  for(i in 1:npar){
  parNames[i] = paste0("beta", i)
  assign(paste0("beta", i), beta_draws[i+1,])
  }
  } else{
    npar=0
    parNames = NULL
}


  CIlength = rep(0, npar+1)
  CIlength[1] = stats::quantile(intercept, probs = 0.975)- c(stats::quantile(intercept, probs = 0.025))
  if(nCov>1){
  for(i in (1:npar)){
    b=paste0("beta", i)
  CIlength[i+1] = (stats::quantile(eval(parse(text = b)), probs = 0.975)- c(stats::quantile(eval(parse(text = b)), probs = 0.025)))
  }
  }


  if(is.null(parNames)){
    if(likelihood =="Gaussian"){
      res = data.frame(parameters = c("range", "sigma", "nuggetSd", "intercept"),
                       estimates = c(range_estimate, sigma_estimate, exp(logNugStd_estimate), beta_estimate),
                       CI_Length = c("", "", "", CIlength))

    }else{
    res = data.frame(parameters = c("range", "sigma", "intercept"),
                     estimates = c(range_estimate, sigma_estimate, beta_estimate),
                     CI_Length = c("", "", CIlength))
    }



  } else {
    if(likelihood =="Gaussian"){
      res = data.frame(parameters = c("range", "sigma", "nuggetSd", "intercept", parNames),
                       estimates = c(range_estimate, sigma_estimate, exp(logNugStd_estimate), beta_estimate),
                       CI_Length = c("", "", "", CIlength))

    }else{

  res = data.frame(parameters = c("range", "sigma", "intercept", parNames),
                    estimates = c(range_estimate, sigma_estimate, beta_estimate),
                   CI_Length = c("", "", CIlength))
    }
}

  est = list(res=res, obj=obj,draws =t.draws, likelihood = likelihood)
  class(est) = "res"

  return(est)

}
