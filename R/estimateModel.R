#' Estimate model parameters
#'
#' @param Data A data input list that is created by prepare_input() function
#' @param options A list containing two components, namely "random" and " covariates", representing the spatial random field and covariates.
#' Values of 1 and 0 turn the accounting for jittering in these components on and off.
#' Example usage : options = list(random = 0, covariates = 1) --> then jittering is accounted for only in the covariates.
#' @param priors A list containing two components, namely "beta" and "range". Beta is a vector of two elements and passes the parameters of the Gaussian prior that will be assigned
#' to the covariates (including the intercept). The first element of it is the mean and the second one is the
#' standard deviation of Gaussian prior. Range is a value representing the median range in kilometers, which will be used for constructing the PC (Penalized-complexity) priors.
#' Example usage : priors = list(beta = c(0,1), range = 114)
#' @return A data frame called "res", containing the estimated model parameters, and "obj" which contains the required components for predictions
#' with the model (if wanted).
#' @examples
#' estimateModel(data, options, priors)
#' @export
#' @import TMB
estimateModel = function(data = NULL, options = NULL, priors = NULL){

  flagRandomField = options[["random"]]
  flagCovariates = options[["covariates"]]

  beta_pri = priors[["beta"]]
  rangeMaternPri = priors[["range"]]
  # USpatial  = priors[["USpatial"]]
  # alphaSpatial  = priors[["alphaSpatial"]]

  matern_pri = c(rangeMaternPri, 0.5, USpatial = USpatial , alphaSpatial = alphaSpatial)

  data[["flagRandomField"]] = flagRandomField
  data[["flagCovariates"]] = flagCovariates
  data[["beta_pri"]] = beta_pri
  data[["matern_pri"]] = matern_pri

  nCov = dim(X_betaUrban)[[2]] # number of covariates in the model, including the intercept

  tmb_params <- list(beta=rep(0, nCov),
                     #log_tau = 5, # Log tau (i.e. log spatial precision, Epsilon)
                     log_tau = 2.74, # Log tau (i.e. log spatial precision, Epsilon)
                     log_kappa = -4, # SPDE parameter related to the range
                     Epsilon_s = rep(0, mesh.s[['n']])#,  RE on mesh vertices
                     #log_nug_std = log(sqrt(0.1))
  )

  # random effects
  rand_effs <- c("Epsilon_s", "beta")

  obj <- MakeADFun(data=data,
                   parameters=tmb_params,
                   random = rand_effs,
                   hessian=TRUE,
                   DLL='myPackage')

  obj <- normalize(obj, flag="flag1", value = 0)

  opt0 = optim(par=obj$par, fn = obj$fn, gr = obj$gr,
               method = c("BFGS"), hessian = FALSE, control=list(parscale=c(.1, .1)))

  par <- obj$env$last.par

  log_tau = par[[7]]
  log_kappa = par[[8]]

  beta_estimate = c(par[-c(7,8)])
  range_estimate = sqrt(8.0)/exp(log_kappa)
  sigma_estimate = 1.0 / sqrt(4.0 * 3.14159265359 *
                          exp(2.0 * log_tau) * exp(2.0 * log_kappa))

  parNames = rep(0, nCov)
  for(i in 1:nCov){
  parNames[i] = paste0("beta", i)
  }

  res = data.frame(parameters = c("range", "sigma", "intercept", parNames),
                    estimates = c(range_estimate, sigma_estimate, beta_estimate))
  return(res, obj)

}
