#' Predicts new outcomes for the response variable
#'
#' @param obj The output object created by the "estimateModel" function that contains the TMB components for predictions
#' @param predCoords A matrix containing the coordinates of the prediction locations in kilometers
#' @return A matrix called "PredictedResponses", containing the mean, median,
#' standard deviation and the lower and the upper bounds of 95% credible intervals
#' @examples
#' predRes(obj, predCoords)
#' @export
#' @import TMB
#' @import SUMMER
predRes = function(obj, predCoords){
  obj <- normalize(obj, flag="flag1", value = 0)
  opt0 = optim(par=obj$par, fn = obj$fn, gr = obj$gr,
               method = c("BFGS"), hessian = FALSE, control=list(parscale=c(.1, .1)))
  par <- obj$env$last.par
  Qtest = obj$env$spHess(par, random = TRUE)

  #Sampling
  mu <- c(par[-c(7,8)])

  # Simulate draws
  rmvnorm_prec <- function(mu, chol_prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- chol_prec #Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z


    prec = Qtest
    L = Cholesky(prec, super = T)

    A.pred = inla.spde.make.A(mesh = mesh.s, loc = predCoords)
    #
    t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = 10000)
    parnames <- c(names(mu))

    epsilon_draws  <- t.draws[parnames == 'Epsilon_s',]
    beta_draws<- t.draws[parnames == 'beta',]

    eta.samples = covariatesAtPred%*%beta_draws + as.matrix(A.pred%*%epsilon_draws)
    #
    trueValues = covariatesAtPred%*%betaSim + u.sim
    trueValues = trueValues[,1]
    #
    nPred = length(eta.samples[,1])
    coverage = 0
    for (m in 1:nPred){
      qnt = quantile(eta.samples[m,], c(0.025, 0.975))
      if ((qnt[[1]]<trueValues[[m]])&(qnt[[2]]>trueValues[[m]])){
        coverage = coverage + 1
      }
    }
    coverage = coverage/nPred
    #
    Logscores = logs_sample(y = trueValues, dat = eta.samples)
    Logscores = mean(Logscores)
    #
    CRPS = crps_sample(y = trueValues, dat = eta.samples, method = "edf")
    CRPS = mean(CRPS)
    #
    BIAS = mean(rowMeans(eta.samples)-trueValues)
    #
    RMSE = sqrt(mean((rowMeans(eta.samples)-trueValues)^2))
    #
    if (flag2 != 0){
      #   # Convert to probability scale
      eta.samples = expit(eta.samples)
    }
    # # Find the median and sd across draws, as well as 95% intervals
    PredictedResponses <- cbind(mean = (apply(eta.samples, 1, mean)),
                                median = (apply(eta.samples, 1, median)),
                                sd     = (apply(eta.samples, 1, sd)),
                                lower = (apply(eta.samples, 1, quantile, .025, na.rm = TRUE)),  #na.rm = TRUE is newly added
                                upper = (apply(eta.samples, 1, quantile, .975, na.rm = TRUE)))

    return(PredictedResponses)

  }





}
