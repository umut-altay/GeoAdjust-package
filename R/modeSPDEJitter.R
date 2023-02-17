# this script fits SPDE model the data and generates predictions

# generate default priors for SPDE model
# from Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
# sigma0: field standard deviation
# NOTE: by default, this constructs a spde prior with unit median marginal variance
#       and median effective range equal to a fifth of the spatial range
# or use inla.spde2.pcmatern (possibly allow (1/4,4) variance here rather than (1/2,2))
# U and alpha are the threshold and probability of crossing the threshold for the variance prior
getSPDEPrior = function(mesh, U=1, alpha=0.05, medianRange=NULL) {
  size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2])))) # 1444.772
  # sigma0=1
  # range0 <- size/5
  # kappa0 <- sqrt(8)/range0
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),
  #                           B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0),
  #                           theta.prior.prec = c(0.1, 1))

  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }

  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }
  }

  if(is.null(medianRange))
    range0 <- size/5
  else
    range0 = medianRange
  spde = INLA::inla.spde2.pcmatern(mesh, prior.range=c(range0, 0.5), prior.sigma = c(U, alpha))
  spde
}

# get a reasonable default mesh triangulation for the SPDE model for the Kenya data
getSPDEMeshKenya = function(locs=NULL, n=5000, max.n=5000, doPlot=FALSE, max.edge=c(7, 200),
                            offset=-.08, cutoff=4, jitterAmount=max.edge[1]/4, seed=123) {

  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }

  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }
  }


  if(is.null(locs)) {
    # jitter the locations used to create the mesh so that they do not always lie on mesh points
    locs=cbind(jitter(dat$east, amount=jitterAmount), jitter(dat$north, amount=jitterAmount))
  }

  # generate mesh on R2
  mesh = INLA::inla.mesh.2d(loc.domain = locs, n=n, max.n=max.n, offset=offset, cutoff=cutoff, max.edge=max.edge)

  # plot the mesh if user wants
  # if(doPlot) {
  #   plot(mesh)
  #   plotMapDat(project=TRUE, border="blue")
  # }

  mesh
}
