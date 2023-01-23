#' Prepares input data list for the model estimation with "estimateMoodel
#' function
#'
#' @param response A list containing the number of trials (ns) and number of
#' successes (ys) for the binomial response, or a list containing the response
#' values (ys) for the Gaussian response.
#' @param locObs A matrix containing the coordinates of the already jittered
#' survey cluster centers in kilometers
#' @param likelihood A value indication which likelihood will be used. (1 for
#' binomial and 0 for Gaussian)
#' @param jScale Jittering scale, where 1 represents the default DHS jittering
#' scheme
#' @param urban A vector containing the urbanization classification of the
#' administrative area that each cluster center is initially located within.
#' (U for urban and R for rural)
#' @param mesh.s A mesh created based on the country borders
#' @param adminMap A shape file containing the borders of the administrative
#' area level, which was respected while the cluster centers were initially
#' being jittered. (can be obtained from https://gadm.org)
#' @param nSubAPerPoint A value representing the number of unique
#' sub-integration point angles per integration point
#' @param nSubRPerPoint A value representing the number of unique
#' sub-integration point radii per integration point
#' @param covariateData A list containing the covariate rasters
#' @return A list containing a list of data inputs for TMB, the corresponding
#' mesh and the corresponding matrix of observation locations
#' @examples
#' \dontrun{
#' inputData <- prepareInput(response = response, locObs = locObs,
#' likelihood = likelihood, jScale = jScale,
#' urban = urban, mesh.s = mesh.s, adminMap = adminMap,
#' nSubAPerPoint = nSubAPerPoint, nSubRPerPoint = nSubRPerPoint,
#' covariateData = covariateData)
#' }
#' @export
#' @import INLA
prepareInput = function(response=NULL, locObs=NULL, likelihood, jScale=NULL,
                         urban=NULL, mesh.s=NULL,
                         adminMap=NULL, nSubAPerPoint=10, nSubRPerPoint=10,
                         covariateData=NULL){

  # extract arguments
  flag2 = likelihood

  # number of observed locations
  nLoc = length(locObs[,1])

  #response variable Gaussian/Binomial
  if (flag2 == 0){
    ys = response[["ys"]]
    ns = rep(1, nLoc)
  } else {
    ys = response[["ys"]]
    ns = response[["ns"]]
  }
  #
  #
  # spde components

  USpatial = 1
  alphaSpatial = 0.05
  #
  #jittering the points a bit just to make a mesh
  spde = getSPDEPrior(mesh.s, U=USpatial, alpha=alphaSpatial)
  A.proj = inla.spde.make.A(mesh = mesh.s, loc = cbind(locObs[,1], locObs[,2]))
  #
  # TMB input for the model that accounts for jittering

  # convert urban U/R into TRUE/FALSE
  for (i in 1:length(urban)){
    if (urban[i]=='U'){
      urban[i]='TRUE'
    }else{
      urban[i]='FALSE'
    }
  }
  urbanVals=as.logical(urban)

  intPointInfo = makeAllIntegrationPoints(coords = cbind(locObs[,1], locObs[,2]), urbanVals,
                                          numPointsUrban=1+15*4, numPointsRural=1+15*9,
                                          scalingFactor = jScale,
                                          JInnerUrban=5, JOuterUrban=0,
                                          JInnerRural=5, JOuterRural=5,
                                          adminMap=adminMap,
                                          nSubAPerPoint=nSubAPerPoint,
                                          nSubRPerPoint=nSubRPerPoint,
                                          testMode=FALSE)


  # if(testMode) {
  #   return(intPointInfo)
  # }


  xUrban = intPointInfo$xUrban
  yUrban = intPointInfo$yUrban
  wUrban = intPointInfo$wUrban
  xRural = intPointInfo$xRural
  yRural = intPointInfo$yRural
  wRural = intPointInfo$wRural

  # get the long set of coordinates
  coordsUrban = cbind(c(xUrban), c(yUrban))
  coordsRural = cbind(c(xRural), c(yRural))

  # separate observations by urbanicity
  ysUrban = ys[urbanVals]
  ysRural = ys[!urbanVals]
  nsUrban = ns[urbanVals]
  nsRural = ns[!urbanVals]

  # Covariate values at jittered locations and integration points:

  # Convert them into degrees and Mollweide
  coordsUrbanDegree = convertKMToDeg(coordsUrban)
  coordsRuralDegree = convertKMToDeg(coordsRural)

  # coordsUrbanDegree = cbind(coordsUrbanDegree[,1], coordsUrbanDegree[,2])
  # coordsRuralDegree = cbind(coordsRuralDegree[,1], coordsRuralDegree[,2])
  # Convert them into SpatialPoints
  coordsUrbanDegree = SpatialPoints(cbind(coordsUrbanDegree[,1], coordsUrbanDegree[,2]), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), bbox = NULL)
  coordsRuralDegree = SpatialPoints(cbind(coordsRuralDegree[,1], coordsRuralDegree[,2]), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), bbox = NULL)

  # Extract the corresponding covariate values seperately for urban/rural
  for (i in 1:length(covariateData)){
    #Extract covariate values from data rasters at dhsLocs
    assign(paste0("covariate_Urban", i), raster::extract(covariateData[[i]], coordsUrbanDegree, ncol=2))
    assign(paste0("covariate_Rural", i), raster::extract(covariateData[[i]], coordsRuralDegree, ncol=2))

  }

  nLoc_urban = length(coordsUrbanDegree@coords[,1])
  nLoc_rural = length(coordsRuralDegree@coords[,1])

  desMatrixJittUrban = rep(1, nLoc_urban)
  desMatrixJittRural = rep(1, nLoc_rural)

  for(i in 1:length(covariateData)){
    desMatrixJittUrban = cbind(desMatrixJittUrban, get(paste0("covariate_Urban", i)))
    desMatrixJittRural = cbind(desMatrixJittRural, get(paste0("covariate_Rural", i)))
  }

  desMatrixJittUrban[is.nan(desMatrixJittUrban)] = 0
  desMatrixJittUrban[is.na(desMatrixJittUrban)] = 0

  desMatrixJittRural[is.nan(desMatrixJittRural)] = 0
  desMatrixJittRural[is.na(desMatrixJittRural)] = 0
  #
  desMatrixJittUrban = as.matrix(desMatrixJittUrban)
  desMatrixJittRural = as.matrix(desMatrixJittRural)


  n_integrationPointsUrban = ncol(wUrban)
  n_integrationPointsRural = ncol(wRural)

  #
  # # Construct projection matrices, and get other relevant info for TMB
  out = makeJitterDataForTMB(intPointInfo, ys , urbanVals, ns, spdeMesh = mesh.s)
  ysUrban = out$ysUrban
  ysRural = out$ysRural
  nsUrban = out$nsUrban
  nsRural = out$nsRural
  AUrban = out$AUrban
  ARural = out$ARural
  #
  # Compile inputs for TMB
  data <- list(num_iUrban = length(ysUrban),  # Total number of urban observations
               num_iRural = length(ysRural),  # Total number of rural observations
               num_s = mesh.s[['n']], # num. of vertices in SPDE mesh
               y_iUrban   = ysUrban, # num. of pos. urban obs in the cluster
               y_iRural   = ysRural, # num. of pos. rural obs in the cluster
               n_iUrban   = nsUrban,  # num. of urban exposures in the cluster
               n_iRural   = nsRural,  # num. of rural exposures in the cluster
               n_integrationPointsUrban = n_integrationPointsUrban,
               n_integrationPointsRural = n_integrationPointsRural,
               wUrban = wUrban,
               wRural = wRural,
               X_betaUrban = desMatrixJittUrban,
               X_betaRural = desMatrixJittRural,
               M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
               M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
               M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
               AprojUrban = AUrban,             # Projection matrix (urban)
               AprojRural = ARural,             # Projection matrix (rural)
               options = c(1, ## if 1, use normalization trick
                           1), ## if 1, run adreport
               # normalization flag.
               flag1 = 1,
               flag2 = flag2 #(0/1 for Gaussian/Binomial)
  )

  return(data)
}