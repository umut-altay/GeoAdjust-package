#' Prepares input data list for the model estimation with estimateModel()
#' function.
#'
#' @param response A list containing the number of trials (ns) and number of
#' successes (ys) for the binomial model, response
#' values (ys) for the Gaussian model or the Poisson counts for the Poisson
#' model.
#' @param locObs An sf class multipoint object containing the coordinates of DHS
#' survey cluster centers. (It should contain the crs information.)
#' @param likelihood A value indicating which likelihood model should be used
#' (0, 1 or 2 for Gaussian, binomial or Poisson, respectively).
#' @param jScale Jittering scale, where 1 represents the default DHS jittering
#' scheme. It allows experimenting with larger jittering scales. Otherwise should
#' remain as 1.
#' @param urban A vector containing the urbanization classification of the
#' administrative area that each cluster center is initially located within
#' (U for urban and R for rural).
#' @param mesh.s A triangulation mesh. It should be constructed in the
#' target_crs.
#' @param adminMap An sf class multipolygon object. It contains the borders of
#' the administrative area level that was respected while the cluster centers
#' were initially being jittered (can be obtained from https://gadm.org). It
#' should contain crs information.
#' @param covariateData A list containing the covariates as SpatRaster objects.
#' Each SpatRaster should contain the crs information.
#' @param target_crs A projection string representing the desired coordinate
#' reference system according to which, the integration rings and integration
#' points will be located. The measurement unit of the target_crs should be in
#' kilometers.
#' @return A list containing the data input for estimateModel() function.
#' @examples
#' path1 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' path2 <- system.file("extdata", "exampleMesh.rda", package = "GeoAdjust")
#' load(path1)
#' load(path2)
#' crs_KM = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
#' crs_Degrees = "+proj=longlat +datum=WGS84"
#' locObs = data.frame(long = surveyData$long, lat = surveyData$lat)
#' locObs = sf::st_as_sf(locObs, coords=c("long","lat"), crs = crs_Degrees)
#' inputData <- prepareInput(response = list(ys = surveyData$ys, ns = surveyData$ns),
#' locObs = locObs, likelihood = 1, jScale = 1,
#' urban = surveyData$urbanRural, mesh.s = exampleMesh, adminMap = adm1,
#' covariateData = NULL, target_crs = crs_KM)
#' @export
prepareInput = function(response=NULL, locObs=NULL, likelihood, jScale=1,
                         urban=NULL, mesh.s=NULL,
                         adminMap=NULL,
                         covariateData=NULL, target_crs){

  locObs = sf::st_transform(locObs, target_crs)
  adminMap = sf::st_transform(adminMap, target_crs)

  # extract arguments
  flag2 = likelihood

  # number of observed locations
  nLoc = length(sf::st_coordinates(locObs)[,1])

  #response variable Gaussian/Binomial/Poisson
  if (flag2 == 0){
    ys = response[["ys"]]
    ns = rep(1, nLoc)
  } else if (flag2 == 1){
    ys = response[["ys"]]
    ns = response[["ns"]]
  } else{
    ys = response[["ys"]]
    ns = rep(1, nLoc)
  }
  #
  #
  # spde components

  USpatial = 1
  alphaSpatial = 0.05
  #
  #jittering the points a bit just to make a mesh
  spdeComponents = fmesher::fm_fem(mesh.s)
  A.mesher = fmesher::fm_evaluator(mesh = mesh.s, loc = cbind(sf::st_coordinates(locObs)[,1], sf::st_coordinates(locObs)[,2]))
  A.proj = A.mesher[["proj"]][["A"]]

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

  #adminMap$OBJECTID = 1:length(adminMap$NAME_1)
  intPointInfo = makeAllIntegrationPoints(coords = locObs, urbanVals = urbanVals,
                                          numPointsUrban=1+15*4, numPointsRural=1+15*9,
                                          scalingFactor = jScale,
                                          JInnerUrban=5, JOuterUrban=0,
                                          JInnerRural=5, JOuterRural=5,
                                          adminMap=adminMap,
                                          nSubAPerPoint=10,
                                          nSubRPerPoint=10,
                                          testMode=FALSE)

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

  UrbanCoor = data.frame(x = coordsUrban[,1], y = coordsUrban[,2])
  UrbanCoor = sf::st_as_sf(UrbanCoor, coords=c("x","y"), crs = sf::st_crs(locObs))

  RuralCoor = data.frame(x = coordsRural[,1], y = coordsRural[,2])
  RuralCoor = sf::st_as_sf(RuralCoor, coords=c("x","y"), crs = sf::st_crs(locObs))

if (!is.null(covariateData)){
  for (i in 1:length(covariateData)){

    crs_CovRaster = sf::st_crs(covariateData[[i]])
    coorVectorUrban = terra::vect(sf::st_transform(UrbanCoor, crs_CovRaster[["wkt"]]))
    coorVectorRural = terra::vect(sf::st_transform(RuralCoor, crs_CovRaster[["wkt"]]))

    #Extract covariate values from data rasters at dhsLocs
    assign(paste0("covariate_Urban", i), terra::extract(covariateData[[i]], coorVectorUrban)[,2])
    assign(paste0("covariate_Rural", i), terra::extract(covariateData[[i]], coorVectorRural)[,2])

  }
  #}
}
  nLoc_urban = length(sf::st_coordinates(UrbanCoor[,1]))
  nLoc_rural = length(sf::st_coordinates(RuralCoor[,1]))

  desMatrixJittUrban = as.matrix(rep(1, nLoc_urban))
  desMatrixJittRural = as.matrix(rep(1, nLoc_rural))

  if(!is.null(covariateData)){

  for(i in 1:length(covariateData)){
    desMatrixJittUrban = cbind(desMatrixJittUrban, get(paste0("covariate_Urban", i)))
    desMatrixJittRural = cbind(desMatrixJittRural, get(paste0("covariate_Rural", i)))
  }
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
               M0 = spdeComponents[["c0"]],
               M1 = spdeComponents[["g1"]],
               M2 = spdeComponents[["g2"]],
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
