# constructs the points over which to integrate the likelihood
# with respect to the jittering distribution, and their weights,
# where the weights are related to the jittering distribution
# Arguments:
# urban: whether or not to generate integration points from the
#        urban or rural jittering distributions.
# numPoints: the number of integration points. 1 goes in the first
#            ring, and the rest are evenly split among the other
#            rings
# scalingFactor: adjust typical DHS jitter distance by a scaling factor
# JInner: the number of `inner' rings within 5km, including the first central ring
# JOuter: the number of `outer' rings beyond 5km. In urban case, only
#         M=JInner+JOuter is used
# integrationPointType: either the integration point is set as the
#                       center of mass within the integration area
#                       or the midpoint of the radius and angle
# verbose: whether to run the function in verbose mode
#' @importFrom stats aggregate
getIntegrationPoints = function(urban=TRUE, numPoints=ifelse(urban, 11, 16),
                                scalingFactor,
                                JInner=3, JOuter=ifelse(urban, 0, 1),
                                integrationPointType=c("mean", "midpoint"),
                                verbose=TRUE) {
  integrationPointType = match.arg(integrationPointType)

  M = JInner + JOuter
  ms = c(1, rep(floor((numPoints - 1) / (M - 1)), M - 1))
  msInner = ms[1:JInner]
  msOuter = ms[(JInner + 1):M]
  if(sum(ms) != numPoints) {
    stop("(numPoints - 1)/M must be a whole number")
  }

  if(urban) {
    # 2 is largest possible displacement distance
    rsInner = 2 * cumsum(ms) / sum(ms)
    rsOuter = NULL
    # r2 = (m2 + 1) * r3 / sum(ms)
    # r1 = r3 / (1 + m2 + m3)
  } else {
    rsInner = 5 * cumsum(msInner) / sum(msInner)
    rsOuter = 5 + 5 * cumsum(msOuter) / sum(msOuter)
    # r3 = 10 # largest possible displacement distance for 1% of data
    # r2 = 5 # largest possible displacement distance for 99% of data
    # r1 = r2 / (1 + m2)
  }
  # r1 = r2 / sqrt(m2 + 1) # for uniform distribution
  rs = c(rsInner, rsOuter)

  # areas for each ring segment
  As = pi * rs^2 # get areas of M discs of radius r
  As = As - c(0, As[1:(M-1)]) # subtract area of the smaller desk to get area of annulus
  As = As / ms # divide annulus area by number of integration areas per annulus
  # A1 = pi * r1^2
  # A2 = (pi * r2^2 - A1) / m2
  # A3 = (pi * r3^2 - pi * r2^2) / m3

  # these probabilities and weights were for a uniform distribution:
  # p1 = p2 = 99/100 * 1/(pi*r2^2) + 1/100 * 1/(pi*r3^2)
  # p3 = 1/100 * 1/(pi*r3^2)
  #
  # w1 = p1 * A1
  # w2 = p2 * A2
  # w3 = p3 * A3

  # helper function for calculating values of pi(s_i | s_ijk^*) for
  # integration points in ring 1, 2, or 3
  densityFunTemp = function(d, Di=rs[M]) {
    out = 1 / (Di * 2 * pi * d)
    out[d >= Di] = 0
    out[d < 0] = NaN
    out
  }
  densityFunFinal = function(d) {
    if(urban) {
      densityFunTemp(d, rs[M])
    } else {
      densityFunTemp(d, rs[JInner]) * 99 / 100 + densityFunTemp(d, rs[M]) * 1 / 100
    }
  }

  # p1 = Inf
  # p2 = densityFunFinal(mean(c(r1, r2)))
  # p3 = densityFunFinal(mean(c(r2, r3)))
  rsIntegrationPointsMidpoint = c(0, rs[-M] + diff(rs) / 2) # midpoint solution
  if(integrationPointType == "mean") {
    aDiff = c(0, 2 * pi / ms[-1])
    shrinkFactor = sqrt(2 * (1 - cos(aDiff))) / aDiff # a scaling factor smaller than one
    shrinkFactor[1] = 1
    rsIntegrationPoints = rsIntegrationPointsMidpoint * shrinkFactor
    tooShrunk = rsIntegrationPoints < c(0, rs[-M])
    if(any(tooShrunk)) {
      warning(paste0("Center of mass is outside integration area for rings ",
                     paste(tooShrunk, collapse=", "),
                     ". Setting integration point to closest within the integration area"))
      rsIntegrationPoints[toShrunk] = c(0, rs[-M])[tooShrunk]
    }
  } else if(integrationPointType == "midpoint") {
    rsIntegrationPoints = rsIntegrationPointsMidpoint
  }
  ps = densityFunFinal(rsIntegrationPoints)

  # these weights are for a distribution uniform in radial distance and in angle, marginally:
  if(urban) {
    # ws[1] = rs[1] / rs[M]
    # w2 = (r2 - r1) / r3 * 1 / m2
    # w3 = (r3 - r2) / r3 * 1 / m3
    ws = (rs - c(0, rs[-M])) / rs[M] * 1 / ms

    if(verbose) {
      print("This should be 1 (the sum of the weights), and should all be equal:")
    }
  } else {
    # w1 = r1 / r2 * (99 / 100) + r1 / r3 / 100
    # w2 = (r2 - r1) / (r2 * m2) * (99 / 100) + (r2 - r1) / (r3 * m2) / 100
    # w3 = (r3 - r2) / (r3 * m3) / 100
    annulusWidths = diff(c(0, rs))
    annulusWidthsInner = annulusWidths[1:JInner]
    annulusWidthsOuter = annulusWidths[-(1:JInner)]
    wsInner = annulusWidthsInner * (
      (99 / 100) / (rsInner[JInner] * msInner) +
        (1 / 100) / (rs[M] * msInner) )
    wsOuter = annulusWidthsOuter * (1 / 100) / (rs[M] * msOuter)
    ws = c(wsInner, wsOuter)

    if(verbose) {
      print(paste0("This should be 1 (the sum of the weights), and the first ", JInner, " (and last ", JOuter, ") should be equal:"))
    }
  }
  if(verbose) {
    print(sum(ms*ws))
    print(matrix(ws, nrow=1))
  }

  # pts1 = matrix(c(0, 0), nrow=1)
  # as2 = seq(0, 2*pi, l=m2+1)[-(m2+1)]
  # pts2 = (r1 + r2) / 2 * cbind(cos(as2), sin(as2))
  # as3 = seq(0, 2*pi, l=m3+1)[-(m3+1)] + pi/m3
  # pts3 = (r2 + r3) / 2 * cbind(cos(as3), sin(as3))
  pts = list()
  as = list()
  for(j in 1:M) {
    if(j == 1) {
      as = list(as = 0)
      pts = list(pts = matrix(c(0, 0), nrow=1))
    } else {
      thisas = seq(0, 2*pi, l=ms[j]+1)[-(ms[j]+1)]
      if(j > 1 && j %% 2 == 1) {
        thisas = thisas + pi / ms[j]
      }
      as = c(as, as=list(thisas))
      thisPts = rsIntegrationPoints[j] * cbind(cos(thisas), sin(thisas))
      pts = c(pts, pts=list(thisPts))
    }
  }
  names(as) = paste0("as", 1:M)
  names(pts) = paste0("pts", 1:M)

  if(scalingFactor != 1) {
    rs = rs * scalingFactor
    As = As * scalingFactor^2
    ps = ps / scalingFactor^2
    pts = lapply(pts, function(x){x * scalingFactor})
    densityFunFinalScaled = function(d) {
      (1 / scalingFactor^2) * densityFunFinal(d / scalingFactor)
    }
  } else {
    densityFunFinalScaled = densityFunFinal
  }

  # list(r1=r1, r2=r2, r3=r3, m1=m1, m2=m2, m3=m3,
  #      A1=A1, A2=A2, A3=A3,
  #      w1=w1, w2=w2, w3=w3,
  #      p1=p1, p2=p2, p3=p3,
  #      as2=as2, as3=as3,
  #      pts1=pts1, pts2=pts2, pts3=pts3,
  #      densityFun=densityFunFinal)
  list(rs=rs, ms=ms, As=As, ws=ws, ps=ps,
       pts=pts, as=as, ptRs=rsIntegrationPoints,
       densityFun=densityFunFinalScaled)
}

# construct integration points as well as weights.
# Output: 3 pairs of matrices of dimension nCoordsInStratum x nIntegrationPointsInStratum,
#         each pair contains one urban matrix and one equivalent rural matrix
#   x: x/easting coordinates
#   y: y/northing coordinates
#   w: integration weights
# Input:
#   coords: 2 column matrix of observation easting/northing coordinates
#   urbanVals: vector of observation urbanicity classifications
#   numPointsUrban: number of urban numerical integration points
#   numPointsRural: number of rural numerical integration points
#   scalingFactor: factor by which to scale the jittering distribution.
#                  1 corresponds to standard DHS jittering, larger than 1
#                  corresponds to more jittering than DHS
#   JInnerUrban: number of inner integration rings for urban points
#   JOuterUrban: number of outer integration rings for urban points
#   JInnerRural: number of inner integration rings for rural points
#   JOuterRural: number of outer integration rings for rural points
#   integrationPointType: 'mean' is center of mass, 'midpoint' is the
#                         median angle and median radius within the
#                         integration area
makeAllIntegrationPoints = function(coords, urbanVals,
                                    numPointsUrban=11, numPointsRural=16,
                                    scalingFactor,
                                    JInnerUrban=3, JOuterUrban=0,
                                    JInnerRural=3, JOuterRural=1,
                                    integrationPointType=c("mean", "midpoint"),
                                    adminMap=NULL, nSubAPerPoint=10, nSubRPerPoint=10,
                                    testMode=FALSE) {

  # calculate integration points and weights relative to individual points
  outUrban = getIntegrationPoints(urban=TRUE, numPointsUrban,
                                  scalingFactor,
                                  JInnerUrban, JOuterUrban,
                                  integrationPointType,
                                  verbose=FALSE)

  outRural = getIntegrationPoints(urban=FALSE, numPointsRural,
                                  scalingFactor,
                                  JInnerRural, JOuterRural,
                                  integrationPointType,
                                  verbose=FALSE)

  # concatenate integration points and weights into a single vector
  xsUrbanVec = unlist(sapply(outUrban$pts, function(x) {x[,1]}))
  ysUrbanVec = unlist(sapply(outUrban$pts, function(x) {x[,2]}))
  wsUrbanVec = rep(outUrban$ws, outUrban$ms)
  xsRuralVec = unlist(sapply(outRural$pts, function(x) {x[,1]}))
  ysRuralVec = unlist(sapply(outRural$pts, function(x) {x[,2]}))
  wsRuralVec = rep(outRural$ws, outRural$ms)

  # separate coordinates in urban and rural


  coordsUrban = sf::st_coordinates(coords[urbanVals,])
  coordsRural = sf::st_coordinates(coords[!urbanVals,])
  nUrban = nrow(coordsUrban)
  nRural = nrow(coordsRural)

  # calculate the six matrices
  xUrban = outer(coordsUrban[,1], xsUrbanVec, "+")
  yUrban = outer(coordsUrban[,2], ysUrbanVec, "+")
  wUrban = matrix(wsUrbanVec, ncol=length(wsUrbanVec), nrow=nUrban, byrow=TRUE)
  xRural = outer(coordsRural[,1], xsRuralVec, "+")
  yRural = outer(coordsRural[,2], ysRuralVec, "+")
  wRural = matrix(wsRuralVec, ncol=length(wsRuralVec), nrow=nRural, byrow=TRUE)

  if(!is.null(adminMap)) {
    # if adminMap is input, integration weights will be adjusted using a much
    # finer "sub" integration grid

    # first subset jittered points that are close enough to the border to have
    # a possibility of being adjusted

    # determine if points are within max distance of admin area:
    # first calculate max distance from admin area
    maxUrbanDistance = 2 * scalingFactor
    maxRuralDistance = 10 * scalingFactor
    maxDist = rep(maxRuralDistance, nrow(coords))
    maxDist[urbanVals] = maxUrbanDistance

    adminMap$OBJECTID = 1:length(adminMap$NAME_1)

    temp = sf::st_join(coords, adminMap)

    adminID = temp$OBJECTID

# NOTE : adminMap[adminID,] below returns the same number of polygons
# with the number of points and then st_distance outputs a square matrix


    geomMap = sf::st_geometry(obj = adminMap[adminID,])
    geomMapMulti = sf::st_cast(geomMap, to = 'MULTILINESTRING')
    dists = sf::st_distance(coords, geomMapMulti)


      units(dists) <- NULL # remove the units from dists
      minDistsKM = apply(dists, 2 , min) # minimum distances in kilometers

    # set whether or not to update weights based on distance to admin boundaries
    updateI = as.numeric(minDistsKM) < maxDist
    urbanUpdateI = updateI[urbanVals]
    ruralUpdateI = updateI[!urbanVals]

    # calculate updated weights for the integration points near the borders
    tempCoords = sf::st_coordinates(coords[updateI,])
    tempUrbanVals = urbanVals[updateI]
    #require(fields)

    if(testMode) {
      # in this case, we take only one set of coords that is very close to border,
      # and adjust its weight, saving relevant results for plotting
      tempDists = dists[updateI]
      closeCoords = which.min(tempDists)
      smallTempCoords = matrix(tempCoords[closeCoords,], nrow=1)
      smallTempUrbanVals = tempUrbanVals[closeCoords]

      tempNewWsForPlotting = updateWeightsByAdminArea(coords=smallTempCoords, urbanVals=smallTempUrbanVals,
                                                      adminMap=adminMap,
                                                      integrationPointsUrban=outUrban,
                                                      integrationPointsRural=outRural,
                                                      nSubAPerPoint=nSubAPerPoint,
                                                      nSubRPerPoint=nSubRPerPoint,
                                                      testMode=testMode)

      subPts = tempNewWsForPlotting$subPts
      goodSubPts = tempNewWsForPlotting$goodPts
      subWs = tempNewWsForPlotting$updatedSubWs

      if(smallTempUrbanVals) {
        thisIntPts = outUrban
        thisIntWs = tempNewWsForPlotting$wUrban
      } else {
        thisIntPts = outRural
        thisIntWs = tempNewWsForPlotting$wRural
      }

      return(list(centerCoords=smallTempCoords, isUrban=smallTempUrbanVals,
                  intPts=thisIntPts, intWs=thisIntWs,
                  subPts=subPts, goodSubPts=goodSubPts, subWs=subWs))
    }


    tempCoords = data.frame(x=tempCoords[,1], y=tempCoords[,2])
    tempCoords = sf::st_as_sf(tempCoords, coords=c("x", "y"), crs = sf::st_crs(coords))

    tempNewWs = updateWeightsByAdminArea(coords=tempCoords, urbanVals=tempUrbanVals,
                                         adminMap=adminMap,
                                         integrationPointsUrban=outUrban,
                                         integrationPointsRural=outRural,
                                         nSubAPerPoint=nSubAPerPoint,
                                         nSubRPerPoint=nSubRPerPoint)

    # update the weights with the new values
    wUrban[urbanUpdateI,] = tempNewWs$wUrban
    wRural[ruralUpdateI,] = tempNewWs$wRural
  }

  # return list of all matrices
  list(xUrban=xUrban, yUrban=yUrban, wUrban=wUrban,
       xRural=xRural, yRural=yRural, wRural=wRural)
}

# constructs the sub-integration points for each integration point. Used
# to adjust the weights associated with each integration point.
# Arguments:
# integrationPoints: output of getIntegrationPoints
# adminPoly: polygon of Admin area
# nSubAPerPoint: number of unique angles of sub-integration
#                points within any given integration point area
# nSubRPerPoint: number of unique radii of sub-integration
#                points within any given integration point area
getSubIntegrationPoints = function(integrationPoints, centerCoords=cbind(0, 0),
                                   nSubAPerPoint=10, nSubRPerPoint=10) {
  rs = integrationPoints$rs
  ms = integrationPoints$ms
  As = integrationPoints$As
  ws = integrationPoints$ws
  pts = integrationPoints$pts
  as = integrationPoints$as
  ptRs = integrationPoints$ptRs
  densityFun = integrationPoints$densityFun

  centerCoords = pts[[1]]

  # generates sub-integration points for any point
  getSubIntegrationPointsForOnePoint = function(minR, maxR, minA, maxA, theseCenterCoords) {
    widthR = (maxR - minR)/(nSubRPerPoint)
    widthA = (maxA - minA)/(nSubAPerPoint)

    # calculate radial coordinates of the sub-integration points

    # get angular coordinates of the centers of mass
    if(minA <= maxA) {
      theseAs = seq(minA + widthA/2, maxA - widthA/2, by=widthA)
    } else {
      theseAs = seq(minA + widthA/2 - 2*pi, maxA - widthA/2, by=widthA)
      theseAs[theseAs < 0] = theseAs[theseAs < 0] + 2*pi
    }

    # now get radial centers of mass
    theseRs = seq(minR + widthR/2, maxR - widthR/2, by=widthR)
    rsIntegrationPointsMidpoint = theseRs # midpoint solution

    aDiff = widthA
    shrinkFactor = sqrt(2 * (1 - cos(aDiff))) / aDiff # a scaling factor smaller than one
    shrinkFactor[1] = 1
    rsIntegrationPoints = rsIntegrationPointsMidpoint * shrinkFactor
    tooShrunk = rsIntegrationPoints < (rsIntegrationPointsMidpoint - widthR/2)
    if(any(tooShrunk)) {
      warning(paste0("Center of mass is outside integration area for rings ",
                     paste(tooShrunk, collapse=", "),
                     ". Setting integration point to closest within the integration area"))
      rsIntegrationPoints[toShrunk] = rsIntegrationPointsMidpoint[toShrunk] - widthR/2
    }
    theseRs = rsIntegrationPoints

    thesePointsRadial = fields::make.surface.grid(list(rs=theseRs, as=theseAs))

    # convert to Euclidean coordinates
    thesePointsEuclidean = cbind(thesePointsRadial[,1]*cos(thesePointsRadial[,2]),
                                 thesePointsRadial[,1]*sin(thesePointsRadial[,2]))

    # translate coordinates based on the jittered observation coordinates
    sweep(thesePointsEuclidean, 2, c(centerCoords), "+")
  }

  # for every ring:
  #   for every point:
  #     get sub-integration points
  subWs = list()
  subPts = list()
  for(i in 1:length(rs)) {
    thisIntW = ws[i]
    theseSubWs = rep(thisIntW/(nSubRPerPoint*nSubAPerPoint),
                     each=nSubRPerPoint*nSubAPerPoint*nrow(pts[[i]]))

    theseas = as[[i]]
    theseMinR = ifelse(i==1, 0, rs[i-1])
    theseMaxR = rs[i]
    if(length(theseas) != 1) {
      aWidth = theseas[2] - theseas[1]
    } else {
      aWidth = 2*pi
    }
    theseSubPts = c()

    for(j in 1:length(theseas)) {
      # determine boundaries of this integration area
      thisMinA = theseas[j] - aWidth/2
      thisMaxA = theseas[j] + aWidth/2
      thisMinR = theseMinR
      thisMaxR = theseMaxR

      # obtain sub-integration points
      thisSubPts = getSubIntegrationPointsForOnePoint(minR=thisMinR, maxR=thisMaxR,
                                                      minA=thisMinA, maxA=thisMaxA)
      theseSubPts = rbind(theseSubPts, thisSubPts)
    }

    subPts = c(subPts, list(theseSubPts))
    subWs = c(subWs, list(theseSubWs))
  }

  list(subPts=subPts, subWs=subWs)
}

updateWeightsByAdminArea = function(coords,
                                    urbanVals, adminMap,
                                    integrationPointsUrban,
                                    integrationPointsRural,
                                    nSubAPerPoint=10, nSubRPerPoint=10,
                                    testMode=FALSE) {


  adminMapPoly = adminMap
  # calculate set of typical sub-integration points for urban and rural clusters
  subIntegrationPointsUrban = getSubIntegrationPoints(integrationPoints=integrationPointsUrban,
                                                      nSubAPerPoint=nSubAPerPoint,
                                                      nSubRPerPoint=nSubRPerPoint)

  subIntegrationPointsRural = getSubIntegrationPoints(integrationPoints=integrationPointsRural,
                                                      nSubAPerPoint=nSubAPerPoint,
                                                      nSubRPerPoint=nSubRPerPoint)

  # get admin areas associated with coordinates

  out = sf::st_join(coords, adminMap)
  adminNames = out$NAME_1
  adminIDs = out$OBJECTID

  # for each jittered coordinate:
  #   for each integration point:
  #     get associated sub-integration points
  #     get proportion in correct admin area
  #     update integration point weight
  wsUrban = matrix(nrow=sum(urbanVals), ncol=sum(sapply(integrationPointsUrban$pts, function(x) {nrow(x)})))
  wsRural = matrix(nrow=sum(!urbanVals), ncol=sum(sapply(integrationPointsRural$pts, function(x) {nrow(x)})))
  iUrban = 1
  iRural = 1


  for(i in 1:nrow(coords)) {
    # time1 = proc.time()[3]
    theseCoords = matrix(sf::st_coordinates(coords[i,]), nrow=1)
    thisArea = adminNames[i]
    thisAreaID = adminIDs[i]
    thisPoly = adminMapPoly[thisAreaID,]

    # get sub-integration points
    isUrban = urbanVals[i]
    if(isUrban) {
      thisSubOut = subIntegrationPointsUrban
    } else {
      thisSubOut = subIntegrationPointsRural
    }
    thisSubWs = thisSubOut$subWs
    thisSubPts = thisSubOut$subPts
    thisSubPts = lapply(thisSubPts, function(x) {sweep(x, 2, c(theseCoords), "+")})

    goodAreas = list()
    for (i in 1:length(thisSubPts)){
      coorSet = as.data.frame(thisSubPts[[i]])
      coorSet = stats::setNames(coorSet, c("x", "y"))

      thisSubPtsSP = sf::st_as_sf(coorSet, coords=c("x","y"), crs = sf::st_crs(coords))
      goodAreasTemp = sf::st_join(thisSubPtsSP, thisPoly)
      goodAreas[[i]] = !is.na(goodAreasTemp$OBJECTID)
    }

    # update weights for sub-integration points
    updatedSubWs = thisSubWs
    updatedSubWs = lapply(1:length(updatedSubWs), function(x) {
      temp = updatedSubWs[[x]]
      temp[!goodAreas[[x]]] = 0
      temp
    })

    # sum sub-integration weights to get new (unnormalized) integration weights
    nSubPts = nSubRPerPoint * nSubAPerPoint
    tempWs = lapply(updatedSubWs, function(x) {
      nIntPts = length(x) / nSubPts
      aggIDs = rep(1:nIntPts, each=nSubPts)
      aggregate(x, by=list(aggIDs), FUN=sum)$x
    })
    tempWs = unlist(tempWs)

    # normalize new integration weights to sum to 1
    finalWs = tempWs/sum(tempWs)

    # update weights matrix
    if(isUrban) {
      wsUrban[iUrban,] = finalWs
      iUrban = iUrban + 1
    } else {
      wsRural[iRural,] = finalWs
      iRural = iRural + 1
    }

    # time2 = proc.time()[3]
    # print(paste0("Iteration ", i, "/", nrow(coords), " took ", round(time2-time1, 2), " seconds"))
  }

  if(!testMode) {
    list(wUrban=wsUrban, wRural=wsRural)
  } else {
    list(wUrban=wsUrban, wRural=wsRural,
         subPts=thisSubPts, goodPts=goodAreas, updatedSubWs=updatedSubWs)
  }
}

# Inputs:
#  integrationPointInfo: outfrom from makeAllIntegrationPoints
#  ys: observation vector
#  urbanicity: vector of TRUE/FALSE depending on urbanicity of the observation
#  ns: observation denominator vector
#  spdeMesh: spde triangular basis function mesh object
# Outputs:
#  dat: data.frame containing y, n, east, north, w, urban for each integration point
#  AUrban: (nObs x nUrbanIntegrationPts) x nMesh sparse spatial projection matrix for
#          urban observations. Every nObsUrban rows is the same observation but a new
#          integration point
#  ARural: (nObs x nRuralIntegrationPts) x nMesh sparse spatial projection matrix for
#          rural observations. Every nObsRural rows is the same observation but a new
#          integration point
makeJitterDataForTMB = function(integrationPointInfo, ys, urbanicity, ns, spdeMesh) {
  # first extract the integration point information


  xUrban = integrationPointInfo$xUrban
  yUrban = integrationPointInfo$yUrban
  wUrban = integrationPointInfo$wUrban
  xRural = integrationPointInfo$xRural
  yRural = integrationPointInfo$yRural
  wRural = integrationPointInfo$wRural

  # get the long set of coordinates
  coordsUrban = cbind(c(xUrban), c(yUrban))
  coordsRural = cbind(c(xRural), c(yRural))

  # separate observations by urbanicity
  ysUrban = ys[urbanicity]
  ysRural = ys[!urbanicity]
  nsUrban = ns[urbanicity]
  nsRural = ns[!urbanicity]


  # construct `A' matrices
  AUrban.mesher = fmesher::fm_evaluator(mesh = spdeMesh,
                               loc = coordsUrban)
  AUrban = AUrban.mesher[["proj"]][["A"]]

  ARural.mesher = fmesher::fm_evaluator(mesh = spdeMesh,
                               loc = coordsRural)
  ARural = ARural.mesher[["proj"]][["A"]]


  list(ysUrban=ysUrban, ysRural=ysRural,
       nsUrban=nsUrban, nsRural=nsRural,
       AUrban=AUrban, ARural=ARural)
}
