#' Plots the predictions and the corresponding uncertainty (coefficient of
#' variation)
#'
#' @param pred A matrix that is the output of predRes() function.
#' @param predRaster The prediction raster that is constructed by the gridCountry() function.
#' @param admin0 A SpatialPolygonsDataFrame representing the national level (admin0) borders of the country.
#' @param admin1 A SpatialPolygonsDataFrame representing the first level (admin1) subnational borders of the country.
#' @param admin2 A SpatialPolygonsDataFrame representing the second level (admin2) subnational borders of the country.
#' @param rmPoly A number referring to the ID number of the admin2 level polygon that needs to be left uncolored. It can be set to NULL as well.
#' @param locObs A data frame containing the coordinates of the observation points (DHS locations) in kilometers.
#' @return A list of two ggplot objects. One of them (ggPred) shows the median predictions and the other one (ggUncertainty) shows the
#' corresponding coefficient of variations across the country, respectively.
#' @examples
#' path1 <- system.file("extdata", "examplePredictionResults.rda", package = "GeoAdjust")
#' path2 <- system.file("extdata", "exampleGrid.rda", package = "GeoAdjust")
#' path3 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' load(path2)
#' load(path3)
#' plots = plotPred(pred = examplePredictionResults,
#' predRaster = exampleGrid[["predRast"]], admin0 = adm0,
#' admin1 = adm1, admin2 = NULL, rmPoly = NULL,
#' locObs = data.frame(East = surveyData$east, North = surveyData$north))
#' @export
plotPred = function(pred = NULL, predRaster = NULL, admin0 = NULL, admin1 = NULL, admin2 = NULL, rmPoly = NULL, locObs = NULL){

  proj = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
  admin0_trnsfrmd = sp::spTransform(admin0,proj)
  admin1_trnsfrmd = sp::spTransform(admin1,proj)
  if(!is.null(admin2)){
    admin2_trnsfrmd = sp::spTransform(admin2,proj)
  }
  idx = 1:raster::ncell(predRaster)
  predCoords = raster::xyFromCell(predRaster, idx)
  predCoords = sp::SpatialPoints(predCoords, proj4string=sp::CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))

  uncertainty = (pred[,3]/pred[,1])*100
  uncertainty = raster::setValues(predRaster, values = uncertainty, index=idx)

  pred = raster::setValues(predRaster, values = pred[,2], index=idx)

  if (!is.null(rmPoly)){
    polyg = admin2_trnsfrmd[rmPoly,] # the lake

    # find which points are inside the polygon that needs to be removed, and assign NA to them
    pointInPolygon = rgeos::gWithin(predCoords, polyg, byid=TRUE)

    inRmPoly = which(pointInPolygon == TRUE)

    pred[inRmPoly] = NA
    uncertainty[inRmPoly] = NA
  }

  pred = raster::mask(raster::crop(pred, raster::extent(admin0_trnsfrmd)), admin0_trnsfrmd, snap = 'out')
  uncertainty = raster::mask(raster::crop(uncertainty, raster::extent(admin0_trnsfrmd)), admin0_trnsfrmd, snap = 'out')


  dfCountry <- ggplot2::fortify(admin1_trnsfrmd, region = "NAME_1")

  locsPred = raster::xyFromCell(predRaster, idx)

  # plotting the predictions
  val = raster::getValues(pred)
  d=data.frame(East = locsPred[,1],
               North = locsPred[,2],val = val)

  ggPred = ggplot2::ggplot(d, ggplot2::aes(East,North)) +
    ggplot2::geom_raster(ggplot2::aes(fill=val)) + ggplot2::theme_bw() +
    ggplot2::geom_path(data = dfCountry, ggplot2::aes(long,lat, group = group),colour = "black", inherit.aes = FALSE)+
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 35), axis.text.x = ggplot2::element_text(size = 35)) +
    ggplot2::theme(axis.title.x=ggplot2::element_text(size = ggplot2::rel(3))) + ggplot2::theme(axis.title.y=ggplot2::element_text(size = ggplot2::rel(3)))+
    ggplot2::theme(legend.title = ggplot2::element_text(size = ggplot2::rel(3))) + ggplot2::coord_fixed() +
    ggplot2::xlab("Easting (km)") +
    ggplot2::ylab("Northing (km)")  +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank()) +
    ggplot2::theme(legend.text=ggplot2::element_text(size=35))+
    ggplot2::scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(min(val), max(val)), na.value="white") +ggplot2::geom_point(data = locObs, color = "red", size=0.001, shape="plus")+
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 2.5, barheight = 25, title = ggplot2::labs("pred."), title.vjust=3) ) +
    ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))

  # plotting the uncertainty (coefficient of variation)
  val = raster::getValues(uncertainty)
  d=data.frame(East = locsPred[,1],
               North = locsPred[,2],val = val)


  ggUncertainty = ggplot2::ggplot(d, ggplot2::aes(East,North)) +
    ggplot2::geom_raster(ggplot2::aes(fill=val)) + ggplot2::theme_bw() +
    ggplot2::geom_path(data = dfCountry, ggplot2::aes(long,lat, group = group),colour = "black", inherit.aes = FALSE)+
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 35), axis.text.x = ggplot2::element_text(size = 35)) +
    ggplot2::theme(axis.title.x=ggplot2::element_text(size = ggplot2::rel(3))) + ggplot2::theme(axis.title.y=ggplot2::element_text(size = ggplot2::rel(3)))+
    ggplot2::theme(legend.title = ggplot2::element_text(size = ggplot2::rel(3))) + ggplot2::coord_fixed() +
    ggplot2::xlab("Easting (km)") +
    ggplot2::ylab("Northing (km)")  +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank()) +
    ggplot2::theme(legend.text=ggplot2::element_text(size=35))+
    ggplot2::scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(min(val), max(val)), na.value="white") +ggplot2::geom_point(data = locObs, color = "red", size=0.001, shape="plus")+
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 2.5, barheight = 25, title = ggplot2::labs("cv (%)"), title.vjust=3) ) +
    ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))

  return(list(ggPred = ggPred, ggUncertainty = ggUncertainty))
}





