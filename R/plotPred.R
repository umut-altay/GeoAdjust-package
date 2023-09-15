#' Plots the predictions and the corresponding uncertainty (coefficient of
#' variation)
#'
#' @param pred A matrix that is the output of predRes() function.
#' @param predRaster The prediction raster that is constructed by the
#' gridCountry() function.
#' @param admin0 An sf class MULTIPOLYGON representing the national level (admin0)
#' borders of the country.
#' @param admin1 An sf class MULTIPOLYGON representing the first level (admin1)
#' subnational borders of the country.
#' @param admin2 An sf class MULTIPOLYGON representing the second level (admin2)
#' subnational borders of the country.
#' @param rmPoly A number referring to the ID number of the admin2 level polygon
#' that needs to be left uncolored. It can be set to NULL as well.
#' @param target_crs A projection string representing the desired coordinate
#' reference system according to which the maps will be created.
#' @return A list of two ggplot objects. One of them (ggPred) shows the median
#' predictions and the other one (ggUncertainty) shows the
#' corresponding coefficient of variations across the country, respectively.
#' @examples
#' \donttest{
#' path1 <- system.file("extdata", "examplePredictionResults.rda", package = "GeoAdjust")
#' path2 <- system.file("extdata", "geoData.rda", package = "GeoAdjust")
#' load(path1)
#' load(path2)
#' crs_KM = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
#' exampleGrid <- gridCountry(admin0 = adm0, res = 5, target_crs = crs_KM)
#' plots = plotPred(pred = examplePredictionResults,
#' predRaster = exampleGrid[["predRast"]], admin0 = adm0,
#' admin1 = adm1, target_crs = crs_KM)
#' }
#' @export
plotPred = function(pred = NULL, predRaster = NULL, admin0 = NULL, admin1 = NULL, admin2 = NULL, rmPoly = NULL, target_crs=NULL){

  admin0_trnsfrmd = sf::st_transform(admin0, target_crs)
  admin1_trnsfrmd = sf::st_transform(admin1, target_crs)

  if(!is.null(admin2)){
    admin2_trnsfrmd = sf::st_transform(admin2, target_crs)
  }

  idx = 1:terra::ncell(predRaster)

  predCoords = terra::xyFromCell(predRaster, idx)
  predCoords = data.frame(east = predCoords[,1], north = predCoords[,2])
  predCoords = sf::st_as_sf(predCoords, coords=c("east","north"), crs = target_crs)

  uncertainty = (pred[,3]/pred[,1])*100
  uncertainty = terra::setValues(predRaster, values = uncertainty)

  pred = terra::setValues(predRaster, values = pred[,2])

  if (!is.null(rmPoly)){
    polyg = admin2_trnsfrmd[rmPoly,] # the polygon to be left uncolored

    # find which points are inside the polygon that needs to be removed, and assign NA to them

    pointInPolygon = sf::st_within(predCoords, polyg, sparse = FALSE, prepared = TRUE)

    inRmPoly = which(pointInPolygon == TRUE)

    pred[inRmPoly] = NA
    uncertainty[inRmPoly] = NA
  }

  pred = terra::mask(pred, admin0_trnsfrmd)
  uncertainty = terra::mask(uncertainty, admin0_trnsfrmd)

  dfCountry <- ggplot2::fortify(admin1_trnsfrmd, region = "NAME_1")

  # plotting the predictions
  val = terra::values(pred)
  val = val[,1]
  d=data.frame(East = sf::st_coordinates(predCoords)[,1],
               North = sf::st_coordinates(predCoords)[,2],
               val=val)
  colnames(d) = c("East", "North", "val")

  ggPred =ggplot2::ggplot(d, ggplot2::aes(East,North)) +
    ggplot2::geom_raster(ggplot2::aes(fill=val)) + ggplot2::theme_bw() +
    ggplot2::geom_sf(data=dfCountry, fill=NA, colour = "lightgrey",inherit.aes = FALSE)+ggplot2::coord_sf(datum=sf::st_crs(target_crs))+
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10), axis.text.x = ggplot2::element_text(size = 10)) +
    ggplot2::theme(axis.title.x=ggplot2::element_text(size = ggplot2::rel(3))) + ggplot2::theme(axis.title.y=ggplot2::element_text(size = ggplot2::rel(3)))+
    ggplot2::theme(legend.title = ggplot2::element_text(size = ggplot2::rel(3))) + #ggplot2::coord_fixed() +
    ggplot2::xlab("Easting (km)") +
    ggplot2::ylab("Northing (km)")  +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank()) +
    ggplot2::theme(legend.text=ggplot2::element_text(size=35))+
    ggplot2::scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(min(val, na.rm=T), max(val, na.rm=T)), na.value="white") +#ggplot2::geom_point(data = locObs, color = "red", size=0.001, shape="plus")+
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 2.5, barheight = 25, title = ggplot2::labs("pred."), title.vjust=3) ) +
    ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))


  # plotting the uncertainty (coefficient of variation)
  val = terra::values(uncertainty)
  val = val[,1]
  d=data.frame(East = sf::st_coordinates(predCoords)[,1],
               North = sf::st_coordinates(predCoords)[,2],
               val=val)
  colnames(d) = c("East", "North", "val")

  ggUncertainty = ggplot2::ggplot(d, ggplot2::aes(East,North)) +
    ggplot2::geom_raster(ggplot2::aes(fill=val)) + ggplot2::theme_bw() +
    ggplot2::geom_sf(data=dfCountry, fill=NA, colour = "lightgrey",inherit.aes = FALSE)+ggplot2::coord_sf(datum=sf::st_crs(target_crs))+
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10), axis.text.x = ggplot2::element_text(size = 10)) +
    ggplot2::theme(axis.title.x=ggplot2::element_text(size = ggplot2::rel(3))) + ggplot2::theme(axis.title.y=ggplot2::element_text(size = ggplot2::rel(3)))+
    ggplot2::theme(legend.title = ggplot2::element_text(size = ggplot2::rel(3))) + #ggplot2::coord_fixed() +
    ggplot2::xlab("Easting (km)") +
    ggplot2::ylab("Northing (km)")  +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank()) +
    ggplot2::theme(legend.text=ggplot2::element_text(size=35))+
    ggplot2::scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(min(val), max(val)), na.value="white") +#ggplot2::geom_point(data = locObs, color = "red", size=0.001, shape="plus")+
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 2.5, barheight = 25, title = ggplot2::labs("cv (%)"), title.vjust=3) ) +
    ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))

  return(list(ggPred = ggPred, ggUncertainty = ggUncertainty))
}





