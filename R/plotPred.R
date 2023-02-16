#' Plots the predictions and the corresponding uncertainty (coefficient of
#' variation), and saves them into the preferred directory as separate .pdf files.
#'
#' @param pred A matrix that is the output of predRes() function.
#' @param predRaster The prediction raster that is constructed by the gridCountry() function.
#' @param admin0 A SpatialPolygonsDataFrame representing the national level (admin0) borders of the country.
#' @param admin1 A SpatialPolygonsDataFrame representing the first level (admin1) subnational borders of the country.
#' @param admin2 A SpatialPolygonsDataFrame representing the second level (admin2) subnational borders of the country.
#' @param rmPoly A number referring to the ID number of the admin2 level polygon that needs to be left uncolored.
#' @param locObs A data frame containing the coordinates of the observation points (DHS locations) in kilometers.
#' @param dir The directory that the plots needs to be saved into.
#' @examples
#' \dontrun{
#' plotPred(x = predcitions, predRaster = predRaster, admin0 = admin0,
#' admin1 = admin1, admin2 = admin2, rmPoly = 160, locObs = locObs,
#' dir = "~/Desktop")
#' }
#' @export
plotPred = function(pred = NULL, predRaster = NULL, admin0 = NULL, admin1 = NULL, admin2 = NULL, rmPoly = NULL, locObs = NULL, dir = NULL){

  proj = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
  admin0_trnsfrmd = sp::spTransform(admin0,proj)
  admin1_trnsfrmd = sp::spTransform(admin1,proj)
  admin2_trnsfrmd = sp::spTransform(admin2,proj)

  idx = 1:raster::ncell(predRaster)
  predCoords = raster::xyFromCell(predRast, idx)
  predCoords = sp::SpatialPoints(predCoords, proj4string=sp::CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))

  uncertainty = (pred[,3]/pred[,1])*100
  uncertainty = raster::setValues(predRast, values = uncertainty, index=idx)

  pred = raster::setValues(predRast, values = pred[,2], index=idx)

  if (rmPoly){
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

  locsPred = raster::xyFromCell(predRast, idx)

  # plotting the predictions
  val = raster::getValues(pred)
  d=data.frame(East = locsPred[,1],
             North = locsPred[,2],val = val)

  ggplot2::ggplot(d, ggplot2::aes(East,North)) +
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
    ggplot2::scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(0, 0.99), na.value="white") +ggplot2::geom_point(data = locObs, color = "red", size=0.001, shape="plus")+
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 2.5, barheight = 25, title = ggplot2::labs("pred."), title.vjust=3) ) +
    ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))

  ggplot2::ggsave("predictions.pdf", path = dir)

  # plotting the uncertainty (coefficient of variation)
  val = raster::getValues(uncertainty)
  d=data.frame(East = locsPred[,1],
               North = locsPred[,2],val = val)


  ggplot2::ggplot(d, ggplot2::aes(East,North)) +
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
    ggplot2::scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(0.68, 233.75), na.value="white") +ggplot2::geom_point(data = locObs, color = "red", size=0.001, shape="plus")+
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 2.5, barheight = 25, title = ggplot2::labs("cv (%)"), title.vjust=3) ) +
    ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))

  ggplot2::ggsave("uncertainty.pdf", path = dir)


}






