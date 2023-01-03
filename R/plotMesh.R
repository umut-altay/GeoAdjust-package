#' Create a mesh based on the country borders.
#'
#' @param admin0 A SpatialPolygonsDataFrame representing the country borders.
#' @param mesh The triangulation mesh
#' @param name The preferred name for the plot
#' @param dir The preferred directory that the plot will be saved into
#' @return Saves a pdf file showing the country (admin0) borders plotted onto the mesh, to the preferred directory
#' @examples
#' plotMesh(admin0 = NULL, mesh = NULL, name = NULL, dir = NULL)
#' @export
#' @import rgeos
#' @import sp
plotMesh = function(admin0, mesh, name, dir){   #INLA is needed to plot the mesh
  dir = dir
  name = name
  proj = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"
  admin0_trnsfrmd = spTransform(admin0,proj)
  simplified = gSimplify(admin0_trnsfrmd, tol=.01, topologyPreserve=TRUE)

  path = file.path(dir, name)
  pdf(paste0(path), width = 25, height = 14) # Open a new pdf file
  plot(mesh)
  plot(simplified, border = "red", lwd=3, add = TRUE)
  dev.off() # Close the file
}




