% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meshCountry.R
\name{meshCountry}
\alias{meshCountry}
\title{Creates a constrained refined Delaunay triangulation mesh based on the country borders.}
\usage{
meshCountry(admin0 = NULL, max.edge = NULL, offset = NULL)
}
\arguments{
\item{admin0}{A SpatialPolygonsDataFrame representing the country borders,
in UTM:zone 37 coordinate system.}

\item{max.edge}{A vector of two values. The first and the second elements of
the vector represent the largest allowed triangle lengths for the inner and outer mesh, respectively.}

\item{offset}{A value representing the extension distance for the inla.mesh.2d object}
}
\value{
A constrained refined Delaunay triangulation mesh created based on the country borders.
}
\description{
Creates a constrained refined Delaunay triangulation mesh based on the country borders.
}
\examples{
if(requireNamespace("INLA")){
path1 <- system.file("extdata", "adm0UTM37.rda", package = "GeoAdjust")
load(path1)
mesh.s <- meshCountry(admin0 = adm0UTM37, max.edge = c(25, 50), offset = -.08)
}
}
