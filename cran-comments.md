## Resubmission
This is a resubmission. In this version I have:

* Removed all \dontrun{} commands.

* Created two small size SpatialPolygonsDataFrame objects representing national 
  and sub-national level borders of an artificially created example country 
  with only four sub-national administrative areas. 
  
* Created an example data frame with only 10 observations, mimicing the 
  Demographic and Health Surveys program (DHS) survey data. These three files 
  are stored under"geoData.rda" file.
  
* Reprojected the national level SpatialPolygonsDataFrame into UTM: zone 37 
  coordinates and stored as "adm0UTM37.rda". 
  
* Created example outputs of gridCountry(), meshCountry(), prepareInput() and
  predRes() functions. These are called exampleGrid, exampleMesh, 
  exampleInputData and examplePredictionResults, respectively.
  
* Stored all example files under inst/extdata subfolder.

* Didn't store an example output from estimateModel() function due to its 
  size (77 MB).
  
* Put the examples that cause more than 10 sec check time in \donttest{}
  command.
  
## R CMD check results

0 errors | 0 warnings | 3 notes 

* checking installed package size ... NOTE
  installed size is 13.2Mb
  sub-directories of 1Mb or more:
  libs  12.4Mb
  
  The package contains small example files under inst/extdata. 
  subfolder. The files are constructed in smallest size possible.
  
* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'INLA'

  The code is structured in a way that the functions either instruct the user 
  about how to install INLA, or if it is already installed, it gets attached to
  the NAMESPACE.

* checking CRAN incoming feasibility ... [9s] NOTE
  Maintainer: 'Umut Altay <altayumut.ua@gmail.com>'

  A note indicating the name of the package maintainer.


  
