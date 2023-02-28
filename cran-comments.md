## Resubmission
This is a resubmission. In this version I have:

* Removed all  \dontrun{} commands.

* Created two small size SpatialPolygonsDataFrame objects representing national 
  and sub-national level borders of a square shaped example country with only 
  four sub-national administrative areas. Created an example data frame with 
  only 10 observations, mimicing the Demographic and Health Surveys program 
  (DHS) survey data. These three files are stored under"geoData.rda" file.
  
* Reprojected the national level SpatialPolygonsDataFrame into UTM: zone 37 
  coordinates and stored as as "adm0UTM37.rda". 
    
  
* Created example outputs of gridCountry(), meshCountry(), prepareInput() and
  predRes() functions.
  
* Didn't store the output from estimateModel() function due to its size (77 MB). 


    
