## Resubmission
This is a resubmission due to major updates/improvements in the package functionality:

1.) The package is migrated from soon to be retired R-spatial packages to the 
    up-to-date ones. Currently the package uses functions from sf and terra to 
    operate on spatial objects and rasters.

2.) Dependency on INLA is removed by using functions from fmesher instead.

3.) All examples are updated according to the new functionality.

4.) The small example data files that are part of the package are remade accordingly,
    in a minimum possible size.
    
5.) As a result of the updates, the exported functions convertDegToKM(), 
    convertDegToMollweide(), convertKMToDeg(), and covMatern() are removed, since
    they are no more needed. The internal function modeSPDEJitter is also no more
    needed, and removed as well.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking installed package size ... NOTE
  installed size is 13.1Mb
  sub-directories of 1Mb or more:
  libs  12.4Mb
  
  The package contains small example files under inst/extdata. 
  subfolder. The files are constructed in smallest size possible.
  



  
