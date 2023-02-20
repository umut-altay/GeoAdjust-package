## Resubmission
This is a resubmission. In this version I have:

* Added two references to the description section of the description file. 
  They are preprint papers that explain the methods behind the package.

* Added \value to .Rd files of functions plotPred() and print.res(). Explained 
  the structure of the output.
  
* Explained the output of estimateModel() function better, since its class is 
  the class that is printed by print.res() function.
  
* Fixed the unexecutable code in man/gridCountry.Rd.

* Removed the \dontrun{} command from three functions:
  -convertDegToKM
  -convertDegToMollweide
  -convertKMToDeg
  
* Couldn't remove the \dontrun{} command from other functions because,
  - They either require external data to run (shape files),
  - Or they need INLA package to be installed, which is not on CRAN, 
  - Or they need output from other functions, which either run on INLA, or require 
    external data, or both.
  - \donttest{} would potentially cause errors due to the reasons above.

* Came accross an "rgdal is required for spTransform methods" error in winbuilder
  and fixed it by adding rgdal package to suggests. 


## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

*   checking installed package size ... NOTE
    installed size is 12.6Mb
    sub-directories of 1Mb or more:
    libs  12.4Mb
    
    The package currently contains only the files that sustain its 
    functionality. There are no extra files (images, vignette, etc.)

*   checking package dependencies ... NOTE
    Package suggested but not available for checking: 'INLA'
    
    The package INLA is not on CRAN, therefore it is registered under additional
    repositories and suggests in the description file. The functions which 
    require INLA in the package are written in a way that they either instruct the 
    user on how to install INLA from its repository, or if it is already 
    installed, they attach the namespace.
    
*   checking CRAN incoming feasibility ... [11s] NOTE
    Maintainer: 'The package maintainer <altayumut.ua@gmail.com>'

    The note indicating the email address of the package maintainer.
    
