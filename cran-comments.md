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
    
