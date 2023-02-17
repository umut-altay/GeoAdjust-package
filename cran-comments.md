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
    repositories in the description file. The functions which requires INLA in
    the package are written in a way that they either warn the user to install
    INLA from its relevant repository, or if it is already installed, they 
    attach the namespace.
    
*   checking CRAN incoming feasibility ... [10s] NOTE
    Maintainer: 'The package maintainer <geir-arne.fuglstad@ntnu.no>'
    
    The package author (Umut Altay) is a PhD student, whose email account
    will expire after graduation. Geir-Arne Fuglstad is the authors' supervisor.
    The package maintainer will be Geir-Arne Fuglstad since he has a permanent 
    email account.
    
    
    
