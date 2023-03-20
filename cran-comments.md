## Resubmission
This is a resubmission to fix a critical error that has been just discovered:

Line 205 of the file compute.cpp was operating on wrong indexing as beta[i+1].
This caused the optimization to throw an error, because the code was trying to 
run on non-existing betas (covariate effect sizes). The indexing should have 
been beta[i], instead. 

This was a critical error preventing the estimateModel() function from working
when there are covariates in the model. It is fixed now. 

* In addition, there are two small changes which improve the functionality:

  1.) An argument called n.sims added to the estimateModel() function. The 
  argument controls the number of samples that should be drawn for each parameter
  in the model. This was previously hard coded as 10.000. Now it is possible to 
  set it to a lower number to speed up the computation.
  
  2.) plotPred() function was plotting two ggplot objects to the plots pane. Now
  the function outputs them within a list, so that the package users can decide
  which object they want to be plotted. 

* Examples are updated according to the small changes mentioned above.

  
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

* checking CRAN incoming feasibility ... [11s] NOTE
  Maintainer: 'Umut Altay <altayumut.ua@gmail.com>'

  A note indicating the name of the package maintainer.


  
