#' Prints the output of estimateModel() function.
#'
#' @param x A list containing the model estimation output, returned by
#' estimateModel() function.
#' @param ... not used
#' @return Prints the model estimation results of class res as a table that
#' shows the estimated model parameters and the corresponding 95% credible
#' interval lengths.
#' @examples
#' \donttest{
#' path1 <- system.file("extdata", "exampleInputData.rda", package = "GeoAdjust")
#' path2 <- system.file("extdata", "exampleMesh.rda", package = "GeoAdjust")
#' load(path1)
#' load(path2)
#' nNodes = exampleMesh[['n']]
#' results <- estimateModel(data = exampleInputData, nNodes = nNodes,
#' options = list(random = 1, covariates = 1), priors = list(beta = c(0,1),
#' range = 114, USpatial = 1, alphaSpatial = 0.05, UNugget = 1, alphaNug = 0.05))
#' print(results)
#' }
#' @export
print.res <- function(x, ...){
  p1 = lengths(x[["res"]]) # number of rows
  likelihood = x[["likelihood"]]
  cat("GeoAdjust::estimateModel()", "\n")
  cat("----------------------------------\n")
  cat("Likelihood :    ", likelihood, sep= "      ")
  cat("\n")
  cat("----------------------------------\n")
  p2 = c("parameter", "estimate")
  cat(cat(p2, sep= " "), "95% CI length", sep = "    ")
  cat("\n")
  cat(cat(paste0("range"), format(round(as.numeric(x[["res"]][["estimates"]][[1]]), 4), nsmall = 4),sep ="     "),
      format(round(as.numeric(x[["res"]][["CI_Length"]][[1]]), 4), nsmall = 4), sep ="    ")
  cat("\n")
  cat(cat(paste0("sigma"), format(round(as.numeric(x[["res"]][["estimates"]][[2]]), 4), nsmall = 4), sep ="     "),
      format(round(as.numeric(x[["res"]][["CI_Length"]][[2]]), 4), nsmall = 4),sep = "      ")
  cat("\n")

  if(likelihood=="Gaussian"){
    cat(cat(paste0("nuggetSD"), format(round(as.numeric(x[["res"]][["estimates"]][[3]]), 4), nsmall = 4),sep ="  "),
        format(round(as.numeric(x[["res"]][["CI_Length"]][[3]]), 4), nsmall = 4), sep ="      ")
    cat("\n")

    cat(cat(paste0("intercept"), format(round(as.numeric(x[["res"]][["estimates"]][[4]]), 4), nsmall = 4),sep =" "),
          format(round(as.numeric(x[["res"]][["CI_Length"]][[4]]), 4), nsmall = 4), sep ="      ")
         cat("\n")

    if(p1[[1]]>4){
    for (i in 5:p1[[1]]){
      estimated = format(round(as.numeric(x[["res"]][["estimates"]][[i]]), 4), nsmall = 4)
      CI = format(round(as.numeric(x[["res"]][["CI_Length"]][[i]]), 4), nsmall = 4)
      cat(cat(paste0(x[["res"]][["parameters"]][[i]]), estimated, sep ="     "),
          CI, sep ="      ")
      cat("\n")
      cat("----------------------------------\n")
    }
    }

  }else{
    cat(cat(paste0("intercept"), format(round(as.numeric(x[["res"]][["estimates"]][[3]]), 4), nsmall = 4),sep =" "),
         format(round(as.numeric(x[["res"]][["CI_Length"]][[3]]), 4), nsmall = 4), sep ="      ")
         cat("\n")
         if(p1[[1]]>3){
    for (i in 4:p1[[1]]){
      estimated = format(round(as.numeric(x[["res"]][["estimates"]][[i]]), 4), nsmall = 4)
      CI = format(round(as.numeric(x[["res"]][["CI_Length"]][[i]]), 4), nsmall = 4)
      cat(cat(paste0(x[["res"]][["parameters"]][[i]]), estimated, sep ="     "),
          CI, sep ="      ")
      cat("\n")
      cat("----------------------------------\n")
    }
         }
  }

  # cat(paste0("beta0"), format(round(as.numeric(x[["res"]][["estimates"]][[3]]), 4), nsmall = 4),
  #     format(round(as.numeric(x[["res"]][["CI_Length"]][[3]]), 4), nsmall = 4), sep ="     ")
  # cat("\n")
  # if(length(x[["res"]][["parameters"]])>3){
  #  for (i in 4:p1[[1]]){
  # estimated = format(round(as.numeric(x[["res"]][["estimates"]][[i]]), 4), nsmall = 4)
  # CI = format(round(as.numeric(x[["res"]][["CI_Length"]][[i]]), 4), nsmall = 4)
  # cat(cat(paste0(x[["res"]][["parameters"]][[i]]), estimated, sep ="     "),
  #     CI, sep ="      ")
  # cat("\n")
  # cat("----------------------------------\n")
  #  }
  # } else{
  #   cat("----------------------------------\n")
  # }

}








