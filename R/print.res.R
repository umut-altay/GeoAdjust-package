#' Prints the output of estimateModel() function.
#'
#' @param x A list containing the model estimation output, returned by
#' estimateModel() function.
#' @param ... ignored.
#' @return Prints the model estimation results of class res as a table that
#' shows the estimated model parameters and the corresponding 95% credible
#' interval lengths.
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
  cat(paste0("range"), format(round(as.numeric(x[["res"]][["estimates"]][[1]]), 4), nsmall = 4),
      format(round(as.numeric(x[["res"]][["CI_Length"]][[1]]), 4), nsmall = 4), sep ="     ")
  cat("\n")
  cat(cat(paste0("sigma"), format(round(as.numeric(x[["res"]][["estimates"]][[2]]), 4), nsmall = 4), sep ="     "),
      format(round(as.numeric(x[["res"]][["CI_Length"]][[2]]), 4), nsmall = 4),sep = "      ")
  cat("\n")
  cat(paste0("beta0"), format(round(as.numeric(x[["res"]][["estimates"]][[3]]), 4), nsmall = 4),
      format(round(as.numeric(x[["res"]][["CI_Length"]][[3]]), 4), nsmall = 4), sep ="     ")
  cat("\n")
   for (i in 4:p1[[1]]){
  estimated = format(round(as.numeric(x[["res"]][["estimates"]][[i]]), 4), nsmall = 4)
  CI = format(round(as.numeric(x[["res"]][["CI_Length"]][[i]]), 4), nsmall = 4)
  cat(cat(paste0(x[["res"]][["parameters"]][[i]]), estimated, sep ="     "),
      CI, sep ="      ")
  cat("\n")
  cat("----------------------------------\n")
   }

}








