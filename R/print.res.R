#' Prints the output of estimateModel() function.
#'
#' @param est Distance matrix.
#' @return Prints the model parameter estimates and the corresponding credible interval lengths.
#' @examples
#' \dontrun{
#' print.res(est = est)
#' }
#' @export
print.res <- function(est){
  p1 = lengths(est[["res"]]) # number of rows
  likelihood = est[["likelihood"]]
  cat("GeoAdjust::estimateModel()", "\n")
  cat("----------------------------------\n")
  cat("Likelihood :    ", likelihood, sep= "      ")
  cat("\n")
  cat("----------------------------------\n")
  p2 = c("parameter", "estimate")
  cat(cat(p2, sep= " "), "95% CI length", sep = "    ")
  cat("\n")
  cat(paste0("range"), format(round(as.numeric(est[["res"]][["estimates"]][[1]]), 4), nsmall = 4),
      format(round(as.numeric(est[["res"]][["CI_Length"]][[1]]), 4), nsmall = 4), sep ="     ")
  cat("\n")
  cat(cat(paste0("sigma"), format(round(as.numeric(est[["res"]][["estimates"]][[2]]), 4), nsmall = 4), sep ="     "),
      format(round(as.numeric(est[["res"]][["CI_Length"]][[2]]), 4), nsmall = 4),sep = "      ")
  cat("\n")
  cat(paste0("beta0"), format(round(as.numeric(est[["res"]][["estimates"]][[3]]), 4), nsmall = 4),
      format(round(as.numeric(est[["res"]][["CI_Length"]][[3]]), 4), nsmall = 4), sep ="     ")
  cat("\n")
   for (i in 4:p1[[1]]){
  estimated = format(round(as.numeric(est[["res"]][["estimates"]][[i]]), 4), nsmall = 4)
  CI = format(round(as.numeric(est[["res"]][["CI_Length"]][[i]]), 4), nsmall = 4)
  cat(cat(paste0(est[["res"]][["parameters"]][[i]]), estimated, sep ="     "),
      CI, sep ="      ")
  cat("\n")
  cat("----------------------------------\n")
   }

}








