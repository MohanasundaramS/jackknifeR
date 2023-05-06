#' @title Delete-d Jackknife Estimate for Correlation between Two Variables
#'
#' @description This function creates jackknife samples from the data by
#' sequentially removing *d* observations from the data,
#' calculates correlation between the two variables using the jackknife samples
#' and estimates the jackknife correlation coefficients, bias standard error,
#' standard error and confidence intervals.
#'
#' @param x,y Numeric vectors of equal length
#' @param d Number of observations to be deleted from data to make jackknife samples. The default is 1 (for delete-1 jackknife).
#' @param conf Confidence level, a positive number < 1. The default is 0.95.
#' @return A list containing a summary data frame of jackknife correlation
#'    coefficient estimates with bias, standard error. t-statistics,
#'    and confidence intervals,correlation estimate of original data and
#'    a data frame with correlation estimates of individual jackknife samples.
#'
#' @references Quenouille, M. H. (1956). Notes on Bias in Estimation.
#' *Biometrika*, *43*(3/4), 353-360.
#' \doi{10.2307/2332914}
#' @references Tukey, J. W. (1958). Bias and Confidence in Not-quite Large Samples.
#' *Annals of Mathematical Statistics*, *29*(2), 614-623.
#' \doi{10.1214/aoms/1177706647}
#' @references Shi, X. (1988). A note on the delete-d jackknife variance estimators.
#' *Statistics & Probability Letters*, *6*(5), 341-347.
#' \doi{10.1016/0167-7152(88)90011-9}
#' @seealso [cor()] which is used to estimate correlation coefficient.
#' @importFrom stats coefficients cor qnorm
#' @importFrom utils combn
#' @export
#' @examples
#' ## library(jackknifeR)
#' j.cor <- jackknife.cor(cars$speed, cars$dist, d = 2)
#' summary(j.cor)
#'
jackknife.cor <- function(data, d = 1, conf = 0.95){
  cl <- match.call()
  fn <- function(data){
    return(cor(data[,1], data[,2]))
  }
  j.cor <- jackknife(statistic = fn, d = d, data =  data, conf = 0.95)
  j.cor$call <- cl
  return(j.cor)
}

