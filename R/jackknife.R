#' @title Delete-d Jackknife for Estimates
#'
#' @description This function creates jackknife samples from the data by
#' sequentially removing *d* observations from the data,
#' and calculates the estimates by the specified function and its bias,
#' standard error, and confidence intervals.
#' @param statistic a function returning a vector of estimates to be passed to jackknife
#' @param d Number of observations to be deleted from data to make jackknife samples. The default is 1 (for delete-1 jackknife).
#' @param data Data frame with dependent and independent independent variables specified in the formula
#' @param conf Confidence level, a positive number < 1. The default is 0.95.
#' @return A list containing a summary data frame of jackknife estimates
#'    with bias, standard error. t-statistics, and confidence intervals,
#'    estimate for the original sample and a data frame with
#'    estimates for jackknife samples.
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
#' @seealso [jackknife.lm()] which is used for jackknifing in linear regression.
#' @importFrom stats coefficients qnorm
#' @importFrom utils combn
#' @export
#' @examples
#' ## library(jackknifeR)
#' fn <- function(data){
#'    mod <- lm(speed~dist, data = data)
#'    return(coef(mod))}
#' jkn <- jackknife(statistic = fn, d = 2, data = cars)
#' jkn$jackknife.summary
#'
#'

jackknife <- function(statistic, d = 1,  data, conf = 0.95){

  n <- nrow(data)
  if(is.numeric(conf)==FALSE||conf>1||conf<0) stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
  if((n*ncol(data))^d > 9e+06) stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
  if((n*ncol(data))^d > 1e+04){message("This may take more time. Please wait...")}

  cmb <- combn(n, d) # Row indexes to be eliminated for jackknife
  N <- ncol(cmb)     # Total number of jackknife samples

  stat <- statistic(data) # Estimate by defined function

  # A data frame to collect the jackknife estimates
  jk <- as.data.frame(matrix(nrow = N, ncol = length(stat)))
  colnames(jk) <- names(stat)

  for (i in 1:N) {
    j <- cmb[,i]
    jk[i,] <- statistic(data[-j,])
  }

  theta_hat <- stat # Regression coefficient estimates
  theta_dot_hat <- colMeans(jk) # Mean of regression coefficient estimates of jackknife samples
  bias <- (n-d) * (theta_dot_hat-theta_hat) # Bias
  est <- theta_hat-bias
  jack_se <- sqrt((n-d)/d  *  rowMeans(apply(jk, 1, function(x) (x-theta_hat)^2))) # Jackknife standard error
  jack_ci_lower <- est-(qnorm(0.5+(conf/2))*jack_se)
  jack_ci_upper <- est+(qnorm(0.5+(conf/2))*jack_se)

  jackknife.summary <- data.frame(Estimate = est,
                                  bias = bias,
                                  se = jack_se,
                                  t = est/jack_se,
                                  ci.lower = jack_ci_lower,
                                  ci.upper = jack_ci_upper)

  return(list(jackknife.summary = jackknife.summary,
              original.estimate = stat,
              Jackknife.samples.est = jk))
}

