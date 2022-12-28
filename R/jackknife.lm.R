#' @title Delete-d Jackknife Estimate for Linear Regression
#'
#' @description This function creates jackknife samples from the data by
#' sequentially removing *d* observations from the data,
#' fits models linear regression model using the jackknife samples
#' as specified in the formula and estimates the jackknife coefficients
#' bias standard error, standard error and confidence intervals.
#'
#' @param formula Simple or multiple  linear regression formula with dependent and independent variables
#' @param d Number of observations to be deleted from data to make jackknife samples. The default is 1 (for delete-1 jackknife).
#' @param data Data frame with dependent and independent independent variables specified in the formula
#' @param conf Confidence level, a positive number < 1. The default is 0.95.
#' @return A list containing a summary data frame of jackknife estimates
#'    with bias, standard error. t-statistics, and confidence intervals,
#'    linear regression model of original data and a data frame with
#'    coefficient estimates of jackknife samples.
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
#' @seealso [lm()] which is used for linear regression.
#' @importFrom stats coefficients lm qnorm
#' @importFrom utils combn
#' @export
#' @examples
#' ## library(jackknifeR)
#' j.lm <- jackknife.lm(dist~speed, d = 2, data = cars)
#' j.lm$jackknife.summary
#' summary(j.lm$lm_mod)$coefficients
#'
jackknife.lm <- function(formula, d = 1,  data, conf = 0.95){

  n <- nrow(data)
  if(is.numeric(conf)==FALSE||conf>1||conf<0) stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
  if((n*ncol(data))^d > 9e+07) stop("The number of jackknife sub-samples will be huge")
  if((n*ncol(data))^d > 1e+04){message("This may take more time. Please wait...")}

  cmb <- combn(n, d) # Row indexes to be eliminated for jackknife
  N <- ncol(cmb)     # Total number of jackknife samples

  lm_mod <- lm(formula, data = data) # Linear regression estimate

  # A data frame to collect the jackknife estimates
  jk <- as.data.frame(matrix(nrow = N, ncol = length(coefficients(lm_mod))))
  colnames(jk) <- names(coefficients(lm_mod))

  for (i in 1:N) {
    j <- cmb[,i]
    mod <- lm(formula, data = data[-j,])
    jk[i,] <- coefficients(mod)
  }

  theta_hat <- coefficients(lm_mod) # Regression coefficient estimates
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

  jk.r <- list(jackknife.summary = jackknife.summary,
               d = d,
               conf.level = conf,
               stat = formula,
               n.jack = N,
               original.estimate = theta_hat,
               lm_mod = lm_mod,
               Jackknife.samples.est = jk)
  class(jk.r) <- "jk"

  return(jk.r)
}

