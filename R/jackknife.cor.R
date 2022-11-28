#' @title Delete-d Jackknife Estimate for Correlation between Two Variables
#'
#' @description This function creates jackknife samples from the data by
#' sequentially removing *d* observations from the data,
#' calculates correlation between the two variables using the jackknife samples
#' and estimates the jackknife correlation coefficients, bias standard error,
#' standard error and confidence intervals.
#'
#' @param x,y Numeric vectors of equal length
#' @param d Number of observations to be deleted from data to make jackknife samples
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
#' @importFrom stats coefficients cor lm qnorm
#' @importFrom utils combn
#' @export
#' @examples
#' j.cor <- jackknife.cor(cars$speed, cars$dist, d = 2)
#' j.cor$jackknife.summary
#' j.cor$biased_cor
#'
jackknife.cor <- function(x, y, d){
  n <- length(x)

  if((n*2)^d > 2e+08){message("This may take more time. Please wait...")}

  cmb <- combn(n, d) # Row indexes to be eliminated for jackknife
  N <- ncol(cmb)     # Total number of jackknife samples
  jk <- numeric(N)   # A numeric vector to collect the jackknife estimates

  for (i in 1:N) {
    j <- cmb[,i]
    jk[i] <- cor(x[-j], y[-j])
  }

  theta_hat <- cor(x, y) # Biased Correlation Estimate
  theta_dot_hat <- mean(jk) # Mean of correlation estimates of jackknife samples
  bias <- (n-d) * (theta_dot_hat-theta_hat) #Bias
  est <- theta_hat - bias
  jack_se <- sqrt(((n-d)/d)  *  mean((jk-theta_hat)^2)) # Jackknife standard error
  jack_ci_lower <- est-(qnorm(0.975)*jack_se)
  jack_ci_upper <- est+(qnorm(0.975)*jack_se)

  jackknife.summary <- data.frame(Estimate = est,
                                  bias = bias,
                                  se = jack_se,
                                  t = est/jack_se,
                                  ci.lower = jack_ci_lower,
                                  ci.upper = jack_ci_upper)

  return(list(jackknife.summary = jackknife.summary,
              biased_cor = theta_hat,
              Jackknife.samples.coeff = jk))
}

