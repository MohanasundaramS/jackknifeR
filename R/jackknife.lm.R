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
#' @param numCores Number of processors to be used
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
#' @importFrom stats lm.fit model.frame model.matrix model.response hatvalues
#' @examples
#' ## library(jackknifeR)
#' jk <- jackknife.lm(mpg ~ wt + hp, d = 2, data = mtcars, numCores = 2)
#' summary(jk)
#' @export
jackknife.lm <- function(formula, d = 1, data, conf = 0.95, numCores = detectCores()) {
  cl <- match.call()

  # Precompute model components
  mf <- model.frame(formula, data)
  tm <- terms(mf)
  X <- model.matrix(tm, mf)
  var_names <- colnames(X)
  p <- length(var_names)

  # Define statistic function
  coef_fun <- function(data_sub) {
    mf_sub <- model.frame(tm, data_sub, na.action = na.omit)
    X_sub <- model.matrix(tm, mf_sub)
    qr_X <- qr(X_sub)
    if (qr_X$rank < p) return(rep(NA, p))
    fit <- lm.fit(X_sub, model.response(mf_sub))
    coefs <- setNames(rep(NA, p), var_names)
    coefs[colnames(X_sub)] <- fit$coefficients
    coefs
  }

  # Call core jackknife function
  result <- jackknife(
    statistic = coef_fun,
    d = d,
    data = data,
    conf = conf,
    numCores = numCores,
    weight = FALSE
  )

  result$call <- cl
  result$formula <- formula
  result
}
