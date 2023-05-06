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
#' @param numCores Number of processors to be used
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
#' @importFrom stats coef qt as.formula
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach %do%
#' @export
#' @examples
#' ## library(jackknifeR)
#' fn <- function(data){
#'    mod <- lm(speed~dist, data = data)
#'    return(coef(mod))}
#' jkn <- jackknife(statistic = fn, d = 2, data = cars, numCores= 2)
#' jkn
#'
#'

jackknife <- function(statistic, d = 1,  data, conf = 0.95, numCores = detectCores()) {
  n <- nrow(data)
  p <- ncol(data)
  npd <- (n * p) ^ d
  qt_df <- n - d

  if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) {
    stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
  }

  if (npd > 9e+12) {
    stop("The number of jackknife sub-samples will be huge")
  }

  if (npd > 1e+05) {
    message("This may take more time. Please wait...")
  }

  fn <- function(data){
    statistic(data)
  }

  st <- fn(data) # Estimate by defined function

  indices <- lapply(1:n, function(i) seq(n)[-i])

  idx = NULL

  cl <- makeCluster(numCores)
  jk <- foreach(idx = indices, .packages = 'parallel') %do% {
    fn(data[unlist(idx),])
  }

  jk <- do.call(rbind, jk)

  theta_hat <- st # Regression coefficient estimates
  theta_dot_hat <- colMeans(jk) # Mean of regression coefficient estimates of jackknife samples
  bias <- (n - d) * (theta_dot_hat - theta_hat) # Bias
  est <- theta_hat - bias

  se_d <- colMeans(sweep(jk, 2, theta_hat) ^ 2, na.rm = TRUE)
  jack_se <- sqrt((n - d) / d * se_d) # Jackknife standard error
  t_val <- qt(p = (1 - conf) / 2, df = qt_df, lower.tail = FALSE)
  jack_ci_lower <- est - t_val * jack_se
  jack_ci_upper <- est + t_val * jack_se

  jackknife.summary <- data.frame(Estimate = est,
                                  bias = bias,
                                  se = jack_se,
                                  t = est / jack_se,
                                  ci.lower = jack_ci_lower,
                                  ci.upper = jack_ci_upper)
  jk.r <- list(call = match.call(),
               jackknife.summary = jackknife.summary,
               d = d,
               conf.level = conf,
               stat = substitute(statistic),
               n.jack = n ^ d,
               original.estimate = theta_hat,
               Jackknife.samples.est = jk)
  class(jk.r) <- "jk"

  stopCluster(cl)
  return(jk.r)
}

#' @export
print.jk <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nConfidence Level: ",
      paste(x$conf.level, sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Jackknife Summary:\n")
  print(format(x$jackknife.summary, digits = digits),
        print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @export
summary.jk <- function (object, ...){
  if(!inherits(object,"jk")){
    stop("Error: The object is not of class 'jk'.")
  }

  ans <- list(call = object$call,
              stat = object$stat,
              d = object$d,
              n.jack = object$n.jack,
              conf.level = object$conf.level,
              jackknife.summary = object$jackknife.summary,
              original.estimate = object$original.estimate,
              Jackknife.samples.est = object$Jackknife.samples.est)

  class(ans) <- "summary.jk"
  ans
}

#' @export
print.summary.jk <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat(cat("\nCall: "),
      deparse(x$call), "\n", sep = "")
  cat("\nConfidence Level: ",
      paste(x$conf.level, sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nd: ",
      paste(x$d, sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Jackknife Summary:\n")
  print(format(x$jackknife.summary, digits = digits),
        print.gap = 2L, quote = FALSE)
  cat("\n")
  cat("Estimate from Original Sample:\n")
  print(format(x$original.estimate, digits = digits),
        print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

