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
#' @param weight Logical, TRUE for weighted jackknife standard error of regression estimates. Default weight = FALSE
#' @param hat_values Vector of hat values (leverages) from the model. Required if `weight = TRUE
#' @param residuals Vector of residuals from the model. Required if `weight = TRUE`.
#' @param X Model matrix. Required if `weight = TRUE`.
#' @param p Number of predictors in the model. Required if `weight = TRUE`.
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
#' @importFrom future plan multisession
#' @importFrom doFuture registerDoFuture
#' @importFrom foreach foreach %dopar%
#' @importFrom utils combn
#' @importFrom stats model.matrix lm qt setNames
#' @importFrom stats complete.cases na.omit pt symnum terms
#' @importFrom future.apply future_apply
#' @export
#' @examples
#' library(future)
#' plan(multisession)  # Initialize once per session
#' # For linear regression coefficients
#' jk_results <- jackknife(
#' statistic = function(sub_data) coef(lm(mpg ~ wt + hp, data = sub_data)),
#' d = 2,
#' data = mtcars,
#' conf = 0.95, numCores = 2)
#' print(jk_results)

jackknife <- function(statistic, d = 1, data, conf = 0.95, numCores = detectCores(),
                      weight = FALSE, hat_values = NULL, residuals = NULL, X = NULL, p = NULL) {

  # Initialize connection cleanup
  on.exit({
    future::plan(future::sequential)  # Reset to sequential plan
    closeAllConnections()            # Close any remaining connections
  })

  idx <- NULL
  n <- nrow(data)
  nd <- choose(n, d)

  # Common setup
  indices <- utils::combn(n, d, simplify = FALSE)

  # Parallel setup
  registerDoFuture()
  plan(multisession, workers = numCores)

  if (weight) {
    # ---- Weighted Jackknife Core ----
    # Compute sum of leverages for deleted observations
    sum_hat <- future.apply::future_sapply(indices, function(idx) sum(hat_values[idx]))

    # Compute subsample coefficients
    jk <- foreach::foreach(idx = indices, .combine = rbind) %dopar% {
      statistic(data[-idx, ])
    }

    # Remove invalid subsamples
    valid <- complete.cases(jk)
    jk <- jk[valid, ]
    sum_hat <- sum_hat[valid]
    nd_valid <- nrow(jk)

    # Compute delete-d scaling
    scaling_factor <- (n - d) / (d * nd)

    # ---- Critical Fix: Leverage-adjusted residuals ----
    # Compute HC3-adjusted residuals
    adjusted_residuals <- residuals^2 / (1 - hat_values)^2  # HC3 adjustment

    # Compute Hinkley's weighted variance with scaling
    sum_Rxx <- crossprod(X, diag(adjusted_residuals) %*% X)
    D0_inv <- solve(crossprod(X))
    V_w <- scaling_factor * (n / (n - p)) * D0_inv %*% sum_Rxx %*% D0_inv
    se <- sqrt(diag(V_w))

    # Compute bias using original estimates
    theta_hat <- statistic(data)
    theta_dot <- colMeans(jk)
    bias <- ((n - d)/d) * (theta_dot - theta_hat)
  } else {
    # ---- Original Unweighted Method (Unchanged) ----
    jk <- foreach::foreach(idx = indices, .combine = rbind) %dopar% {
      res <- statistic(data[-idx, ])
      matrix(res, nrow = 1)  # Ensure output is matrix
    }
    valid <- complete.cases(jk)
    jk <- jk[valid, , drop = FALSE]  # Keep as matrix
    nd_valid <- nrow(jk)

    theta_hat <- statistic(data)
    theta_dot <- colMeans(jk)
    bias <- ((n - d)/d) * (theta_dot - theta_hat)
    ssq <- colSums(sweep(jk, 2, theta_dot, "-")^2)
    se <- sqrt(((n - d)/(d * nd_valid)) * ssq)
  }

  # Explicitly stop clusters at end
  future::plan(future::sequential)
  closeAllConnections()

  # Common post-processing
  df <- n - d - ifelse(weight, p, 0)
  t_val <- qt((1 + conf)/2, df = df)
  ci <- theta_hat - bias + t_val * se %o% c(-1, 1)

  structure(
    list(
      estimates = theta_hat,
      bias = bias,
      df = df,
      se = se,
      ci = ci,
      conf = conf,
      d = d,
      call = match.call()
    ),
    class = "jackknife"
  )
}

#' @export
print.jackknife <- function(x, ...) {
  cat("Delete-", x$d, " Jackknife Results\n", sep = "")
  cat("Confidence Level: ", x$conf, "\n\n")
  print(data.frame(
    Original_Estimate = x$estimates,
    Bias = x$bias,
    SE = x$se,
    CI_Lower = x$ci[,1],
    CI_Upper = x$ci[,2]
  ))
}

#' @export
summary.jackknife <- function(object, ...) {
  # Calculate significance stars and p-values
  t_stat <- object$estimates / object$se
  p_value <- 2 * pt(-abs(t_stat), df = object$df)
  stars <- symnum(p_value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))

  # Create formatted data frame
  summary_df <- data.frame(
    Estimate = object$estimates,
    Bias = object$bias,
    Std.Error = object$se,
    `t value` = t_stat,
    `Pr(>|t|)` = format.pval(p_value),
    ` ` = stars,
    CI.Lower = object$ci[,1],
    CI.Upper = object$ci[,2],
    check.names = FALSE
  )

  structure(
    list(
      summary = summary_df,
      d = object$d,
      conf = object$conf,
      call = object$call
    ),
    class = "summary.jackknife"
  )
}

#' @export
print.summary.jackknife <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nDelete-", x$d, " Jackknife Results (", x$conf*100, "% CI)\n", sep="")
  print(x$summary, digits = 4, na.print = "-")
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}
