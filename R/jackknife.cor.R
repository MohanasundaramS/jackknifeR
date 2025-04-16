#' @title Delete-d Jackknife Estimate for Correlation between Two Variables
#'
#' @description This function creates jackknife samples from the data by sequentially removing *d* observations, calculates the correlation, and estimates bias, standard error, and confidence intervals.
#'
#' @param data A data frame with two numeric columns.
#' @param d Number of observations to delete (default: 1).
#' @param conf Confidence level (default: 0.95).
#' @param numCores Number of processors (default: `detectCores()`).
#' @return A list of class "jackknife" containing estimates, bias, standard error, and confidence intervals.
#' @references Quenouille (1956), Tukey (1958), Shi (1988).
#' @seealso [cor()], [jackknife()]
#' @importFrom parallel detectCores
#' @importFrom stats cor
#' @export
#' @examples
#' j.cor <- jackknife.cor(cars, d = 2, numCores = 2)
#' summary(j.cor)
jackknife.cor <- function(data, d = 1, conf = 0.95, numCores = parallel::detectCores()) {
  # Input validation
  if (!is.data.frame(data) || ncol(data) != 2 || !all(vapply(data, is.numeric, logical(1)))) {
    stop("data must be a 2-column numeric data frame")
  }

  cl <- match.call()

  # Force scalar output with proper formatting
  result <- jackknife(
    statistic = function(subdata) {
      cor_val <- cor(subdata[, 1], subdata[, 2])
      unname(cor_val)  # Remove any names
    },
    d = d,
    data = data,
    conf = conf,
    numCores = numCores
  )

  result$call <- cl
  return(result)
}
