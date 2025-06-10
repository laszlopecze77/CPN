#' Simulate Data from a Compound Poisson-Normal Model
#'
#' This function generates synthetic data based on a Compound
#' Poisson-Normal (CPN) model. The number of events for each
#' observation is drawn from a Poisson distribution, and the
#' outcome is the sum of normally distributed values for each event.
#'
#' The predictors include a binary categorical variable (`x1`) and
#' a continuous variable (`x2`). A linear model with coefficients
#' `beta` is used to model the log of the Poisson rate.
#'
#' @param n Integer. Number of observations to simulate. Default is 100.
#' @param beta Numeric vector. Coefficients for the linear predictor,
#'   including intercept. Default is \code{c(0.5, -0.3, 0.7)}. Must
#'   match the number of columns in the model matrix: intercept, x1B, x2.
#' @param mu Numeric. Mean of the Normal distribution for each event.
#'   Default is 1.
#' @param sigma Numeric. Standard deviation of the Normal distribution
#'   for each event. Default is 2.
#'
#' @return A \code{data.frame} with \code{n} rows and 3 columns:
#'   \describe{
#'     \item{y}{Numeric response variable generated from the CPN model.}
#'     \item{x1}{Categorical predictor with levels "A" and "B".}
#'     \item{x2}{Continuous predictor drawn from a standard normal
#'       distribution.}
#'   }
#'
#' @examples
#' set.seed(42)
#' simulated_data <- simulate_cpn_data(n = 200)
#' head(simulated_data)
#'
#' @export
simulate_cpn_data <- function(n = 100,
                              beta = c(0.5, -0.3, 0.7),
                              mu = 1,
                              sigma = 2) {
  set.seed(123)  # For reproducibility

  # Simulate predictors
  x1 <- factor(sample(c("A", "B"), size = n, replace = TRUE))  # Categorical
  x2 <- rnorm(n)  # Continuous

  # Create model matrix (includes intercept and dummy variable for x1)
  X <- model.matrix(~ x1 + x2) # nolint

  # Compute linear predictor and Poisson rates
  eta <- X %*% beta
  lambda <- exp(eta)

  # Simulate response variable
  y <- numeric(n)
  for (i in 1:n) {
    k <- stats::rpois(1, lambda[i])
    if (k > 0) {
      y[i] <- sum(stats::rnorm(k, mean = mu, sd = sigma))
    } else {
      y[i] <- 0
    }
  }

  # Return as a data.frame
  data.frame(y = y, x1 = x1, x2 = x2)
}
