#' Bayesian Compound Poisson-Normal (bCPN) Regression Model
#'
#' Fits a Bayesian Compound Poisson-Normal regression model using JAGS.
#' Supports intercept-only or single binary factor models.
#'
#' @param formula A formula specifying the model. Must be either intercept-only
#'   (`y ~ 1`) or contain exactly one binary factor (`y ~ factor`).
#' @param data A data frame containing the variables in the model.
#' @param n_iter Integer. Number of MCMC iterations to run after burn-in.
#'   Default is 20,000.
#' @param burn_in Integer. Number of MCMC burn-in iterations. Default is 5,000.
#' @param n_chains Integer. Number of MCMC chains. Default is 1.
#' @param k_max Numeric scalar or NULL. Maximum truncation value for the Poisson
#'   counts. Must be between 10 and 100. Default is 10.
#' @param priors A named list specifying prior distributions as JAGS model
#'   strings. Defaults are:
#'   \itemize{
#'     \item \code{beta_0}: prior for intercept (default:
#'           \code{"dnorm(0, 0.000001)"})
#'     \item \code{beta_1}: prior for factor levels (default:
#'           \code{"dnorm(0, 0.00001)"})
#'     \item \code{mu}: prior for normal mean (default:
#'           \code{"dnorm(0, 0.000000001)"})
#'     \item \code{sigma}: prior for normal sd (default:
#'           \code{"dgamma(0.01, 0.01)"})
#'   }
#'
#' @return An \code{mcmc.list} object from the \code{coda} package containing
#'   posterior samples for monitored parameters (including \code{mu},
#'   \code{sigma}, \code{Expected}, and \code{Diff} if applicable).
#'
#' @details
#' The model assumes the response follows a Compound Poisson-Normal
#' distribution where the count variable \(X\) is Poisson-distributed with a
#' log-linear model for its rate parameter, and the observations \(y\) are
#' normal with mean \(\eqn{\mu} X\) and variance \(\eqn{\sigma}^2 / X\) (with a
#' large variance if \(X=0\)).
#'
#' The parameter \code{k_max} truncates the Poisson rate at a maximum value to
#' avoid numerical issues in JAGS.
#'
#' Only models with an intercept or a single binary factor are supported.
#'
#' @examples
#' testdata <- simulate_cpn_data()
#' fit <- bcpn(y ~ x1, data = testdata, n_iter = 20000)
#' summary(fit)
#'
#' @importFrom rjags jags.model coda.samples
#' @importFrom stats as.formula model.frame model.response terms update
#' @export
bcpn <- function(formula,
                 data,
                 n_iter = 20000,
                 burn_in = 5000,
                 n_chains = 1,
                 k_max = 10,
                 priors = list(
                   beta_0 = "dnorm(0, 0.000001)",
                   beta_1 = "dnorm(0, 0.00001)",
                   mu = "dnorm(0, 0.000000001)",
                   sigma = "dgamma(0.01, 0.01)"
                 )) {

  null_or_numeric_1l <- function(x) {
    is.null(x) || (is.numeric(x) && length(x) == 1L)
  }

  stopifnot(
    inherits(formula, "formula"),
    is.data.frame(data),
    null_or_numeric_1l(k_max)
  )

  if (k_max < 10 || k_max > 100) {
    stop("k_max should be between 10 and 100")
  }

  formula <- stats::as.formula(formula, env = parent.frame())
  mf <- stats::model.frame(formula = formula, data = data)
  y <- stats::model.response(mf)

  terms_rhs <- attr(terms(formula), "term.labels")

  if (length(terms_rhs) == 0) {
    ## Intercept-only model
    N <- length(y) #nolint
    model_cpn <- sprintf("
    model {
      beta_0 ~ %s
      mu ~ %s
      sigma ~ %s
      tau <- 1 / (sigma * sigma)

      lambda <- exp(beta_0)
      Expected <- mu * lambda

      for (i in 1:N) {
        X[i] ~ dpois(censored_lambda[i])
        censored_lambda[i] <- min(lambda, k_max)
        tau_Y[i] <- ifelse(X[i] > 0, tau / X[i], 1e10)
        y[i] ~ dnorm(mu * X[i], tau_Y[i])
      }
    }
    ", priors$beta_0, priors$mu, priors$sigma)

    data_list <- list(N = N, y = y, k_max = k_max)
    inits <- function() list(beta_0 = 0, mu = mean(y), sigma = 1)
    monitor_params <- c("mu", "sigma", "Expected")

  } else if (length(terms_rhs) == 1) {
    ## Single binary factor model
    trt_var <- terms_rhs[1]
    trt <- factor(data[[trt_var]])
    if (nlevels(trt) != 2) {
      stop("Treatment factor must have exactly 2 levels.")
    }

    N <- length(y) #nolint

    model_cpn <- sprintf("
    model {
      for (k in 1:2) {
        beta_1[k] ~ %s
      }
      beta_0 ~ %s
      mu ~ %s
      sigma ~ %s
      tau <- 1 / (sigma * sigma)

      for (k in 1:2) {
        reslambda[k] <- exp(beta_0 + beta_1[k])
        Expected[k] <- mu * reslambda[k]
      }

      Diff <- Expected[2] - Expected[1]

      for (i in 1:N) {
        X[i] ~ dpois(censored_lambda[i])
        censored_lambda[i] <- min(lambda[i], k_max)
        log(lambda[i]) <- beta_0 + beta_1[trt[i]]
        tau_Y[i] <- ifelse(X[i] > 0, tau / X[i], 1e10)
        y[i] ~ dnorm(mu * X[i], tau_Y[i])
      }
    }
    ", priors$beta_1, priors$beta_0, priors$mu, priors$sigma)

    data_list <- list(
      N = N,
      y = y,
      trt = as.numeric(trt),
      k_max = k_max
    )
    inits <- function() {
      list(
        beta_0 = 0,
        beta_1 = c(0, 0),
        mu = mean(y),
        sigma = 1
      )
    }

    monitor_params <- c("mu", "sigma", "Expected", "Diff")

  } else {
    stop("Only intercept-only or single binary factor models are supported.")
  }

  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' required.")
  }

  model <- rjags::jags.model(textConnection(model_cpn),
                             data = data_list,
                             inits = inits
                             , n.chains = n_chains)
  stats::update(model, burn_in)

  rjags::coda.samples(model, variable.names = monitor_params, n.iter = n_iter)
}


#' Plot Posterior Distributions of Expected Means from bCPN Model
#'
#' Plots posterior density estimates of the \code{Expected} parameter(s)
#' from an MCMC output of the \code{bcpn} function.
#'
#' @param x An \code{mcmc.list} object containing samples from \code{bcpn}.
#' @param param_names Character vector of parameter names to plot. Defaults to
#'   all \code{Expected} parameters found.
#' @param col Colors for the density lines.
#'    Default is \code{c("skyblue", "salmon")}.
#' @param lwd Line width for density lines. Default is 2.
#' @param add_medians Logical indicating whether to add median vertical lines.
#'    Default TRUE.
#' @param ... Additional graphical parameters passed to \code{plot}.
#'
#' @importFrom graphics abline legend lines text
#' @export
plot.mcmc.list <- function(x,
                           param_names = NULL,
                           col = c("skyblue", "salmon"),
                           lwd = 2,
                           add_medians = TRUE,
                           ...) {

  mcmc_mat <- as.matrix(x)

  # Detect both scalar and vector 'Expected'
  available_params <- grep("^Expected(\\[|$)", colnames(mcmc_mat), value = TRUE)

  if (length(available_params) == 0) {
    stop("No 'Expected' parameters found in MCMC output.")
  }

  if (is.null(param_names)) {
    param_names <- available_params
  } else {
    if (!all(param_names %in% available_params)) {
      stop("Some param_names not found in MCMC output.")
    }
  }

  n_params <- length(param_names)

  if (n_params == 1) {
    # Single Expected parameter (scalar or Expected[1])
    param1 <- mcmc_mat[, param_names[1]]
    dens1 <- stats::density(param1)
    median1 <- stats::median(param1)

    plot(dens1, col = col[1], lwd = lwd,
         xlab = "Expected Value",
         ylab = "Density",
         main = "Posterior Distribution of the Mean",
         ...)

    if (add_medians) {
      graphics::abline(v = median1, col = col[1], lwd = 2, lty = 2)
      graphics::text(median1, 0, labels = round(median1, 2),
                     col = col[1], pos = 4, srt = 45)
    }

    graphics::legend("topright", legend = param_names[1],
                     col = col[1], lwd = lwd, bty = "n")

  } else if (n_params >= 2) {
    # At least two Expected parameters (e.g., two-level factor model)
    param1 <- mcmc_mat[, param_names[1]]
    param2 <- mcmc_mat[, param_names[2]]

    dens1 <- stats::density(param1)
    dens2 <- stats::density(param2)

    median1 <- stats::median(param1)
    median2 <- stats::median(param2)

    ymax <- max(dens1$y, dens2$y)

    plot(dens1, col = col[1], lwd = lwd,
         xlab = "Expected Value",
         ylab = "Density",
         main = "Posterior Distributions of the Means",
         ylim = c(0, ymax * 1.1),
         ...)

    lines(dens2, col = col[2], lwd = lwd)

    if (add_medians) {
      graphics::abline(v = median1, col = col[1], lwd = 2, lty = 2)
      graphics::abline(v = median2, col = col[2], lwd = 2, lty = 2)

      graphics::text(median1, 0, labels = round(median1, 2),
                     col = col[1], pos = 4, srt = 45)
      graphics::text(median2, 0, labels = round(median2, 2),
                     col = col[2], pos = 4, srt = 45)
    }

    graphics::legend("topright", legend = param_names[1:2],
                     col = col[1:2], lwd = lwd, bty = "n")
  } else {
    stop("Unexpected number of 'Expected' parameters found.")
  }
}
