set.seed(123)
# Function to simulate a Compound Poisson-Normal dataset
simulate_cpn <- function(n, lambda, mu, sigma) {
  counts <- rpois(n, lambda)  # Poisson counts
  data <- sapply(counts, function(k) sum(rnorm(k, mean = mu, sd = sigma)))
  data
}

##### Negative mu - Example 1 #####
set.seed(123)

sim_data1 <- simulate_cpn(n = 100, lambda = 1.3, mu = -400, sigma = 400)
sim_data2 <- simulate_cpn(n = 100, lambda = 1, mu = -400, sigma = 400)

# Combine data and group labels
sim_data_pooled <- c(sim_data1, sim_data2)
group <- factor(c(rep(1, length(sim_data1)), rep(2, length(sim_data2))))
df_sim <- data.frame(y = sim_data_pooled, group = group)

# Fit model
cpn_fit <- cpn(y ~ group, data = df_sim)

library(emmeans)
test_that(
  "emmeans output matches manual calculation for cpn model",
  {
    # Run emmeans on link and response (rate) scales
    emm_link <- emmeans(cpn_fit, ~group)
    emm_resp <- emmeans(cpn_fit, ~group, type = "response")

    # Create reference grid and design matrix
    grid <- ref_grid(cpn_fit)@grid
    X <- model.matrix(delete.response(terms(cpn_fit)), grid) # nolint

    # Extract relevant quantities
    bhat <- cpn_fit$coefficients[
      setdiff(names(cpn_fit$coefficients), c("mu", "sigma"))
    ]
    eta_hat <- as.vector(X[, names(bhat), drop = FALSE] %*% bhat)

    V_beta <- vcov(cpn_fit)[names(bhat), names(bhat)]   # nolint
    var_eta <- sapply(seq_len(nrow(X)), function(i) {
      xi <- X[i, , drop = FALSE]
      as.numeric(xi %*% V_beta %*% t(xi))
    })
    se_eta <- sqrt(var_eta)

    z <- qnorm(0.975)
    ci_link_lower <- eta_hat - z * se_eta
    ci_link_upper <- eta_hat + z * se_eta

    # Rate scale (exp(eta)), matching emmeans(..., type = "response")
    response_hat <- exp(eta_hat)
    se_response <- sqrt((exp(eta_hat))^2 * var_eta)
    ci_resp_lower <- response_hat - z * se_response
    ci_resp_upper <- response_hat + z * se_response

    # Extract emmeans summaries
    summ_link <- summary(emm_link)
    summ_resp <- summary(emm_resp)

    # Set tolerance
    tol <- 0.1

    for (i in seq_along(grid$group)) {
      # Link scale checks
      expect_equal(summ_link$emmean[i], eta_hat[i], tolerance = tol)
      expect_equal(summ_link$SE[i], se_eta[i], tolerance = tol)
      expect_equal(summ_link$asymp.LCL[i], ci_link_lower[i], tolerance = tol)
      expect_equal(summ_link$asymp.UCL[i], ci_link_upper[i], tolerance = tol)

      # Rate scale checks
      expect_equal(summ_resp$response[i], response_hat[i], tolerance = tol)
      expect_equal(summ_resp$SE[i], se_response[i], tolerance = tol)
      expect_equal(summ_resp$asymp.LCL[i], ci_resp_lower[i], tolerance = tol)
      expect_equal(summ_resp$asymp.UCL[i], ci_resp_upper[i], tolerance = tol)
    }
  }
)
