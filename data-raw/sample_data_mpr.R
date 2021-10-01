# Sample Data for MPR-SIC Method ------------------------------------------
# * Generate heteroscedastic data -----------------------------------------
beta_true <- c(0, 1, 0.5, 0.5, 1, 0.5, 1, 0, 0, 0, 0,
               0, 0)
alpha_true <- c(0, 0.5, 1, 0.5, 1, 0, 0, 0.5, 1, 0, 0,
                0, 0)
n <- 500
set.seed(100)
  p <- length(beta_true) - 1

  x <- cbind(
    rep(1, n),
    matrix(rnorm(n * p), n, p)
  )
  # Generate error vector
  phi_true <- as.vector(exp(x %*% alpha_true))
  err_vec <- rnorm(n = n, mean = 0, sd = (sqrt(phi_true)))

  y <- (x %*% beta_true) + err_vec


sample_data_mpr <- list(
  "x" = x,
  "y" = y,
  "beta_true" = beta_true,
  "alpha_true" = alpha_true
)
sample_data_mpr$x <- sample_data_mpr$x[, -1] # remove column of 1's

# save(sample_data_mpr,
#      file = here::here("data", "sample_data_mpr.RData"))
