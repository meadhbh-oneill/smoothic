


# All functions to pull from ----------------------------------------------
# Smooth Generalized Normal Distribution ----------------------------------
# Optimize kappa as a single parameter (constant)
# * Generate Data ---------------------------------------------------------
kappa_to_nu <- function(kappa,
                        kappa_omega) {
  log(kappa - kappa_omega)
}

nu_to_kappa <- function(nu,
                        kappa_omega) {
  exp(nu) + kappa_omega
}

nu_to_kappa_matrix <- function(nu,
                               x3,
                               kappa_omega) {
  exp(x3 %*% nu) + kappa_omega
}

# Simulating data ---------------------------------------------------------
# Density -----------------------------------------------------------------
f_density <- function(x,
                      kappa,
                      tau) {
  c_tilde_val <- c_tilde_integrate_smoothgnd(
    kappa = kappa,
    tau = tau
  )
  density_gnd_approx_standard(
    x = x,
    kappa_shape = kappa,
    tau = tau,
    c_tilde = c_tilde_val
  )
}


# CDF ---------------------------------------------------------------------
F_cdf <- function(from,
                  to,
                  by,
                  kappa,
                  tau) {
  x <- seq(from, to, by = by)
  f <- f_density(
    x = x,
    kappa = kappa,
    tau = tau
  )
  f <- f / sum(f * by)
  F <- cumsum(f * by)
  cbind(x = x, F = F)
}

# Inverse CDF -------------------------------------------------------------
F2x <- function(u, F) { # converts uniform (0, 1) to values from distribution
  svec <- 1 - u
  F <- as.data.frame(F)
  x <- F$x
  x <- c(min(x), x, max(x))
  S <- c(1, 1 - F$F, 0)
  xS <- cbind(x, S)[length(x):1, ]
  alls <- c(svec, xS[, 2])
  inset <- rep(c(0, 1), c(length(svec), length(xS[, 2])))[order(alls)]
  index <- cumsum(inset)[inset == 0][order(order(svec))]
  index <- ifelse(index == 0, 1, index) + (svec %in% S)
  xS[index, 1]
}

# Quantile function ----
F_cdf_full <- function(from, # can have other mu and s values
                       to,
                       by,
                       kappa,
                       tau,
                       mu,
                       s) {
  x <- seq(from, to, by = by)
  f <- smooth_gnd_density(
    x = x,
    mu = mu,
    s_scale = s,
    kappa = kappa,
    tau = tau,
    list_general = list_families$smoothgnd$general,
    method_c_tilde = "integrate"
  )
  f <- f / sum(f * by)
  F <- cumsum(f * by)
  cbind(x = x, F = F)
}

qgndnorm <- function(p,
                     mu,
                     s,
                     kappa,
                     tau) {
  F <- F_cdf_full(
    from = -20, # cdf
    to = 20,
    by = 0.001,
    kappa = kappa,
    tau = tau,
    mu = mu,
    s = s
  )

  F2x(u = p,
      F = F)
}

# Simulate from distribution ----------------------------------------------
data_sim <- function(n,
                     beta_true,
                     alpha_true,
                     kappa_true,
                     tau,
                     from = -20,
                     to = 20,
                     by = 0.001) {
  stopifnot(length(beta_true) == length(alpha_true))

  p <- length(beta_true) - 1

  x <- cbind(
    rep(1, n),
    matrix(rnorm(n * p), n, p) # generate X as normal, can also be binary
  )

  mu_i <- as.numeric(x %*% beta_true)
  phi_i <- as.numeric(exp(x %*% alpha_true))
  s_i <- sqrt(phi_i)

  u <- runif(n) # generate uniform (0, 1)

  F <- F_cdf(
    from = from, # cdf
    to = to,
    by = by,
    kappa = kappa_true,
    tau = tau
  )
  z <- F2x(
    u = u, # convert uniform values to values from dist (inverse cdf)
    F = F
  )

  y <- mu_i + (s_i * z)

  return(list(
    x = x,
    y = y,
    beta_true = beta_true,
    alpha_true = alpha_true,
    kappa_true = kappa_true
  ))
}

data_sim_binary <- function(n,
                            beta_true,
                            alpha_true,
                            kappa_true,
                            tau,
                            from = -20,
                            to = 20,
                            by = 0.001,
                            binary_positions) {
  stopifnot(length(beta_true) == length(alpha_true))

  p <- length(beta_true) - 1

  x <- cbind(
    rep(1, n),
    matrix(rnorm(n * p), n, p) # generate X as normal, can also be binary
  )

  # Generate binary columns
  bin_col <- matrix((2 * rbinom(n = n * length(binary_positions), size = 1, prob = 0.5) - 1), n, length(binary_positions)) # transform 2(X)-1 to map from 0, 1 to -1, 1

  x[, binary_positions] <- bin_col # replace x cols with binary columns

  mu_i <- as.numeric(x %*% beta_true)
  phi_i <- as.numeric(exp(x %*% alpha_true))
  s_i <- sqrt(phi_i)

  u <- runif(n) # generate uniform (0, 1)

  F <- F_cdf(
    from = from, # cdf
    to = to,
    by = by,
    kappa = kappa_true,
    tau = tau
  )
  z <- F2x(
    u = u, # convert uniform values to values from dist (inverse cdf)
    F = F
  )

  y <- mu_i + (s_i * z)

  return(list(
    x = x,
    y = y,
    beta_true = beta_true,
    alpha_true = alpha_true,
    kappa_true = kappa_true,
    binary_positions = binary_positions
  ))
}

# =========================================================================#
# Fitting Functions --------------------------------------------------------
# Iterative Loop -----------------------------------------------------------
# For both MPR & SPR
it_loop <- function(theta,
                    FUN,
                    tol,
                    max_it,
                    ...) {
  theta1 <- theta + 1
  it <- 0
  max_step_reached <- NULL
  steps <- NULL
  while (max(abs(theta1 - theta)) > tol & it < max_it) {
    theta1 <- theta
    fit <- FUN(theta1, ...)
    theta <- fit$estimate
    it <- it + 1

    max_step_reached[it] <- fit$steps[1] # store logical 1 if max step count reached
    steps[it] <- fit$steps[2] # store steps taken
  }
  max_iteration_reached <- ifelse(it == max_it, TRUE, FALSE)
  list(
    "estimate" = as.vector(theta), "maximum" = fit$maximum,
    "max_step_reached" = max_step_reached, "steps" = steps,
    "iterations" = c(max_iteration_reached, it)
  )
} # steps & iterations = 1st value 1: true, 0: false for max reached

# =========================================================================#
# ** Smooth GND -----------------------------------------------------------
# *** c_tilde integral at each step ---------------------------------------
density_gnd_approx <- function(x,
                               mu_loc,
                               s_scale,
                               kappa_shape,
                               tau,
                               c_tilde) { # c_tilde = c_ek * c_k
  c_tilde * (1 / s_scale) * exp(-((sqrt((x - mu_loc)^2 + tau^2) - tau) / s_scale)^kappa_shape)
}

density_gnd_approx_standard <- function(x, # mu = 0, s = 1
                                        kappa_shape,
                                        tau,
                                        c_tilde) { # c_tilde = c_ek * c_k
  c_tilde * exp(-((sqrt((x)^2 + tau^2) - tau))^kappa_shape)
}

# ** Get variance of smooth gnd distribution ------------------------------
xj_smooth_gnd <- function(x,
                          j,
                          kappa_shape,
                          tau) {
  c_tilde_val <- c_tilde_integrate_smoothgnd(
    kappa = kappa_shape,
    tau = tau
  )

  density_func_val <- density_gnd_approx_standard(x,
                                                  kappa_shape = kappa_shape,
                                                  tau = tau,
                                                  c_tilde = c_tilde_val
  )
  (x^j) * density_func_val # integrate over this
}

variance_standard_smooth_gnd <- function(kappa_shape,
                                         tau) {
  ex <- integrate(xj_smooth_gnd,
                  lower = -Inf,
                  upper = Inf,
                  j = 1,
                  kappa_shape = kappa_shape,
                  tau = tau
  )$value # should be zero for symmetric distributions
  ex2 <- integrate(xj_smooth_gnd,
                   lower = -Inf,
                   upper = Inf,
                   j = 2,
                   kappa_shape = kappa_shape,
                   tau = tau
  )$value
  # variance exp(x^2) - exp(x)^2
  ex2 - (ex^2)
}



# *** standard integrate function -----------------------------------------
c_tilde_integrate_smoothgnd <- function(kappa,
                                        tau) {
  stopifnot(length(tau) == 1)
  length_kappa <- length(kappa)

  if (length(unique(kappa)) == 1) {
    kappa <- kappa[1]
  }

  output <- NULL
  for (i in kappa) {
    kappa_pos <- which(kappa == i)
    integral <- integrate(density_gnd_approx_standard,
                          lower = -Inf,
                          upper = Inf,
                          # rel.tol = 1e-15,
                          kappa_shape = i,
                          tau = tau,
                          c_tilde = 1 # 1 as unknown here
    )$value
    output[kappa_pos] <- (1 / integral) # reciprocal
  }

  if (length(output) != length_kappa) {
    output <- rep(output, length_kappa)
  } # return length of input kappa

  output # output
}

# *** transform integral --------------------------------------------------
# integrate between 0 and 1
density_gnd_approx_standard_transform <- function(t,
                                                  kappa_shape,
                                                  tau,
                                                  c_tilde) {
  (c_tilde * exp(-((sqrt(((t / (1 - t)))^2 + tau^2) - tau))^kappa_shape)) * (1 / ((1 - t)^2))
}

c_tilde_integrate_transform_smoothgnd <- function(kappa,
                                                  tau) {
  stopifnot(length(tau) == 1)
  length_kappa <- length(kappa)

  if (length(unique(kappa)) == 1) {
    kappa <- kappa[1]
  }

  output <- NULL
  for (i in kappa) {
    kappa_pos <- which(kappa == i)
    integral <- 2 * integrate(density_gnd_approx_standard_transform,
                              lower = 0,
                              upper = 1,
                              kappa_shape = i,
                              tau = tau,
                              c_tilde = 1 # 1 as unknown here
    )$value
    output[kappa_pos] <- (1 / integral) # reciprocal
  }

  if (length(output) != length_kappa) {
    output <- rep(output, length_kappa)
  } # return length of input kappa

  output # output
}

# *** trapezoidal approximation -------------------------------------------
trapezoidal_approximation <- function(FUN,
                                      a,
                                      b,
                                      N,
                                      ...) {
  delta_x <- ((b - a) / N)

  f_x_N <- FUN(
    b,
    ...
  )
  f_x_0 <- FUN(
    a,
    ...
  )
  rhs <- (f_x_N + f_x_0) / 2

  vec <- seq(a, b, length.out = (N + 1))[-1]
  f_x_k <- FUN(
    vec,
    ...
  )

  delta_x * (sum(f_x_k) + rhs)
}

c_tilde_trapezoid_smoothgnd <- function(kappa,
                                        tau,
                                        N = 1000) {
  stopifnot(length(tau) == 1)
  length_kappa <- length(kappa)

  if (length(unique(kappa)) == 1) {
    kappa <- kappa[1]
  }

  output <- NULL
  for (i in kappa) {
    kappa_pos <- which(kappa == i)
    integral <- 2 * trapezoidal_approximation(density_gnd_approx_standard_transform,
                                              a = 0,
                                              b = 1 - (1e-10),
                                              N = N,
                                              kappa_shape = i,
                                              tau = tau,
                                              c_tilde = 1
    ) # 2 * integral between 0 and 1
    output[kappa_pos] <- (1 / integral) # reciprocal
  }

  if (length(output) != length_kappa) {
    output <- rep(output, length_kappa)
  } # return length of input kappa

  output # output
}

# *** rmutil int function -------------------------------------------------
density_gnd_approx_standard_single_input <- function(kappa_shape,
                                                     tau,
                                                     c_tilde) { # c_tilde = c_ek * c_k
  # return function with single input x
  function(x) c_tilde * exp(-((sqrt((x)^2 + tau^2) - tau))^kappa_shape)
}

c_tilde_int_smoothgnd <- function(kappa,
                                  tau) {
  stopifnot(length(tau) == 1)
  length_kappa <- length(kappa)

  if (length(unique(kappa)) == 1) {
    kappa <- kappa[1]
  }

  output <- NULL
  for (i in kappa) {
    kappa_pos <- which(kappa == i)

    dens_func <- density_gnd_approx_standard_single_input(
      kappa_shape = i,
      tau = tau,
      c_tilde = 1
    ) # 1 as unknown here
    integral <- 2 * rmutil::int(dens_func, # only takes functions with a single input
                                a = 0,
                                b = Inf
    )
    output[kappa_pos] <- (1 / integral) # reciprocal
  }

  if (length(output) != length_kappa) {
    output <- rep(output, length_kappa)
  } # return length of input kappa

  output # output
}

# * Finite Difference Method ----------------------------------------------
# *** get values of c_tilde and approximate derivative with finite difference method ----
vals_finite_diff <- function(kappa,
                             tau,
                             h_kappa,
                             FUN_c_tilde, # c_tilde_smoothgnd
                             ...) {
  # f(x) --
  c_tilde_vec <- FUN_c_tilde(
    kappa = kappa,
    tau = tau,
    ...
  )
  # f(x + h) --
  c_tilde_plus_vec <- FUN_c_tilde(
    kappa = (kappa + h_kappa),
    tau = tau,
    ...
  )
  # f(x - h) --
  c_tilde_minus_vec <- FUN_c_tilde(
    kappa = (kappa - h_kappa),
    tau = tau,
    ...
  )
  list(
    "c_tilde_vec" = c_tilde_vec,
    "c_tilde_plus_vec" = c_tilde_plus_vec,
    "c_tilde_minus_vec" = c_tilde_minus_vec
  )
}

# *** finite difference method --------------------------------------------
c_tilde_prime_func_finite_diff <- function(c_tilde_plus_vec,
                                           c_tilde_minus_vec,
                                           h_kappa) {
  (c_tilde_plus_vec - c_tilde_minus_vec) / (2 * h_kappa)
}

c_tilde_double_prime_func_finite_diff <- function(c_tilde_vec,
                                                  c_tilde_plus_vec,
                                                  c_tilde_minus_vec,
                                                  h_kappa) {
  (c_tilde_plus_vec + c_tilde_minus_vec - (2 * c_tilde_vec)) / ((h_kappa)^2)
}

# Using numDeriv package and integral -------------------------------------
c_tilde_prime_func_numDeriv <- function(kappa,
                                        tau,
                                        FUN_c_tilde,
                                        ...) {
  stopifnot(length(tau) == 1)
  length_kappa <- length(kappa)

  if (length(unique(kappa)) == 1) {
    kappa <- kappa[1]
  }

  c_tilde_prime_vec <- numDeriv::grad(FUN_c_tilde,
                                      x = kappa,
                                      tau = tau,
                                      ...
  )

  if (length(c_tilde_prime_vec) != length_kappa) {
    c_tilde_prime_vec <- rep(c_tilde_prime_vec, length_kappa)
  } # return length of input kappa

  c_tilde_prime_vec
}

hessian_vectorized <- Vectorize(numDeriv::hessian, vectorize.args = "x")
c_tilde_double_prime_func_numDeriv <- function(kappa,
                                               tau,
                                               FUN_c_tilde,
                                               ...) {
  stopifnot(length(tau) == 1)
  length_kappa <- length(kappa)

  if (length(unique(kappa)) == 1) {
    kappa <- kappa[1]
  }

  c_tilde_double_prime_vec <- hessian_vectorized(FUN_c_tilde,
                                                 x = kappa,
                                                 tau = tau,
                                                 ...
  )

  if (length(c_tilde_double_prime_vec) != length_kappa) {
    c_tilde_double_prime_vec <- rep(c_tilde_double_prime_vec, length_kappa)
  } # return length of input kappa

  c_tilde_double_prime_vec
}

# Actual derivative of rho function (rho = exp(-abs(y)^k)) ----------------
# rho = exp(-((sqrt((((y)))^2 + tau^2) - tau)^(exp(x_3 * nu) + omega)))
rho_prime <- function(y,
                      kappa,
                      kappa_no_omega,
                      tau) {
  (log(sqrt(y^2 + tau^2) - tau)) * ((-kappa_no_omega) * exp(-(sqrt(y^2 + tau^2) - tau)^(kappa))) * ((sqrt(y^2 + tau^2) - tau)^(kappa))
}

rho_double_prime <- function(y,
                             kappa,
                             kappa_no_omega,
                             tau) {
  ((log(sqrt(y^2 + tau^2) - tau)^2) * (-(kappa_no_omega)^2) * (exp(-(sqrt(y^2 + tau^2) - tau)^(kappa))) * ((sqrt(y^2 + tau^2) - tau)^(kappa))) - ((log(sqrt(y^2 + tau^2) - tau) * (kappa_no_omega) * (exp(-(sqrt(y^2 + tau^2) - tau)^(kappa))) * ((sqrt(y^2 + tau^2) - tau)^(kappa))) * (1 - (kappa_no_omega) * log(sqrt(y^2 + tau^2) - tau) * (sqrt(y^2 + tau^2) - tau)^(kappa)))
}

# get values of integrals of these functions
rho_integrals_smoothgnd <- function(kappa,
                                    kappa_no_omega,
                                    tau) {
  stopifnot(length(tau) == 1)
  length_kappa <- length(kappa)

  if (length(unique(kappa)) == 1) {
    kappa <- kappa[1]
    kappa_no_omega <- kappa_no_omega[1]
  }

  rho_prime_vec <- NULL # initialize vectors
  rho_double_prime_vec <- NULL

  for (i in kappa) {
    kappa_pos <- which(kappa == i)

    integral_prime <- integrate(rho_prime,
                                lower = -Inf,
                                upper = Inf,
                                kappa = i,
                                kappa_no_omega = kappa_no_omega[kappa_pos],
                                tau = tau
    )$value

    integral_double_prime <- integrate(rho_double_prime,
                                       lower = -Inf,
                                       upper = Inf,
                                       kappa = i,
                                       kappa_no_omega = kappa_no_omega[kappa_pos],
                                       tau = tau
    )$value

    rho_prime_vec[kappa_pos] <- integral_prime
    rho_double_prime_vec[kappa_pos] <- integral_double_prime
  }

  stopifnot(length(rho_prime_vec) == length(rho_double_prime_vec))

  if (length(rho_prime_vec) != length_kappa) {
    rho_prime_vec <- rep(rho_prime_vec, length_kappa)
    rho_double_prime_vec <- rep(rho_double_prime_vec, length_kappa)
  } # return length of input kappa

  list(
    "rho_prime_vec" = rho_prime_vec,
    "rho_double_prime_vec" = rho_double_prime_vec
  )
}

# get all integral values - c_tilde, rho and rho prime
vals_actual <- function(kappa,
                        kappa_no_omega,
                        tau,
                        FUN_c_tilde # function to get c_tilde
) {
  # c_tilde --
  c_tilde_vec <- FUN_c_tilde(
    kappa = kappa,
    tau = tau
  )

  # rho and rho_prime --
  rho_list <- rho_integrals_smoothgnd(
    kappa,
    kappa_no_omega,
    tau
  )
  list(
    "c_tilde_vec" = c_tilde_vec,
    "rho_prime_vec" = rho_list$rho_prime_vec,
    "rho_double_prime_vec" = rho_list$rho_double_prime_vec
  )
}

c_tilde_prime_func_actual <- function(c_tilde_vec,
                                      rho_prime_vec) {
  (-rho_prime_vec) / ((c_tilde_vec)^2)
}

c_tilde_double_prime_func_actual <- function(c_tilde_vec,
                                             rho_prime_vec,
                                             rho_double_prime_vec) {
  ((c_tilde_vec * (-rho_double_prime_vec)) + (2 * ((rho_prime_vec)^2))) / ((c_tilde_vec)^3)
}

# *** c_tilde z & W ----------------------------------------------------------
z_c_tilde_smoothgnd <- function(c_tilde_vec,
                                c_tilde_prime_vec) {
  (c_tilde_prime_vec / c_tilde_vec) # c_tilde_prime/c_tilde
}

W_c_tilde_smoothgnd <- function(c_tilde_vec,
                                c_tilde_prime_vec,
                                c_tilde_double_prime_vec) {
  -(((c_tilde_vec * c_tilde_double_prime_vec) - ((c_tilde_prime_vec)^2)) / ((c_tilde_vec)^2)) # negative 2nd deriv
}

# *** function with smooth absolute value ---------------------------------
g_smoothgnd <- function(z,
                        kappa,
                        tau) {
  (sqrt((z)^2 + tau^2) - tau)^kappa
}

# ** Derivatives ----------------------------------------------------------
z_beta_smoothgnd <- function(mu,
                             y,
                             phi,
                             kappa,
                             tau) {
  ((((kappa) * (1 / phi)) * (y - mu) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / sqrt((1 / phi) * ((y - mu)^2) + tau^2))
}

z_alpha_smoothgnd <- function(mu,
                              y,
                              phi,
                              kappa,
                              tau) {
  (((((kappa) * (1 / phi)) * ((y - mu)^2) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / (2 * sqrt((1 / phi) * ((y - mu)^2) + tau^2))) - (1 / 2))
}

z_nu_smoothgnd <- function(mu, # need to include c_tilde
                           y,
                           phi,
                           kappa,
                           kappa_no_omega,
                           tau) {
  ((-kappa_no_omega) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa) * (log(sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)))
}

W_beta_smoothgnd <- function(mu,
                             y,
                             phi,
                             kappa,
                             tau) {
  -((-((kappa - 1) * ((kappa) * ((1 / phi)^2)) * ((y - mu)^2) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 2)) / ((1 / phi) * ((y - mu)^2) + tau^2)) - ((((kappa) * (1 / phi)) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / sqrt((1 / phi) * ((y - mu)^2) + tau^2)) + ((((kappa) * ((1 / phi)^2)) * ((y - mu)^2) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / ((1 / phi) * ((y - mu)^2) + tau^2)^(3 / 2))) # negative second derivative
}

W_alpha_smoothgnd <- function(mu,
                              y,
                              phi,
                              kappa,
                              tau) {
  -((-((kappa - 1) * ((kappa) * ((1 / phi)^2)) * ((y - mu)^4) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 2)) / (4 * ((1 / phi) * ((y - mu)^2) + tau^2))) - ((((kappa) * (1 / phi)) * ((y - mu)^2) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / (2 * sqrt((1 / phi) * ((y - mu)^2) + tau^2))) + ((((kappa) * ((1 / phi)^2)) * ((y - mu)^4) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / (4 * ((1 / phi) * ((y - mu)^2) + tau^2)^(3 / 2)))) # negative second derivative
}

W_nu_smoothgnd <- function(mu, # need to include c_tilde
                           y,
                           phi,
                           kappa,
                           kappa_no_omega,
                           tau) {
  -(((-(kappa_no_omega^2)) * (sqrt((1 / phi) * (y - mu)^2 + tau^2) - tau)^(kappa) * (log(sqrt((1 / phi) * (y - mu)^2 + tau^2) - tau)^2)) - (kappa_no_omega * (sqrt((1 / phi) * (y - mu)^2 + tau^2) - tau)^(kappa) * (log(sqrt((1 / phi) * (y - mu)^2 + tau^2) - tau)))) # negative second derivative
}

W_beta_alpha_smoothgnd <- function(mu,
                                   y,
                                   phi,
                                   kappa,
                                   tau) {
  -((-((kappa - 1) * ((kappa) * ((1 / phi)^2)) * ((y - mu)^3) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 2)) / (2 * ((1 / phi) * ((y - mu)^2) + tau^2))) - ((((kappa) * (1 / phi)) * (y - mu) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / sqrt((1 / phi) * ((y - mu)^2) + tau^2)) + ((((kappa) * ((1 / phi)^2)) * ((y - mu)^3) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / (2 * ((1 / phi) * ((y - mu)^2) + tau^2)^(3 / 2)))) # negative second deriv
}

W_beta_nu_smoothgnd <- function(mu,
                                y,
                                phi,
                                kappa,
                                kappa_no_omega,
                                tau) {
  -(((((kappa_no_omega) * (1 / phi)) * (y - mu) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / sqrt((1 / phi) * ((y - mu)^2) + tau^2)) + ((((kappa_no_omega * kappa) * (1 / phi)) * (y - mu) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1) * log(sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)) / sqrt((1 / phi) * ((y - mu)^2) + tau^2))) # negative second derivative
}

W_alpha_nu_smoothgnd <- function(mu,
                                 y,
                                 phi,
                                 kappa,
                                 kappa_no_omega,
                                 tau) {
  -(((((kappa_no_omega) * (1 / phi)) * ((y - mu)^2) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1)) / (2 * sqrt((1 / phi) * ((y - mu)^2) + tau^2))) + ((((kappa_no_omega * kappa) * (1 / phi)) * ((y - mu)^2) * (sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)^(kappa - 1) * log(sqrt((1 / phi) * ((y - mu)^2) + tau^2) - tau)) / (2 * sqrt((1 / phi) * ((y - mu)^2) + tau^2)))) # negative second derivative
}

# *** list of families ----------------------------------------------------
list_families <- list()
list_families$smoothgnd <- list(
  general = list(
    g = g_smoothgnd, # smooth absolute value
    c_tilde_integrate = c_tilde_integrate_smoothgnd, # use integral to get normalizing constant value
    c_tilde_integrate_transform = c_tilde_integrate_transform_smoothgnd, # standard integrate function but transform integral between 0 & 1
    c_tilde_trapezoid = c_tilde_trapezoid_smoothgnd,
    c_tilde_int = c_tilde_int_smoothgnd # rmutil function
  ),
  derivatives = list(
    z_beta = z_beta_smoothgnd,
    z_alpha = z_alpha_smoothgnd,
    z_nu = z_nu_smoothgnd,
    W_beta = W_beta_smoothgnd,
    W_alpha = W_alpha_smoothgnd,
    W_nu = W_nu_smoothgnd,
    W_beta_alpha = W_beta_alpha_smoothgnd,
    W_beta_nu = W_beta_nu_smoothgnd,
    W_alpha_nu = W_alpha_nu_smoothgnd,
    # z and W for c_tilde
    z_c_tilde = z_c_tilde_smoothgnd,
    W_c_tilde = W_c_tilde_smoothgnd,
    # finite difference functions for approx deriv of c_tilde
    vals_finite_diff = vals_finite_diff,
    c_tilde_prime_func_finite_diff = c_tilde_prime_func_finite_diff,
    c_tilde_double_prime_func_finite_diff = c_tilde_double_prime_func_finite_diff,
    # grad and hessian functions from numDeriv
    c_tilde_prime_func_numDeriv = c_tilde_prime_func_numDeriv,
    c_tilde_double_prime_func_numDeriv = c_tilde_double_prime_func_numDeriv,
    # actual derivatives for c_tilde
    vals_actual = vals_actual,
    c_tilde_prime_func_actual = c_tilde_prime_func_actual,
    c_tilde_double_prime_func_actual = c_tilde_double_prime_func_actual
  )
)

# *** density function ----------------------------------------------------
smooth_gnd_density <- function(x,
                               mu,
                               s_scale, # s parameter
                               kappa,
                               tau,
                               list_general, # list containing general family functions
                               method_c_tilde # "table" or "integral"
) {
  # Get family specific functions
  c_tilde_func <- list_general[[paste0("c_tilde_", method_c_tilde)]] # either table or integral

  c_tilde <- c_tilde_func(
    kappa = kappa,
    tau = tau
  ) # look up normalizing constant value
  g <- list_general$g(
    z = ((x - mu) / s_scale),
    kappa = kappa,
    tau = tau
  ) # smooth absolute value

  # Density function --
  c_tilde * (1 / s_scale) * exp(-g)
}

# *** loglikelihood function ----------------------------------------------
likelihood <- function(theta,
                       x1,
                       x2,
                       x3,
                       y,
                       tau,
                       list_general, # list containing general family functions
                       method_c_tilde, # "table" or "integral"
                       kappa_omega # offset value to restrict kappa
) {
  n <- length(y)

  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1

  beta <- theta[1:(p1 + 1)]
  alpha <- theta[(p1 + 2):(p1 + 2 + p2)]
  nu <- theta[(p1 + 3 + p2):length(theta)]

  mu <- x1 %*% beta
  phi <- exp(x2 %*% alpha) # phi = s^2
  kappa <- nu_to_kappa_matrix(
    nu = nu,
    x3 = x3,
    kappa_omega = kappa_omega
  ) # restrict kappa

  z <- ((y - mu) / sqrt(phi))

  # Get family specific functions --
  c_tilde_func <- list_general[[paste0("c_tilde_", method_c_tilde)]] # table or integral

  c_tilde <- c_tilde_func(
    kappa = kappa,
    tau = tau
  ) # look up normalizing constant value

  g <- list_general$g(
    z = z,
    kappa = kappa,
    tau = tau
  ) # smooth absolute value

  # Likelihood function --
  like <- (sum(log(c_tilde))) - ((0.5) * (sum(x2 %*% alpha))) - sum(g)
  like
}

likelihood_neg <- function(...) -likelihood(...)

likelihood_neg_fixed <- function(theta_estimate, # just coefs to be estimated
                                 fix_indices, # positions of coefs that have been removed already
                                 fix_values, # match fix_indices
                                 x1,
                                 x2,
                                 x3,
                                 y,
                                 tau,
                                 list_general,
                                 method_c_tilde,
                                 kappa_omega) {
  theta <- theta_estimate
  for (i in fix_indices) {
    pos <- which(fix_indices == i)
    theta <- append(
      x = theta,
      values = fix_values[pos],
      after = i - 1
    ) # insert value after fix_index - 1 position
  }

  likelihood_neg(theta = theta,
                 x1 = x1,
                 x2 = x2,
                 x3 = x3,
                 y = y,
                 tau = tau,
                 list_general = list_general,
                 method_c_tilde = method_c_tilde,
                 kappa_omega = kappa_omega)

}

# *** penalized likelihood ------------------------------------------------
likelihood_penalized <- function(theta,
                                 x1,
                                 x2,
                                 x3,
                                 y,
                                 tau,
                                 epsilon,
                                 list_general, # list containing general family functions
                                 method_c_tilde, # "table" or "integral"
                                 kappa_omega, # offset for restricting kappa
                                 lambda_beta,
                                 lambda_alpha,
                                 lambda_nu) {
  # maximize this
  n <- length(y)
  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1

  beta <- theta[1:(p1 + 1)]
  alpha <- theta[(p1 + 2):(p1 + 2 + p2)]
  nu <- theta[(p1 + 3 + p2):length(theta)]

  lambda_beta <- (eval(parse(text = lambda_beta)))
  lambda_alpha <- (eval(parse(text = lambda_alpha)))
  lambda_nu <- (eval(parse(text = lambda_nu)))

  # remove intercepts as won't be penalized
  beta1_p <- beta[-1]
  alpha1_p <- alpha[-1]
  nu1_p <- nu[-1]

  # penalized likelihood
  p_like <- likelihood(
    theta,
    x1,
    x2,
    x3,
    y,
    tau,
    list_general,
    method_c_tilde,
    kappa_omega
  ) - ((lambda_beta / 2) * (1 + (sum((beta1_p^2) / (beta1_p^2 + epsilon^2))))) - ((lambda_alpha / 2) * (1 + (sum((alpha1_p^2) / (alpha1_p^2 + epsilon^2))))) - ((lambda_nu / 2) * (1 + (sum((nu1_p^2) / (nu1_p^2 + epsilon^2)))))
  p_like # divided by two for penalized setup, +1 for estimation of intercepts
}

likelihood_penalized_neg <- function(...) -likelihood_penalized(...)

# Newton Raphson ----------------------------------------------------------
nr <- function(theta,
               fix_index, # = NA if don't want any values fixed
               x1,
               x2,
               x3,
               y,
               tau,
               list_family, # list containing all functions (general and derivatives) for a specific family
               method_c_tilde, # "table" or "integral"
               method_c_tilde_deriv, # "finite_diff" or "actual"
               kappa_omega, # offset to restrict kappa
               initial_step,
               max_step_it,
               algorithm, # RS for cross-derivatives set to zero, CG for true derivatives
               h_kappa, # depends on how fine a grid of kappa
               ...) {
  # Extract list of general & derivative functions --
  list_general <- list_family$general
  list_deriv <- list_family$derivatives

  # Setup --
  n <- length(y)

  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1

  beta <- theta[1:(p1 + 1)]
  alpha <- theta[(p1 + 2):(p1 + 2 + p2)]
  nu <- theta[(p1 + 3 + p2):length(theta)]

  mu <- x1 %*% beta
  phi <- exp(x2 %*% alpha)
  kappa <- nu_to_kappa_matrix(
    nu = nu,
    x3 = x3,
    kappa_omega = kappa_omega
  )
  kappa_no_omega <- kappa - kappa_omega # needed for derivative related to nu

  tx1 <- t(x1) # save transpose of x1 as an object
  tx2 <- t(x2) # save transpose of x2 as an object
  tx3 <- t(x3) # save transpose of x3 as an object

  # Extract derivative functions --
  z_beta <- list_deriv$z_beta(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  z_alpha <- list_deriv$z_alpha(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  z_nu <- list_deriv$z_nu(
    mu,
    y,
    phi,
    kappa,
    kappa_no_omega,
    tau
  )
  W_beta <- list_deriv$W_beta(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  W_alpha <- list_deriv$W_alpha(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  W_nu <- list_deriv$W_nu(
    mu,
    y,
    phi,
    kappa,
    kappa_no_omega,
    tau
  )

  switch(algorithm,
         "CG" = {
           W_beta_alpha <- list_deriv$W_beta_alpha(
             mu,
             y,
             phi,
             kappa,
             tau
           ) # true derivatives

           W_beta_nu <- list_deriv$W_beta_nu(
             mu,
             y,
             phi,
             kappa,
             kappa_no_omega,
             tau
           ) # true derivatives

           W_alpha_nu <- list_deriv$W_alpha_nu(
             mu,
             y,
             phi,
             kappa,
             kappa_no_omega,
             tau
           ) # true derivatives
         },
         "cross_zero" = {
           W_beta_alpha <- rep(0, n) # cross-derivatives = 0 for optimization

           W_beta_nu <- rep(0, n) # cross-derivatives = 0 for optimization

           W_alpha_nu <- rep(0, n) # cross-derivatives = 0 for optimization
         },
         "RS" = {
           W_beta_alpha <- rep(0, n) # cross-derivatives = 0 for optimization

           W_beta_nu <- rep(0, n) # cross-derivatives = 0 for optimization

           W_alpha_nu <- list_deriv$W_alpha_nu(
             mu,
             y,
             phi,
             kappa,
             kappa_no_omega,
             tau
           ) # true derivatives
         }
  )

  # c_tilde --------------------
  c_tilde_func <- list_general[[paste0("c_tilde_", method_c_tilde)]] # table or integral

  # get derivatives of c_tilde ----
  c_tilde_prime_func <- list_deriv[[paste0("c_tilde_prime_func_", method_c_tilde_deriv)]]
  c_tilde_double_prime_func <- list_deriv[[paste0("c_tilde_double_prime_func_", method_c_tilde_deriv)]]

  # derivative of log likelihood for c_tilde
  # get values of c_tilde and derivatives using finite_difference approx
  switch(method_c_tilde_deriv,
         "finite_diff" = {
           vals_func <- list_deriv[[paste0("vals_", method_c_tilde_deriv)]] # either using finite difference or actual derivatives
           c_vals <- vals_func( # look up values in table c_tilde
             kappa = kappa,
             tau = tau,
             h_kappa = h_kappa,
             FUN_c_tilde = c_tilde_func
           ) # c_tilde, c_tilde_plus, c_tilde_minus for finite difference

           c_tilde_vec <- c_vals$c_tilde_vec

           # c_tilde_prime
           c_tilde_prime_vec <- c_tilde_prime_func(
             c_tilde_plus_vec = c_vals$c_tilde_plus_vec,
             c_tilde_minus_vec = c_vals$c_tilde_minus_vec,
             h_kappa = h_kappa
           )

           # c_tilde_double_prime
           c_tilde_double_prime_vec <- c_tilde_double_prime_func(
             c_tilde_vec = c_tilde_vec,
             c_tilde_plus_vec = c_vals$c_tilde_plus_vec,
             c_tilde_minus_vec = c_vals$c_tilde_minus_vec,
             h_kappa = h_kappa
           )
         },
         "numDeriv" = { # get values using grad and hessian functions from numDeriv
           c_tilde_vec <- c_tilde_func(
             kappa = kappa,
             tau = tau
           )

           c_tilde_prime_vec <- c_tilde_prime_func(
             kappa = kappa,
             tau = tau,
             FUN_c_tilde = c_tilde_func
           )
           c_tilde_double_prime_vec <- c_tilde_double_prime_func(
             kappa = kappa,
             tau = tau,
             FUN_c_tilde = c_tilde_func
           )
         },
         "actual" = { # get values of c_tilde and derivatives using actual
           vals_func <- list_deriv[[paste0("vals_", method_c_tilde_deriv)]] # either using finite difference or actual derivatives
           c_vals <- vals_func(
             kappa = kappa,
             kappa_no_omega = kappa_no_omega,
             tau = tau,
             FUN_c_tilde = c_tilde_func
           ) # return c_tilde_vec, rho_prime_vec, rho_double_prime_vec

           c_tilde_vec <- c_vals$c_tilde_vec

           c_tilde_prime_vec <- c_tilde_prime_func(
             c_tilde_vec = c_tilde_vec,
             rho_prime_vec = c_vals$rho_prime_vec
           )

           c_tilde_double_prime_vec <- c_tilde_double_prime_func(
             c_tilde_vec = c_tilde_vec,
             rho_prime_vec = c_vals$rho_prime_vec,
             rho_double_prime_vec = c_vals$rho_double_prime_vec
           )
         }
  )

  # Get z and W ----
  # z_c_tilde
  z_c_tilde <- list_deriv$z_c_tilde(
    c_tilde_vec = c_tilde_vec,
    c_tilde_prime_vec = c_tilde_prime_vec
  )

  # W_c_tilde
  W_c_tilde <- list_deriv$W_c_tilde(
    c_tilde_vec = c_tilde_vec,
    c_tilde_prime_vec = c_tilde_prime_vec,
    c_tilde_double_prime_vec = c_tilde_double_prime_vec
  )

  ## Score --
  ## Beta
  score_beta <- tx1 %*% z_beta

  ## Alpha
  score_alpha <- tx2 %*% z_alpha

  ## Nu
  score_nu <- tx3 %*% (z_nu + z_c_tilde) # addition of derivative of c_tilde

  score <- c(score_beta, score_alpha, score_nu)

  ## Information Matrix --
  ## Beta
  I_beta <- (t(x1 * c(W_beta))) %*% x1

  ## Alpha
  I_alpha <- (t(x2 * c(W_alpha))) %*% x2

  ## Nu
  I_nu <- (t(x3 * c(W_nu + W_c_tilde))) %*% x3 # addition of derivative of c_tilde

  ## Beta, Alpha, Nu
  I_beta_alpha <- t(x1 * c(W_beta_alpha)) %*% x2
  I_alpha_beta <- t(I_beta_alpha)

  I_beta_nu <- t(x1 * c(W_beta_nu)) %*% x3
  I_nu_beta <- t(I_beta_nu)

  I_alpha_nu <- t(x2 * c(W_alpha_nu)) %*% x3
  I_nu_alpha <- t(I_alpha_nu)

  information <- rbind(
    cbind(I_beta, I_beta_alpha, I_beta_nu),
    cbind(I_alpha_beta, I_alpha, I_alpha_nu),
    cbind(I_nu_beta, I_nu_alpha, I_nu)
  )

  # Remove fix_index from score and information matrix ---
  if (!any(is.na(fix_index))) { # if given fix_index
    score <- score[-fix_index]
    information <- information[-fix_index, -fix_index] # remove row and col corresponding to fix index

    delta <- solve(information, score)

    # include 0 where in fix_index position so delta vector correct length
    for (i in fix_index) {
      delta <- append(
        x = delta,
        values = 0,
        after = i - 1
      ) # insert 0 after fix_index - 1 position
    }
  } else if (any(is.na(fix_index))) { # if fix_index == NA, then get delta as normal
    delta <- solve(information, score)
  }

  ## Step halving
  like_current <- likelihood(
    theta,
    x1,
    x2,
    x3,
    y,
    tau,
    list_general,
    method_c_tilde,
    kappa_omega
  )
  like_new <- like_current - 1
  j <- 0
  while (like_new < like_current & j < max_step_it) {
    theta_new <- theta + ((initial_step * delta) / (2^j))
    like_new <- likelihood(
      theta_new,
      x1,
      x2,
      x3,
      y,
      tau,
      list_general,
      method_c_tilde,
      kappa_omega
    )
    j <- j + 1
  }
  max_step_reached <- ifelse(j == max_step_it, TRUE, FALSE)
  return(list(
    "estimate" = theta_new, "maximum" = like_new,
    "steps" = c(max_step_reached, j),
    "fix_index" = fix_index
  ))
}

# Penalized Newton-Raphson ------------------------------------------------
nr_penalty <- function(theta,
                       fix_index, # = NA if don't want it fixed
                       x1,
                       x2,
                       x3,
                       y,
                       tau,
                       epsilon,
                       list_family, # list containing all functions (general and derivatives) for a specific family
                       method_c_tilde, # "table" or "integral"
                       method_c_tilde_deriv, # "finite_diff" or "actual"
                       kappa_omega, # offset to restrict kappa
                       lambda_beta,
                       lambda_alpha,
                       lambda_nu,
                       initial_step,
                       max_step_it,
                       algorithm,
                       h_kappa, # depends on how fine a grid of kappa
                       ...) {
  # Extract list of general & derivative functions --
  list_general <- list_family$general
  list_deriv <- list_family$derivatives

  # Setup --
  n <- length(y)

  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1
  p3 <- ncol(x3) - 1

  beta <- theta[1:(p1 + 1)]
  alpha <- theta[(p1 + 2):(p1 + 2 + p2)]
  nu <- theta[(p1 + 3 + p2):length(theta)]

  beta1_p <- beta[-1] # remove intercepts as won't be penalized
  alpha1_p <- alpha[-1]
  nu1_p <- nu[-1]

  mu <- x1 %*% beta
  phi <- exp(x2 %*% alpha)
  kappa <- nu_to_kappa_matrix(
    nu = nu,
    x3 = x3,
    kappa_omega = kappa_omega
  )
  kappa_no_omega <- kappa - kappa_omega # needed for derivatives related to nu

  tx1 <- t(x1) # save transpose of x1 as an object
  tx2 <- t(x2) # save transpose of x2 as an object
  tx3 <- t(x3) # save transpose of x3 as an object

  lambda_beta <- (eval(parse(text = lambda_beta)))
  lambda_alpha <- (eval(parse(text = lambda_alpha)))
  lambda_nu <- (eval(parse(text = lambda_nu)))

  # Extract derivative functions --
  z_beta <- list_deriv$z_beta(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  z_alpha <- list_deriv$z_alpha(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  z_nu <- list_deriv$z_nu(
    mu,
    y,
    phi,
    kappa,
    kappa_no_omega,
    tau
  )
  W_beta <- list_deriv$W_beta(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  W_alpha <- list_deriv$W_alpha(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  W_nu <- list_deriv$W_nu(
    mu,
    y,
    phi,
    kappa,
    kappa_no_omega,
    tau
  )

  switch(algorithm,
         "CG" = {
           W_beta_alpha <- list_deriv$W_beta_alpha(
             mu,
             y,
             phi,
             kappa,
             tau
           ) # true derivatives

           W_beta_nu <- list_deriv$W_beta_nu(
             mu,
             y,
             phi,
             kappa,
             kappa_no_omega,
             tau
           ) # true derivatives

           W_alpha_nu <- list_deriv$W_alpha_nu(
             mu,
             y,
             phi,
             kappa,
             kappa_no_omega,
             tau
           ) # true derivatives
         },
         "cross_zero" = {
           W_beta_alpha <- rep(0, n) # cross-derivatives = 0 for optimization

           W_beta_nu <- rep(0, n) # cross-derivatives = 0 for optimization

           W_alpha_nu <- rep(0, n) # cross-derivatives = 0 for optimization
         },
         "RS" = {
           W_beta_alpha <- rep(0, n) # cross-derivatives = 0 for optimization

           W_beta_nu <- rep(0, n) # cross-derivatives = 0 for optimization

           W_alpha_nu <- list_deriv$W_alpha_nu(
             mu,
             y,
             phi,
             kappa,
             kappa_no_omega,
             tau
           ) # true derivatives
         }
  )

  # c_tilde --------------------
  # derivative of log likelihood for c_tilde
  c_tilde_func <- list_general[[paste0("c_tilde_", method_c_tilde)]] # table or integral

  # get derivatives of c_tilde ----
  c_tilde_prime_func <- list_deriv[[paste0("c_tilde_prime_func_", method_c_tilde_deriv)]]
  c_tilde_double_prime_func <- list_deriv[[paste0("c_tilde_double_prime_func_", method_c_tilde_deriv)]]

  # derivative of log likelihood for c_tilde
  # get values of c_tilde and approximate derivatives using finite diff
  switch(method_c_tilde_deriv,
         "finite_diff" = {
           vals_func <- list_deriv[[paste0("vals_", method_c_tilde_deriv)]] # either using finite difference or actual derivatives
           c_vals <- vals_func( # look up values in table c_tilde
             kappa = kappa,
             tau = tau,
             h_kappa = h_kappa,
             FUN_c_tilde = c_tilde_func
           ) # c_tilde, c_tilde_plus, c_tilde_minus for finite difference

           c_tilde_vec <- c_vals$c_tilde_vec

           # c_tilde_prime
           c_tilde_prime_vec <- c_tilde_prime_func(
             c_tilde_plus_vec = c_vals$c_tilde_plus_vec,
             c_tilde_minus_vec = c_vals$c_tilde_minus_vec,
             h_kappa = h_kappa
           )

           # c_tilde_double_prime
           c_tilde_double_prime_vec <- c_tilde_double_prime_func(
             c_tilde_vec = c_tilde_vec,
             c_tilde_plus_vec = c_vals$c_tilde_plus_vec,
             c_tilde_minus_vec = c_vals$c_tilde_minus_vec,
             h_kappa = h_kappa
           )
         },
         "numDeriv" = { # get values using grad and hessian functions from numDeriv
           c_tilde_vec <- c_tilde_func(
             kappa = kappa,
             tau = tau
           )

           c_tilde_prime_vec <- c_tilde_prime_func(
             kappa = kappa,
             tau = tau,
             FUN_c_tilde = c_tilde_func
           )
           c_tilde_double_prime_vec <- c_tilde_double_prime_func(
             kappa = kappa,
             tau = tau,
             FUN_c_tilde = c_tilde_func
           )
         },
         "actual" = { # get values of c_tilde and derivatives using actual
           vals_func <- list_deriv[[paste0("vals_", method_c_tilde_deriv)]] # either using finite difference or actual derivatives
           c_vals <- vals_func(
             kappa = kappa,
             kappa_no_omega = kappa_no_omega,
             tau = tau,
             FUN_c_tilde = c_tilde_func
           ) # return c_tilde_vec, rho_prime_vec, rho_double_prime_vec

           c_tilde_vec <- c_vals$c_tilde_vec

           c_tilde_prime_vec <- c_tilde_prime_func(
             c_tilde_vec = c_tilde_vec,
             rho_prime_vec = c_vals$rho_prime_vec
           )

           c_tilde_double_prime_vec <- c_tilde_double_prime_func(
             c_tilde_vec = c_tilde_vec,
             rho_prime_vec = c_vals$rho_prime_vec,
             rho_double_prime_vec = c_vals$rho_double_prime_vec
           )
         }
  )

  # Get z and W ----
  # z_c_tilde
  z_c_tilde <- list_deriv$z_c_tilde(
    c_tilde_vec = c_tilde_vec,
    c_tilde_prime_vec = c_tilde_prime_vec
  )

  # W_c_tilde
  W_c_tilde <- list_deriv$W_c_tilde(
    c_tilde_vec = c_tilde_vec,
    c_tilde_prime_vec = c_tilde_prime_vec,
    c_tilde_double_prime_vec = c_tilde_double_prime_vec
  )

  ## Score --
  ## Beta
  v_beta <- c(0, (lambda_beta / 2) * ((2 * beta1_p * (epsilon^2)) / ((beta1_p^2 + epsilon^2)^2)))
  score_beta <- (tx1 %*% z_beta) - v_beta

  ## Alpha
  v_alpha <- c(0, (lambda_alpha / 2) * ((2 * alpha1_p * (epsilon^2)) / ((alpha1_p^2 + epsilon^2)^2)))
  score_alpha <- (tx2 %*% z_alpha) - v_alpha

  ## Nu
  v_nu <- c(0, (lambda_nu / 2) * ((2 * nu1_p * (epsilon^2)) / ((nu1_p^2 + epsilon^2)^2)))
  score_nu <- (tx3 %*% (z_nu + z_c_tilde)) - v_nu # addition of derivative of c_tilde

  score <- c(score_beta, score_alpha, score_nu)

  ## Information Matrix --
  ## Beta
  Sigma_beta <- diag(c(
    0,
    (lambda_beta / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
  ),
  nrow = p1 + 1,
  ncol = p1 + 1
  )
  I_beta <- ((t(x1 * c(W_beta))) %*% x1) + Sigma_beta

  ## Alpha
  Sigma_alpha <- diag(c(
    0,
    (lambda_alpha / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (alpha1_p^2))) / ((alpha1_p^2 + epsilon^2)^3))
  ),
  nrow = p2 + 1,
  ncol = p2 + 1
  )
  I_alpha <- ((t(x2 * c(W_alpha))) %*% x2) + Sigma_alpha

  ## Nu
  Sigma_nu <- diag(c(
    0,
    (lambda_nu / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (nu1_p^2))) / ((nu1_p^2 + epsilon^2)^3))
  ),
  nrow = p3 + 1,
  ncol = p3 + 1
  )
  I_nu <- ((t(x3 * c(W_nu + W_c_tilde))) %*% x3) + Sigma_nu # addition of derivative of c_tilde

  ## Beta, Alpha, Nu
  I_beta_alpha <- t(x1 * c(W_beta_alpha)) %*% x2
  I_alpha_beta <- t(I_beta_alpha)

  I_beta_nu <- t(x1 * c(W_beta_nu)) %*% x3
  I_nu_beta <- t(I_beta_nu)

  I_alpha_nu <- t(x2 * c(W_alpha_nu)) %*% x3
  I_nu_alpha <- t(I_alpha_nu)

  information <- rbind(
    cbind(I_beta, I_beta_alpha, I_beta_nu),
    cbind(I_alpha_beta, I_alpha, I_alpha_nu),
    cbind(I_nu_beta, I_nu_alpha, I_nu)
  )

  # Remove fix_index from score and information matrix ---
  if (!any(is.na(fix_index))) { # if given fix_index
    score <- score[-fix_index]
    information <- information[-fix_index, -fix_index] # remove row and col corresponding to fix index

    delta <- solve(information, score)

    # include 0 where in fix_index position so delta vector correct length
    for (i in fix_index) {
      delta <- append(
        x = delta,
        values = 0,
        after = i - 1
      ) # insert 0 after fix_index - 1 position
    }
  } else if (any(is.na(fix_index))) { # if fix_index == NA, then get delta as normal
    delta <- solve(information, score)
  }

  ## Step halving
  like_current <- likelihood_penalized(
    theta,
    x1,
    x2,
    x3,
    y,
    tau,
    epsilon,
    list_general,
    method_c_tilde,
    kappa_omega,
    lambda_beta,
    lambda_alpha,
    lambda_nu
  )
  like_new <- like_current - 1
  j <- 0
  while (like_new < like_current & j < max_step_it) {
    theta_new <- theta + ((initial_step * delta) / (2^j))
    like_new <- likelihood_penalized( # this is the BIC type likelihood to be maximized
      theta_new,
      x1,
      x2,
      x3,
      y,
      tau,
      epsilon,
      list_general,
      method_c_tilde,
      kappa_omega,
      lambda_beta,
      lambda_alpha,
      lambda_nu
    )
    j <- j + 1
  }
  max_step_reached <- ifelse(j == max_step_it, TRUE, FALSE)
  return(list(
    "estimate" = theta_new, "maximum" = like_new,
    "steps" = c(max_step_reached, j),
    "fix_index" = fix_index
  ))
}


# Information Matrices ----------------------------------------------------
information_matrices_fullopt <- function(theta,
                                         x1,
                                         x2,
                                         x3,
                                         y,
                                         tau,
                                         epsilon,
                                         list_family,
                                         method_c_tilde,
                                         method_c_tilde_deriv,
                                         h_kappa,
                                         kappa_omega,
                                         lambda_beta,
                                         lambda_alpha,
                                         lambda_nu) {
  # Extract list of general & derivative functions --
  list_general <- list_family$general
  list_deriv <- list_family$derivatives

  # Setup --
  n <- length(y)

  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1
  p3 <- ncol(x3) - 1

  beta <- theta[1:(p1 + 1)]
  alpha <- theta[(p1 + 2):(p1 + 2 + p2)]
  nu <- theta[(p1 + 3 + p2):length(theta)]

  beta1_p <- beta[-1] # remove intercepts as won't be penalized
  alpha1_p <- alpha[-1]
  nu1_p <- nu[-1]

  mu <- x1 %*% beta
  phi <- exp(x2 %*% alpha)
  kappa <- nu_to_kappa_matrix(
    nu = nu,
    x3 = x3,
    kappa_omega = kappa_omega
  )
  kappa_no_omega <- kappa - kappa_omega # needed for derivative related to nu

  tx1 <- t(x1) # save transpose of x1 as an object
  tx2 <- t(x2) # save transpose of x2 as an object
  tx3 <- t(x3) # save transpose of x3 as an object

  lambda_beta <- (eval(parse(text = lambda_beta)))
  lambda_alpha <- (eval(parse(text = lambda_alpha)))
  lambda_nu <- (eval(parse(text = lambda_nu)))

  W_beta <- list_deriv$W_beta(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  W_alpha <- list_deriv$W_alpha(
    mu,
    y,
    phi,
    kappa,
    tau
  )
  W_nu <- list_deriv$W_nu(
    mu,
    y,
    phi,
    kappa,
    kappa_no_omega,
    tau
  )

  # True Cross Derivatives
  W_beta_alpha <- list_deriv$W_beta_alpha(
    mu,
    y,
    phi,
    kappa,
    tau
  ) # true derivatives

  W_beta_nu <- list_deriv$W_beta_nu(
    mu,
    y,
    phi,
    kappa,
    kappa_no_omega,
    tau
  ) # true derivatives

  W_alpha_nu <- list_deriv$W_alpha_nu(
    mu,
    y,
    phi,
    kappa,
    kappa_no_omega,
    tau
  ) # true derivatives

  # c_tilde --------------------
  c_tilde_func <- list_general[[paste0("c_tilde_", method_c_tilde)]] # table or integral

  # get derivatives of c_tilde ----
  c_tilde_prime_func <- list_deriv[[paste0("c_tilde_prime_func_", method_c_tilde_deriv)]]
  c_tilde_double_prime_func <- list_deriv[[paste0("c_tilde_double_prime_func_", method_c_tilde_deriv)]]

  # derivative of log likelihood for c_tilde
  # get values of c_tilde and derivatives using finite_difference approx
  switch(method_c_tilde_deriv,
         "finite_diff" = {
           vals_func <- list_deriv[[paste0("vals_", method_c_tilde_deriv)]] # either using finite difference or actual derivatives
           c_vals <- vals_func( # look up values in table c_tilde
             kappa = kappa,
             tau = tau,
             h_kappa = h_kappa,
             FUN_c_tilde = c_tilde_func
           ) # c_tilde, c_tilde_plus, c_tilde_minus for finite difference

           c_tilde_vec <- c_vals$c_tilde_vec

           # c_tilde_prime
           c_tilde_prime_vec <- c_tilde_prime_func(
             c_tilde_plus_vec = c_vals$c_tilde_plus_vec,
             c_tilde_minus_vec = c_vals$c_tilde_minus_vec,
             h_kappa = h_kappa
           )

           # c_tilde_double_prime
           c_tilde_double_prime_vec <- c_tilde_double_prime_func(
             c_tilde_vec = c_tilde_vec,
             c_tilde_plus_vec = c_vals$c_tilde_plus_vec,
             c_tilde_minus_vec = c_vals$c_tilde_minus_vec,
             h_kappa = h_kappa
           )
         },
         "numDeriv" = { # get values using grad and hessian functions from numDeriv
           c_tilde_vec <- c_tilde_func(
             kappa = kappa,
             tau = tau
           )

           c_tilde_prime_vec <- c_tilde_prime_func(
             kappa = kappa,
             tau = tau,
             FUN_c_tilde = c_tilde_func
           )
           c_tilde_double_prime_vec <- c_tilde_double_prime_func(
             kappa = kappa,
             tau = tau,
             FUN_c_tilde = c_tilde_func
           )
         },
         "actual" = { # get values of c_tilde and derivatives using actual
           vals_func <- list_deriv[[paste0("vals_", method_c_tilde_deriv)]] # either using finite difference or actual derivatives
           c_vals <- vals_func(
             kappa = kappa,
             kappa_no_omega = kappa_no_omega,
             tau = tau,
             FUN_c_tilde = c_tilde_func
           ) # return c_tilde_vec, rho_prime_vec, rho_double_prime_vec

           c_tilde_vec <- c_vals$c_tilde_vec

           c_tilde_prime_vec <- c_tilde_prime_func(
             c_tilde_vec = c_tilde_vec,
             rho_prime_vec = c_vals$rho_prime_vec
           )

           c_tilde_double_prime_vec <- c_tilde_double_prime_func(
             c_tilde_vec = c_tilde_vec,
             rho_prime_vec = c_vals$rho_prime_vec,
             rho_double_prime_vec = c_vals$rho_double_prime_vec
           )
         }
  )

  # W_c_tilde
  W_c_tilde <- list_deriv$W_c_tilde(
    c_tilde_vec = c_tilde_vec,
    c_tilde_prime_vec = c_tilde_prime_vec,
    c_tilde_double_prime_vec = c_tilde_double_prime_vec
  )


  # Penalized Observed Information Matrix ----
  ## Beta
  Sigma_beta <- diag(c(
    0,
    (lambda_beta / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
  ),
  nrow = p1 + 1,
  ncol = p1 + 1
  )
  I_beta <- ((t(x1 * c(W_beta))) %*% x1) + Sigma_beta

  ## Alpha
  Sigma_alpha <- diag(c(
    0,
    (lambda_alpha / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (alpha1_p^2))) / ((alpha1_p^2 + epsilon^2)^3))
  ),
  nrow = p2 + 1,
  ncol = p2 + 1
  )
  I_alpha <- ((t(x2 * c(W_alpha))) %*% x2) + Sigma_alpha

  ## Nu
  Sigma_nu <- diag(c(
    0,
    (lambda_nu / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (nu1_p^2))) / ((nu1_p^2 + epsilon^2)^3))
  ),
  nrow = p3 + 1,
  ncol = p3 + 1
  )
  I_nu <- ((t(x3 * c(W_nu + W_c_tilde))) %*% x3) + Sigma_nu # addition of derivative of c_tilde

  ## Beta, Alpha, Nu
  I_beta_alpha <- t(x1 * c(W_beta_alpha)) %*% x2
  I_alpha_beta <- t(I_beta_alpha)

  I_beta_nu <- t(x1 * c(W_beta_nu)) %*% x3
  I_nu_beta <- t(I_beta_nu)

  I_alpha_nu <- t(x2 * c(W_alpha_nu)) %*% x3
  I_nu_alpha <- t(I_alpha_nu)

  observed_information_penalized <- rbind(
    cbind(I_beta, I_beta_alpha, I_beta_nu),
    cbind(I_alpha_beta, I_alpha, I_alpha_nu),
    cbind(I_nu_beta, I_nu_alpha, I_nu)
  )

  # Unpenalized Information evaluated at penalized estimates
  ## Beta
  I_beta_unpen <- (t(x1 * c(W_beta))) %*% x1

  ## Alpha
  I_alpha_unpen <- (t(x2 * c(W_alpha))) %*% x2

  ## Nu
  I_nu_unpen <- (t(x3 * c(W_nu + W_c_tilde))) %*% x3 # addition of derivative of c_tilde

  observed_information_unpenalized <- rbind(
    cbind(I_beta_unpen, I_beta_alpha, I_beta_nu),
    cbind(I_alpha_beta, I_alpha_unpen, I_alpha_nu),
    cbind(I_nu_beta, I_nu_alpha, I_nu_unpen)
  )

  list(
    "observed_information_penalized" = observed_information_penalized,
    "observed_information_unpenalized" = observed_information_unpenalized
  )
}

information_matrices_fullopt_remove <- function(theta,
                                                x1,
                                                x2,
                                                x3,
                                                y,
                                                tau,
                                                epsilon,
                                                list_family,
                                                method_c_tilde,
                                                method_c_tilde_deriv,
                                                h_kappa,
                                                kappa_omega,
                                                lambda_beta,
                                                lambda_alpha,
                                                lambda_nu,
                                                less_than_value) {
  info_list <- information_matrices_fullopt(
    theta = theta,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    tau = tau,
    epsilon = epsilon,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    h_kappa = h_kappa,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu
  )
  # p <- ncol(x1) - 1
  # int_pos <- c(1, p + 2, (2 * p) + 3)
  # index_keep <- which(!(abs(theta) < less_than_value))

  p <- ncol(x1) - 1
  int_pos <- c(1, p + 2)
  index_keep <- which(!abs(theta) < less_than_value)
  index_keep <- index_keep[which(index_keep != ((2 * p) + 3))] # remove kappa pos

  int_incl <- int_pos[(which(!(int_pos %in% index_keep)))]

  index_keep <- sort(c(index_keep, int_incl)) # include intercept positions if required

  list(
    "observed_information_penalized" = info_list$"observed_information_penalized"[index_keep, index_keep],
    "observed_information_unpenalized" = info_list$"observed_information_unpenalized"[index_keep, index_keep],
    "index_keep" = index_keep
  )
}

information_matrices_numDeriv <- function(theta,
                                          x1,
                                          x2,
                                          x3,
                                          y,
                                          tau,
                                          epsilon,
                                          list_family,
                                          method_c_tilde,
                                          kappa_omega,
                                          lambda_beta,
                                          lambda_alpha,
                                          lambda_nu) {
  list_general <- list_family$general

  observed_information_penalized <- tryCatch(-numDeriv::hessian(likelihood_penalized, # negative second deriv
                                                                theta,
                                                                x1 = x1,
                                                                x2 = x2,
                                                                x3 = x3,
                                                                y = y,
                                                                tau = tau,
                                                                epsilon = epsilon,
                                                                list_general = list_general,
                                                                method_c_tilde = method_c_tilde,
                                                                kappa_omega = kappa_omega,
                                                                lambda_beta = lambda_beta,
                                                                lambda_alpha = lambda_alpha,
                                                                lambda_nu = lambda_nu
  ),
  error = function(err) NA
  )

  observed_information_unpenalized <- tryCatch(-numDeriv::hessian(likelihood, # negative second deriv, unpen
                                                                  theta,
                                                                  x1 = x1,
                                                                  x2 = x2,
                                                                  x3 = x3,
                                                                  y = y,
                                                                  tau = tau,
                                                                  list_general = list_general,
                                                                  method_c_tilde = method_c_tilde,
                                                                  kappa_omega = kappa_omega
  ),
  error = function(err) NA
  )

  list(
    "observed_information_penalized" = observed_information_penalized,
    "observed_information_unpenalized" = observed_information_unpenalized
  )
}

information_matrices_numDeriv_remove <- function(theta,
                                                 x1,
                                                 x2,
                                                 x3,
                                                 y,
                                                 tau,
                                                 epsilon,
                                                 list_family,
                                                 method_c_tilde,
                                                 kappa_omega,
                                                 lambda_beta,
                                                 lambda_alpha,
                                                 lambda_nu,
                                                 less_than_value) {
  info_list <- information_matrices_numDeriv(
    theta = theta,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    tau = tau,
    epsilon = epsilon,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu
  )
  p <- ncol(x1) - 1
  int_pos <- c(1, p + 2, (2 * p) + 3)
  index_keep <- which(!(abs(theta) < less_than_value))

  # int_pos <- c(1, p + 2)
  # index_keep <- which(!abs(theta) < less_than_value)
  # index_keep <- index_keep[which(index_keep != ((2 * p) + 3))] # remove kappa pos

  int_incl <- int_pos[(which(!(int_pos %in% index_keep)))]

  index_keep <- sort(c(index_keep, int_incl)) # include intercept positions if required

  list(
    "observed_information_penalized" = info_list$"observed_information_penalized"[index_keep, index_keep],
    "observed_information_unpenalized" = info_list$"observed_information_unpenalized"[index_keep, index_keep],
    "index_keep" = index_keep
  )
}

information_matrices_numDeriv_parallel_single <- function(list_input, # will contain theta_est, x, y
                                                          tau,
                                                          epsilon,
                                                          list_family,
                                                          method_c_tilde,
                                                          kappa_omega,
                                                          lambda_beta,
                                                          lambda_alpha,
                                                          lambda_nu) {
  n_now <- nrow(list_input$x)

  information_matrices_numDeriv(
    theta = list_input$theta_est,
    x1 = list_input$x,
    x2 = list_input$x,
    x3 = as.matrix(rep(1, n_now)),
    y = list_input$y,
    tau = tau,
    epsilon = epsilon,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu
  )
}

information_matrices_numDeriv_parallel_multiple <- function(list_input_all, # will contain lists of theta_est, x, y
                                                            tau,
                                                            epsilon,
                                                            list_family,
                                                            method_c_tilde,
                                                            kappa_omega,
                                                            lambda_beta,
                                                            lambda_alpha,
                                                            lambda_nu) {
  parLapply(
    cl,
    list_input_all,
    information_matrices_numDeriv_parallel_single,
    tau = tau,
    epsilon = epsilon,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu
  )
}

information_matrices_numDeriv_remove_parallel_single <- function(list_input, # will contain theta_est, x, y
                                                                 tau,
                                                                 epsilon,
                                                                 list_family,
                                                                 method_c_tilde,
                                                                 kappa_omega,
                                                                 lambda_beta,
                                                                 lambda_alpha,
                                                                 lambda_nu,
                                                                 less_than_value) {
  n_now <- nrow(list_input$x)

  information_matrices_numDeriv_remove(
    theta = list_input$theta_est,
    x1 = list_input$x,
    x2 = list_input$x,
    x3 = as.matrix(rep(1, n_now)),
    y = list_input$y,
    tau = tau,
    epsilon = epsilon,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    less_than_value = less_than_value
  )
}

information_matrices_numDeriv_remove_parallel_multiple <- function(list_input_all, # will contain lists of theta_est, x, y
                                                                   tau,
                                                                   epsilon,
                                                                   list_family,
                                                                   method_c_tilde,
                                                                   kappa_omega,
                                                                   lambda_beta,
                                                                   lambda_alpha,
                                                                   lambda_nu,
                                                                   less_than_value) {
  parLapply(
    cl,
    list_input_all,
    information_matrices_numDeriv_remove_parallel_single,
    tau = tau,
    epsilon = epsilon,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    less_than_value = less_than_value
  )
}


information_matrices_choice <- function(optimizer, # change between analytical and nlm
                                        theta,
                                        x1,
                                        x2,
                                        x3,
                                        y,
                                        tau,
                                        epsilon,
                                        list_family,
                                        method_c_tilde,
                                        method_c_tilde_deriv,
                                        h_kappa,
                                        kappa_omega,
                                        lambda_beta,
                                        lambda_alpha,
                                        lambda_nu) {
  switch(optimizer,
         "fullopt" = {
           info_list <- information_matrices_fullopt(
             theta = theta,
             x1 = x1,
             x2 = x2,
             x3 = x3,
             y = y,
             tau = tau,
             epsilon = epsilon,
             list_family = list_family,
             method_c_tilde = method_c_tilde,
             method_c_tilde_deriv = method_c_tilde_deriv,
             h_kappa = h_kappa,
             kappa_omega = kappa_omega,
             lambda_beta = lambda_beta,
             lambda_alpha = lambda_alpha,
             lambda_nu = lambda_nu
           )
         },
         "nlm" = {
           info_list <- information_matrices_numDeriv(
             theta = theta,
             x1 = x1,
             x2 = x2,
             x3 = x3,
             y = y,
             tau = tau,
             epsilon = epsilon,
             list_family = list_family,
             method_c_tilde = method_c_tilde,
             kappa_omega = kappa_omega,
             lambda_beta = lambda_beta,
             lambda_alpha = lambda_alpha,
             lambda_nu = lambda_nu
           )
         }
  )
  info_list
}

information_matrices_choice_remove <- function(optimizer, # change between analytical and nlm
                                               theta,
                                               x1,
                                               x2,
                                               x3,
                                               y,
                                               tau,
                                               epsilon,
                                               list_family,
                                               method_c_tilde,
                                               method_c_tilde_deriv,
                                               h_kappa,
                                               kappa_omega,
                                               lambda_beta,
                                               lambda_alpha,
                                               lambda_nu,
                                               less_than_value) {
  switch(optimizer,
         "fullopt" = {
           info_list <- information_matrices_fullopt_remove(
             theta = theta,
             x1 = x1,
             x2 = x2,
             x3 = x3,
             y = y,
             tau = tau,
             epsilon = epsilon,
             list_family = list_family,
             method_c_tilde = method_c_tilde,
             method_c_tilde_deriv = method_c_tilde_deriv,
             h_kappa = h_kappa,
             kappa_omega = kappa_omega,
             lambda_beta = lambda_beta,
             lambda_alpha = lambda_alpha,
             lambda_nu = lambda_nu,
             less_than_value = less_than_value
           )
         },
         "nlm" = {
           info_list <- information_matrices_numDeriv_remove(
             theta = theta,
             x1 = x1,
             x2 = x2,
             x3 = x3,
             y = y,
             tau = tau,
             epsilon = epsilon,
             list_family = list_family,
             method_c_tilde = method_c_tilde,
             kappa_omega = kappa_omega,
             lambda_beta = lambda_beta,
             lambda_alpha = lambda_alpha,
             lambda_nu = lambda_nu,
             less_than_value = less_than_value
           )
         }
  )
  info_list
}

# * Telescope -------------------------------------------------------------
telescope <- function(theta_init,
                      FUN_tele,
                      fix_index, # = NA if not required
                      x1,
                      x2,
                      x3,
                      y,
                      tau_tele_vec,
                      eps_tele_vec,
                      list_family, # list containing all functions (general and derivatives) for a specific family
                      method_c_tilde, # "table" or "integral"
                      method_c_tilde_deriv, # "finite_diff" or "actual"
                      kappa_omega, # offset value to restrict kappa
                      lambda_beta,
                      lambda_alpha,
                      lambda_nu,
                      initial_step,
                      max_step_it,
                      algorithm,
                      h_kappa,
                      tol,
                      max_it) {
  p1 <- ncol(x1) - 1 # important for getting dimensions of the output matrix
  p2 <- ncol(x2) - 1
  p3 <- ncol(x3) - 1

  t_initial_guess <- theta_init
  t_res_mat <- matrix(NA,
                      nrow = length(eps_tele_vec),
                      ncol = (length(t_initial_guess) + 6)
  )
  colnames(t_res_mat) <- c(
    "epsilon",
    "tau",
    paste0("beta_", 0:p1),
    paste0("alpha_", 0:p2),
    paste0("nu_", 0:p3),
    "plike_val",
    "step_full",
    "it_full",
    "it"
  )
  ## Telescope
  # tau and epsilon telescope
  steps_vec <- length(eps_tele_vec)
  for (i in 1:steps_vec) {
    eps_val <- eps_tele_vec[i]
    tau_val <- tau_tele_vec[i]
    pos <- i

    t_fit <- it_loop(
      theta = t_initial_guess,
      FUN = FUN_tele, # nr_penalty
      tol = tol,
      max_it = max_it,
      fix_index = fix_index,
      x1 = x1,
      x2 = x2,
      x3 = x3,
      y = y,
      tau = tau_val, # telescope
      epsilon = eps_val, # telescope
      list_family = list_family,
      method_c_tilde = method_c_tilde,
      method_c_tilde_deriv = method_c_tilde_deriv,
      kappa_omega = kappa_omega,
      lambda_beta = lambda_beta,
      lambda_alpha = lambda_alpha,
      lambda_nu = lambda_nu,
      initial_step = initial_step,
      max_step_it = max_step_it,
      algorithm = algorithm,
      h_kappa = h_kappa
    )
    t_initial_guess <- t_fit$estimate
    steps_full <- ifelse(any(t_fit$max_step_reached == 1), 1, 0) # 1 if true, 0 if false

    t_res_mat[pos, ] <- c(eps_val, tau_val, t_fit$estimate, t_fit$maximum, steps_full, t_fit$iterations)
  }
  as.data.frame(t_res_mat) # return data.frame instead of matrix
}

# Profile over kappa (nu) ------------------------------------------------------
## nu warm starts (onesided)
nu_warm_starts <- function(theta_largeeps, # initial values to kick off warm_starts - should be for kappa = 2 and largest epsilon in telescope
                           nu_vec_oneside, # either lhs or rhs values, needs to be sorted in order required, nu = log(kappa_true - kappa_omega)
                           x1,
                           x2,
                           x3,
                           y,
                           tau_tele_vec,
                           eps_tele_vec,
                           list_family, # list containing all functions (general and derivatives) for a specific family
                           method_c_tilde,
                           method_c_tilde_deriv,
                           kappa_omega,
                           lambda_beta,
                           lambda_alpha,
                           lambda_nu,
                           initial_step,
                           max_step_it,
                           algorithm,
                           h_kappa,
                           tol,
                           max_it) {
  theta_warmstart <- theta_largeeps # kick off warmstarts, should be for kappa = 2 and largest epsilon in telescope

  length_theta <- length(theta_warmstart)
  nu_index <- length_theta # position of nu

  telefit_list <- list()
  for (i in nu_vec_oneside) {
    pos <- which(nu_vec_oneside == i)
    theta_warmstart[nu_index] <- i # value that is fixed in theta

    telefit <- tryCatch(telescope(
      theta_init = theta_warmstart, # contains fixed value
      FUN_tele = nr_penalty,
      fix_index = nu_index, # nu fixed
      x1,
      x2,
      x3,
      y,
      tau_tele_vec,
      eps_tele_vec,
      list_family,
      method_c_tilde,
      method_c_tilde_deriv,
      kappa_omega,
      lambda_beta,
      lambda_alpha,
      lambda_nu,
      initial_step,
      max_step_it,
      algorithm,
      h_kappa,
      tol,
      max_it
    ),
    error = function(err) NA
    )

    if (class(telefit)[1] != "logical") { # if not NA, then update initial estimates
      # update initial values
      # careful with index - depends on shape of telescope output
      theta_warmstart <- unlist(telefit[1, 2:(length_theta + 1)]) # large epsilon value
      # can use grep to detect as an alternative
      # pos_extract <- c(grep("beta_", colnames(telefit_init)),
      #                  grep("alpha_", colnames(telefit_init)))
      # theta_warmstart <- unlist(telefit_init[1, pos_extract])
    }
    telefit_list[[pos]] <- telefit # save to list
  }
  telefit_list
}

# ** wrapper fit function warm starts
nu_profile_over <- function(theta_init, # should include first nu to optimize
                            nu_vec, # vec to search over for nu, needs to include nu_init
                            x1,
                            x2,
                            x3,
                            y,
                            tau_tele_vec,
                            eps_tele_vec,
                            list_family,
                            method_c_tilde,
                            method_c_tilde_deriv,
                            kappa_omega,
                            lambda_beta,
                            lambda_alpha,
                            lambda_nu,
                            initial_step,
                            max_step_it,
                            algorithm,
                            h_kappa,
                            tol,
                            max_it) {
  nu_index <- length(theta_init)
  nu_init <- theta_init[nu_index]

  stopifnot(nu_init %in% nu_vec) # nu_vec needs to contain nu_init

  # Initial optimization using inital values and nu_init ---
  telefit_init <- telescope(
    theta_init = theta_init, # contains fixed value
    FUN_tele = nr_penalty,
    fix_index = nu_index, # nu fixed
    x1,
    x2,
    x3,
    y,
    tau_tele_vec,
    eps_tele_vec,
    list_family,
    method_c_tilde,
    method_c_tilde_deriv,
    kappa_omega,
    lambda_beta,
    lambda_alpha,
    lambda_nu,
    initial_step,
    max_step_it,
    algorithm,
    h_kappa,
    tol,
    max_it
  )

  ## Extract initial values at large epsilon (row 1) --
  colnames_telefit <- colnames(telefit_init)
  pos_extract <- c(
    grep("beta_", colnames_telefit),
    grep("alpha_", colnames_telefit),
    grep("nu_", colnames_telefit)
  )

  theta_largeeps <- unlist(telefit_init[1, pos_extract]) # large epsilon

  ## Setup LHS and RHS nu vector --
  nu_lhs <- sort(nu_vec[nu_vec < nu_init], decreasing = TRUE)
  nu_rhs <- sort(nu_vec[nu_vec > nu_init])

  ## LHS --
  list_lhs <- nu_warm_starts(
    theta_largeeps = theta_largeeps,
    nu_vec_oneside = nu_lhs,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    tau_tele_vec = tau_tele_vec,
    eps_tele_vec = eps_tele_vec,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    algorithm = algorithm,
    h_kappa = h_kappa,
    tol = tol,
    max_it = max_it
  )

  ## RHS --
  list_rhs <- nu_warm_starts(
    theta_largeeps = theta_largeeps,
    nu_vec_oneside = nu_rhs,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    tau_tele_vec = tau_tele_vec,
    eps_tele_vec = eps_tele_vec,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    algorithm = algorithm,
    h_kappa = h_kappa,
    tol = tol,
    max_it = max_it
  )
  ## Output --
  fit_list <- append(
    list_lhs,
    list_rhs
  )

  fit_list <- append(
    fit_list,
    list(telefit_init)
  )

  ## Remove any NAs --
  na_remove <- which(is.na(fit_list))
  fit_list[na_remove] <- NULL

  do.call(rbind, fit_list) # return data frame of results
}

## Wrapper for fitting and selecting optimum profile over nu --
nu_profile_over_extract <- function(df_profile) {
  # dataframe of estimates for smallest epsilon
  min_eps <- min(df_profile$epsilon)

  df_opt <- subset(df_profile, epsilon == min_eps) # nrow = length(nu_vec)

  # optimum estimates
  max_pos <- which.max(df_opt$plike_val)

  colnames_telefit <- colnames(df_opt)
  pos_extract <- c(
    grep("beta_", colnames_telefit),
    grep("alpha_", colnames_telefit),
    grep("nu_", colnames_telefit)
  )

  theta <- unlist(df_opt[max_pos, pos_extract]) # scaled so may need to unscale

  # return opt_df and theta
  list(
    theta = theta,
    df_opt = df_opt
  )
}

# nlm ---------------------------------------------------------------------
telescope_nlm <- function(theta_init,
                          FUN_nlm,
                          x1,
                          x2,
                          x3,
                          y,
                          tau_tele_vec,
                          eps_tele_vec,
                          list_family, # list containing all functions (general and derivatives) for a specific family
                          method_c_tilde,
                          kappa_omega, # offset value to restrict kappa
                          lambda_beta,
                          lambda_alpha,
                          lambda_nu,
                          iterlim,
                          stepmax) {
  p1 <- ncol(x1) - 1 # important for getting dimensions of the output matrix
  p2 <- ncol(x2) - 1
  p3 <- ncol(x3) - 1

  t_initial_guess <- theta_init
  t_res_mat <- matrix(NA,
                      nrow = length(eps_tele_vec),
                      ncol = (length(t_initial_guess) + 6)
  )
  colnames(t_res_mat) <- c(
    "epsilon",
    "tau",
    paste0("beta_", 0:p1),
    paste0("alpha_", 0:p2),
    paste0("nu_", 0:p3),
    "plike_val",
    "step_full",
    "it_full",
    "it"
  )

  stepmax_lgl <- ifelse(is.na(stepmax), "no_stepmax", "yes_stepmax")
  ## Telescope
  steps_vec <- length(eps_tele_vec)
  for (i in 1:steps_vec) {
    eps_val <- eps_tele_vec[i]
    tau_val <- tau_tele_vec[i]
    pos <- i

    # run depending on stepmax
    switch(stepmax_lgl,
           "yes_stepmax" = {
             t_fit <- suppressWarnings({
               nlm(FUN_nlm, # nlm minimizes -> likelihood_penalized_neg
                   p = t_initial_guess,
                   x1 = x1,
                   x2 = x2,
                   x3 = x3,
                   y = y,
                   tau = tau_val,
                   epsilon = eps_val,
                   list_general = list_family$general,
                   method_c_tilde = method_c_tilde,
                   kappa_omega = kappa_omega,
                   lambda_beta = lambda_beta,
                   lambda_alpha = lambda_alpha,
                   lambda_nu = lambda_nu,
                   iterlim = iterlim,
                   stepmax = stepmax # include stepmax option
               )
             })
           },
           "no_stepmax" = {
             t_fit <- suppressWarnings({
               nlm(FUN_nlm, # nlm minimizes -> likelihood_penalized_neg
                   p = t_initial_guess,
                   x1 = x1,
                   x2 = x2,
                   x3 = x3,
                   y = y,
                   tau = tau_val,
                   epsilon = eps_val,
                   list_general = list_family$general,
                   method_c_tilde = method_c_tilde,
                   kappa_omega = kappa_omega,
                   lambda_beta = lambda_beta,
                   lambda_alpha = lambda_alpha,
                   lambda_nu = lambda_nu,
                   iterlim = iterlim
               )
             })
           }
    )

    t_initial_guess <- t_fit$estimate
    steps_full <- 0 # 1 if true, 0 if false (just set to zero for nlm as not dealing with step halving)

    it_full <- ifelse(t_fit$iterations == iterlim, 1, 0)

    t_res_mat[pos, ] <- c(eps_val, tau_val, t_fit$estimate, -t_fit$minimum, steps_full, it_full, t_fit$iterations) # multiply plike_val by -1 (ie turn minumum into maximum)
  }
  as.data.frame(t_res_mat) # return data.frame instead of matrix
}

# nlm fixed indices -------------------------------------------------------
likelihood_penalized_neg_fixed <- function(theta_estimate,
                                           fix_indices,
                                           fix_values, # match fix_indices
                                           x1,
                                           x2,
                                           x3,
                                           y,
                                           tau,
                                           epsilon,
                                           list_general, # list containing general family functions
                                           method_c_tilde, # "table" or "integral"
                                           kappa_omega, # offset for restricting kappa
                                           lambda_beta,
                                           lambda_alpha,
                                           lambda_nu) {
  theta <- theta_estimate
  for (i in fix_indices) {
    pos <- which(fix_indices == i)
    theta <- append(
      x = theta,
      values = fix_values[pos],
      after = i - 1
    ) # insert value after fix_index - 1 position
  }

  likelihood_penalized_neg(
    theta,
    x1,
    x2,
    x3,
    y,
    tau,
    epsilon,
    list_general, # list containing general family functions
    method_c_tilde, # "table" or "integral"
    kappa_omega, # offset for restricting kappa
    lambda_beta,
    lambda_alpha,
    lambda_nu
  )
}

nlm_fixed <- function(f, # likelihood_penalized_neg_fixed
                      p, # "initial values" should have the fixed values included here
                      fix_indices,
                      x1,
                      x2,
                      x3,
                      y,
                      tau,
                      epsilon,
                      list_general, # list containing general family functions
                      method_c_tilde, # "table" or "integral"
                      kappa_omega, # offset for restricting kappa
                      lambda_beta,
                      lambda_alpha,
                      lambda_nu,
                      iterlim_nlm,
                      stepmax) {
  theta_estimate <- p[-fix_indices] # remove fix_indices value from the vector

  fix_values <- p[fix_indices] # the values that are fixed

  stepmax_lgl <- ifelse(is.na(stepmax), "no_stepmax", "yes_stepmax")

  switch (stepmax_lgl,
          "yes_stepmax" = {
            res <- nlm(f, # likelihood_penalized_neg_fixed - nlm minimizes
                       p = theta_estimate,
                       fix_indices = fix_indices,
                       fix_values = fix_values,
                       x1 = x1,
                       x2 = x2,
                       x3 = x3,
                       y = y,
                       tau = tau,
                       epsilon = epsilon,
                       list_general = list_general,
                       method_c_tilde = method_c_tilde,
                       kappa_omega = kappa_omega,
                       lambda_beta = lambda_beta,
                       lambda_alpha = lambda_alpha,
                       lambda_nu = lambda_nu,
                       iterlim = iterlim_nlm,
                       stepmax = stepmax # include stepmax option
            )
          },
          "no_stepmax" = {
            res <- nlm(f, # likelihood_penalized_neg_fixed - nlm minimizes
                       p = theta_estimate,
                       fix_indices = fix_indices,
                       fix_values = fix_values,
                       x1 = x1,
                       x2 = x2,
                       x3 = x3,
                       y = y,
                       tau = tau,
                       epsilon = epsilon,
                       list_general = list_general,
                       method_c_tilde = method_c_tilde,
                       kappa_omega = kappa_omega,
                       lambda_beta = lambda_beta,
                       lambda_alpha = lambda_alpha,
                       lambda_nu = lambda_nu,
                       iterlim = iterlim_nlm
            )
          }
  )

  # fill in the fixed values so have complete vector returned
  theta <- res$estimate
  for (i in fix_indices) {
    pos <- which(fix_indices == i)
    theta <- append(
      x = theta,
      values = fix_values[pos],
      after = i - 1
    ) # insert value after fix_index - 1 position
  }

  res$estimate <- theta # replace nlm estimates with complete vector (includes fixed)
  res
}

nlm_fixed_unpenalized <- function(p, # "initial values" should have the fixed values included here
                                  fix_indices,
                                  x1,
                                  x2,
                                  x3,
                                  y,
                                  tau,
                                  list_general, # list containing general family functions
                                  method_c_tilde, # "table" or "integral"
                                  kappa_omega, # offset for restricting kappa
                                  iterlim_nlm,
                                  stepmax) {
  theta_estimate <- p[-fix_indices] # remove fix_indices value from the vector

  fix_values <- p[fix_indices] # the values that are fixed

  stepmax_lgl <- ifelse(is.na(stepmax), "no_stepmax", "yes_stepmax")

  switch(stepmax_lgl,
         "yes_stepmax" = {
           res <- nlm(f = likelihood_neg_fixed, # likelihood_neg_fixed - nlm minimizes, unpenalized
                      p = theta_estimate,
                      fix_indices = fix_indices,
                      fix_values = fix_values,
                      x1 = x1,
                      x2 = x2,
                      x3 = x3,
                      y = y,
                      tau = tau,
                      list_general = list_general,
                      method_c_tilde = method_c_tilde,
                      kappa_omega = kappa_omega,
                      iterlim = iterlim_nlm,
                      stepmax = stepmax # include stepmax option
           )
         },
         "no_stepmax" = {
           res <- nlm(f = likelihood_neg_fixed, # likelihood_neg_fixed - nlm minimizes, unpenalized
                      p = theta_estimate,
                      fix_indices = fix_indices,
                      fix_values = fix_values,
                      x1 = x1,
                      x2 = x2,
                      x3 = x3,
                      y = y,
                      tau = tau,
                      list_general = list_general,
                      method_c_tilde = method_c_tilde,
                      kappa_omega = kappa_omega,
                      iterlim = iterlim_nlm
           )
         }
  )

  # fill in the fixed values so have complete vector returned
  theta <- res$estimate
  for (i in fix_indices) {
    pos <- which(fix_indices == i)
    theta <- append(
      x = theta,
      values = fix_values[pos],
      after = i - 1
    ) # insert value after fix_index - 1 position
  }

  res$estimate <- theta # replace nlm estimates with complete vector (includes fixed)
  res
}

telescope_nlm_fixed <- function(theta_init, # contains fixed values
                                fix_indices,
                                x1,
                                x2,
                                x3,
                                y,
                                tau_tele_vec,
                                eps_tele_vec,
                                list_family, # list containing all functions (general and derivatives) for a specific family
                                method_c_tilde,
                                kappa_omega, # offset value to restrict kappa
                                lambda_beta,
                                lambda_alpha,
                                lambda_nu,
                                iterlim,
                                stepmax) {
  p1 <- ncol(x1) - 1 # important for getting dimensions of the output matrix
  p2 <- ncol(x2) - 1
  p3 <- ncol(x3) - 1

  t_initial_guess <- theta_init
  t_res_mat <- matrix(NA,
                      nrow = length(eps_tele_vec),
                      ncol = (length(t_initial_guess) + 6)
  )
  colnames(t_res_mat) <- c(
    "epsilon",
    "tau",
    paste0("beta_", 0:p1),
    paste0("alpha_", 0:p2),
    paste0("nu_", 0:p3),
    "plike_val",
    "step_full",
    "it_full",
    "it"
  )
  ## Telescope
  steps_vec <- length(eps_tele_vec)

  for (i in 1:steps_vec) {
    eps_val <- eps_tele_vec[i]
    tau_val <- tau_tele_vec[i]
    pos <- i

    t_fit <- suppressWarnings({
      nlm_fixed(likelihood_penalized_neg_fixed, # nlm minimizes -> likelihood_penalized_neg
                p = t_initial_guess,
                fix_indices = fix_indices,
                x1 = x1,
                x2 = x2,
                x3 = x3,
                y = y,
                tau = tau_val,
                epsilon = eps_val,
                list_general = list_family$general,
                method_c_tilde = method_c_tilde,
                kappa_omega = kappa_omega,
                lambda_beta = lambda_beta,
                lambda_alpha = lambda_alpha,
                lambda_nu = lambda_nu,
                iterlim = iterlim,
                stepmax = stepmax
      )
    })

    t_initial_guess <- t_fit$estimate
    steps_full <- 0 # 1 if true, 0 if false (just set to zero for nlm as not dealing with step halving)

    it_full <- ifelse(t_fit$iterations == iterlim, 1, 0)

    t_res_mat[pos, ] <- c(eps_val, tau_val, t_fit$estimate, -t_fit$minimum, steps_full, it_full, t_fit$iterations) # multiply plike_val by -1 (ie turn minumum into maximum)
  }
  as.data.frame(t_res_mat) # return data.frame instead of matrix
}

# Fitting Function --------------------------------------------------------
fitting_func_base <- function(x1, # data should be scaled
                              x2,
                              x3,
                              y,
                              optimizer,
                              iterlim_nlm,
                              stepmax_nlm,
                              nu_profile_vec,
                              theta_init, # initial values for scaled data, should contain fixed value for nu if taking that option
                              tau_1,
                              tau_T,
                              epsilon_1,
                              epsilon_T,
                              steps_T,
                              list_family,
                              method_c_tilde,
                              method_c_tilde_deriv,
                              algorithm,
                              h_kappa,
                              kappa_omega,
                              lambda_beta,
                              lambda_alpha,
                              lambda_nu,
                              initial_step,
                              max_step_it,
                              tol,
                              max_it) {
  tau_tele_vec <- rev(exp(seq(log(tau_T), log(tau_1), length = steps_T)))
  eps_tele_vec <- rev(exp(seq(log(epsilon_T), log(epsilon_1), length = steps_T)))

  switch(optimizer,
         "fullopt" = telescope(
           theta_init = theta_init,
           FUN_tele = nr_penalty,
           fix_index = NA, # don't fix any parameters
           x1 = x1,
           x2 = x2,
           x3 = x3,
           y = y,
           tau_tele_vec = tau_tele_vec,
           eps_tele_vec = eps_tele_vec,
           list_family = list_family,
           method_c_tilde = method_c_tilde,
           method_c_tilde_deriv = method_c_tilde_deriv,
           kappa_omega = kappa_omega,
           lambda_beta = lambda_beta,
           lambda_alpha = lambda_alpha,
           lambda_nu = lambda_nu,
           initial_step = initial_step,
           max_step_it = max_step_it,
           algorithm = algorithm,
           h_kappa = h_kappa,
           tol = tol,
           max_it = max_it
         ),
         "nlm" = telescope_nlm(
           theta_init = theta_init,
           FUN_nlm = likelihood_penalized_neg,
           x1 = x1,
           x2 = x2,
           x3 = x3,
           y = y,
           tau_tele_vec = tau_tele_vec,
           eps_tele_vec = eps_tele_vec,
           list_family = list_family,
           method_c_tilde = method_c_tilde,
           kappa_omega = kappa_omega,
           lambda_beta = lambda_beta,
           lambda_alpha = lambda_alpha,
           lambda_nu = lambda_nu,
           iterlim = iterlim_nlm,
           stepmax = stepmax_nlm
         ),
         "fixed_nu" = telescope(
           theta_init = theta_init, # contains true value for nu
           FUN_tele = nr_penalty,
           fix_index = length(theta_init), # fix nu
           x1 = x1,
           x2 = x2,
           x3 = x3,
           y = y,
           tau_tele_vec = tau_tele_vec,
           eps_tele_vec = eps_tele_vec,
           list_family = list_family,
           method_c_tilde = method_c_tilde,
           method_c_tilde_deriv = method_c_tilde_deriv,
           kappa_omega = kappa_omega,
           lambda_beta = lambda_beta,
           lambda_alpha = lambda_alpha,
           lambda_nu = lambda_nu,
           initial_step = initial_step,
           max_step_it = max_step_it,
           algorithm = algorithm,
           h_kappa = h_kappa,
           tol = tol,
           max_it = max_it
         ),
         "profile_nu" = nu_profile_over(
           theta_init = theta_init,
           nu_vec = nu_profile_vec, # nu_vec needs to contain nu_init
           x1 = x1,
           x2 = x2,
           x3 = x3,
           y = y,
           tau_tele_vec = tau_tele_vec,
           eps_tele_vec = eps_tele_vec,
           list_family = list_family,
           method_c_tilde = method_c_tilde,
           method_c_tilde_deriv = method_c_tilde_deriv,
           kappa_omega = kappa_omega,
           lambda_beta = lambda_beta,
           lambda_alpha = lambda_alpha,
           lambda_nu = lambda_nu,
           initial_step = initial_step,
           max_step_it = max_step_it,
           algorithm = algorithm,
           h_kappa = h_kappa,
           tol = tol,
           max_it = max_it
         )
  )
}

fitting_func_nlmtry <- function(x1, # data should be scaled
                                x2,
                                x3,
                                y,
                                optimizer,
                                iterlim_nlm,
                                stepmax_nlm,
                                nu_profile_vec,
                                theta_init, # initial values for scaled data, should contain fixed value for nu if taking that option
                                tau_1,
                                tau_T,
                                epsilon_1,
                                epsilon_T,
                                steps_T,
                                list_family,
                                method_c_tilde,
                                method_c_tilde_deriv,
                                algorithm,
                                h_kappa,
                                kappa_omega,
                                lambda_beta,
                                lambda_alpha,
                                lambda_nu,
                                initial_step,
                                max_step_it,
                                tol,
                                max_it) {
  tele_nlm <- tryCatch(fitting_func_base(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    optimizer = "nlm",
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    nu_profile_vec = NA,
    theta_init = theta_init, # initial values
    tau_1 = tau_1,
    tau_T = tau_T,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = "integrate",
    method_c_tilde_deriv = NA,
    algorithm = NA,
    h_kappa = NA,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  ), error = function(err) NA)

  nlm_method_c_tilde <- "integrate"

  if (class(tele_nlm) == "logical") { # if failed, then try with trapezoid
    tele_nlm <- tryCatch(fitting_func_base(
      x1 = x1,
      x2 = x2,
      x3 = x3,
      y = y,
      optimizer = "nlm",
      iterlim_nlm = iterlim_nlm,
      stepmax_nlm = stepmax_nlm,
      nu_profile_vec = NA,
      theta_init = theta_init, # initial values
      tau_1 = tau_1,
      tau_T = tau_T,
      epsilon_1 = epsilon_1,
      epsilon_T = epsilon_T,
      steps_T = steps_T,
      list_family = list_family,
      method_c_tilde = "trapezoid",
      method_c_tilde_deriv = NA,
      algorithm = NA,
      h_kappa = NA,
      kappa_omega = kappa_omega,
      lambda_beta = lambda_beta,
      lambda_alpha = lambda_alpha,
      lambda_nu = lambda_nu,
      initial_step = initial_step,
      max_step_it = max_step_it,
      tol = tol,
      max_it = max_it
    ), error = function(err) NA)

    nlm_method_c_tilde <- "trapezoid"
  }
  list(
    "tele_nlm" = tele_nlm, # return, can tele data frame or NA
    "nlm_method_c_tilde" = nlm_method_c_tilde
  )
}

fitting_func_fullopttry <- function(x1, # data should be scaled
                                    x2,
                                    x3,
                                    y,
                                    optimizer,
                                    iterlim_nlm,
                                    stepmax_nlm,
                                    nu_profile_vec,
                                    theta_init, # initial values for scaled data, should contain fixed value for nu if taking that option
                                    tau_1,
                                    tau_T,
                                    epsilon_1,
                                    epsilon_T,
                                    steps_T,
                                    list_family,
                                    method_c_tilde,
                                    method_c_tilde_deriv,
                                    algorithm,
                                    h_kappa,
                                    kappa_omega,
                                    lambda_beta,
                                    lambda_alpha,
                                    lambda_nu,
                                    initial_step,
                                    max_step_it,
                                    tol,
                                    max_it) {
  p <- ncol(x1) - 1

  tele_fullopt <- tryCatch(fitting_func_base(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    optimizer = "fullopt",
    iterlim_nlm = NA,
    stepmax_nlm = NA,
    nu_profile_vec = NA,
    theta_init = theta_init, # initial values
    tau_1 = tau_1,
    tau_T = tau_T,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  ), error = function(err) NA)

  method_fit <- "fullopt"

  tele_nu <- tryCatch(unname(extract_theta_plike_val(
    tele_fullopt,
    "fullopt"
  )$theta[(2 * p) + 3]),
  error = function(err) 20
  ) # return 20 if NA so will enter loop below (20 > 10)

  if (tele_nu > kappa_to_nu(10, kappa_omega)) { # if failed or kappa estimated to be greater than 10, then try with nlm
    tele_fullopt <- tryCatch(fitting_func_base(
      x1 = x1,
      x2 = x2,
      x3 = x3,
      y = y,
      optimizer = "nlm", # use nlm
      iterlim_nlm = iterlim_nlm,
      stepmax_nlm = stepmax_nlm,
      nu_profile_vec = NA,
      theta_init = theta_init, # initial values
      tau_1 = tau_1,
      tau_T = tau_T,
      epsilon_1 = epsilon_1,
      epsilon_T = epsilon_T,
      steps_T = steps_T,
      list_family = list_family,
      method_c_tilde = method_c_tilde,
      method_c_tilde_deriv = NA,
      algorithm = NA,
      h_kappa = NA,
      kappa_omega = kappa_omega,
      lambda_beta = lambda_beta,
      lambda_alpha = lambda_alpha,
      lambda_nu = lambda_nu,
      initial_step = initial_step,
      max_step_it = max_step_it,
      tol = tol,
      max_it = max_it
    ), error = function(err) NA)

    method_fit <- "nlm"
  }
  list(
    "tele_fullopt" = tele_fullopt, # return, can tele data frame or NA
    "method_fit" = method_fit
  )
}

fitting_func_twostep_nlm <- function(x1, # data should be scaled
                                     x2,
                                     x3,
                                     y,
                                     optimizer,
                                     iterlim_nlm,
                                     stepmax_nlm,
                                     theta_init, # initial values for scaled data, should contain fixed value for nu if taking that option
                                     tau_1,
                                     tau_T,
                                     epsilon_1,
                                     epsilon_T,
                                     steps_T,
                                     list_family,
                                     method_c_tilde,
                                     kappa_omega,
                                     lambda_beta,
                                     lambda_alpha,
                                     lambda_nu,
                                     less_than_value) {
  # Fit penalized nlm telescope, then using these estimates to get unpenalized with variables fixed to zero
  initial_fit <- fitting_func_base(x1 = x1,
                                   x2 = x2,
                                   x3 = x3,
                                   y = y,
                                   optimizer = "nlm",
                                   iterlim_nlm = iterlim_nlm,
                                   stepmax_nlm = stepmax_nlm,
                                   nu_profile_vec = NA,
                                   theta_init = theta_init, # initial values
                                   tau_1 = tau_1,
                                   tau_T = tau_T,
                                   epsilon_1 = epsilon_1,
                                   epsilon_T = epsilon_T,
                                   steps_T = steps_T,
                                   list_family = list_family,
                                   method_c_tilde = method_c_tilde,
                                   method_c_tilde_deriv = NA,
                                   algorithm = NA,
                                   h_kappa = NA,
                                   kappa_omega = kappa_omega,
                                   lambda_beta = lambda_beta,
                                   lambda_alpha = lambda_alpha,
                                   lambda_nu = lambda_nu,
                                   initial_step = NA,
                                   max_step_it = NA,
                                   tol = NA,
                                   max_it = NA)
  initial_est <- extract_theta_plike_val(initial_fit,
                                         "nlm")$theta

  initial_round <- initial_est
  initial_round[abs(initial_round) < less_than_value] <- 0

  p <- ncol(x1) - 1
  int_remove <- c(1, (p + 2), (2 * p) + 3) # position of intercepts

  fix_indices <- which(initial_round == 0)
  fix_indices <- fix_indices[!(fix_indices %in% int_remove)]

  final_res <- nlm_fixed_unpenalized(p = initial_round, # includes 0s in positions
                                     fix_indices = fix_indices, # index of coefs set to zero
                                     x1 = x1,
                                     x2 = x2,
                                     x3 = x3,
                                     y = y,
                                     tau = tau_T,
                                     list_general = list_family$general,
                                     method_c_tilde = method_c_tilde,
                                     kappa_omega = kappa_omega,
                                     iterlim_nlm = iterlim_nlm,
                                     stepmax = stepmax_nlm)
  final_res$estimate # return theta
}

second_step_nlm_func <- function(theta_tele, # theta from telescope, unscaled
                                 x, # raw, unscaled
                                 x_sd, # for conversion
                                 y,
                                 iterlim_nlm,
                                 stepmax_nlm,
                                 tau,
                                 list_general,
                                 method_c_tilde,
                                 kappa_omega,
                                 less_than_value) {
  x_scale <- map2_dfc(
    .x = as_tibble(x),
    .y = c(1, x_sd),
    ~ {
      .x / .y
    }
  ) %>%
    as.matrix()

  x_sd_theta <- c(1, x_sd, 1, x_sd, 1)

  theta_tele_scale <- theta_tele * x_sd_theta # scaled estimates

  theta_tele_scale_round <- theta_tele_scale
  theta_tele_scale_round[abs(theta_tele_scale_round) < less_than_value] <- 0

  p <- ncol(x_scale) - 1
  n <- nrow(x_scale)

  int_remove <- c(1, (p + 2), (2 * p) + 3)

  fix_indices <- which(theta_tele_scale_round == 0)
  fix_indices <- fix_indices[!(fix_indices %in% int_remove)]

  fitres <- nlm_fixed_unpenalized(
    p = theta_tele_scale_round,
    fix_indices = fix_indices,
    x1 = x_scale,
    x2 = x_scale,
    x3 = as.matrix(rep(1, n)),
    y = y,
    tau = tau,
    list_general = list_general,
    method_c_tilde = method_c_tilde,
    kappa_omega = kappa_omega,
    iterlim_nlm = iterlim_nlm,
    stepmax = stepmax_nlm
  )
  theta_return <- fitres$estimate / x_sd_theta
  unname(theta_return) # return unscale estimates
}

second_step_nlm_func_parallel_single <- function(list_input, # contains theta, x, x_sd, y
                                                 iterlim_nlm,
                                                 stepmax_nlm,
                                                 tau,
                                                 list_general,
                                                 method_c_tilde,
                                                 kappa_omega,
                                                 less_than_value) {
  theta_est <- tryCatch(second_step_nlm_func(
    theta_tele = list_input$theta_est,
    x = list_input$x,
    x_sd = list_input$x_sd,
    y = list_input$y,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    tau = tau,
    list_general = list_general,
    method_c_tilde = method_c_tilde,
    kappa_omega = kappa_omega,
    less_than_value = less_than_value
  ),
  error = function(err) NA)

  list(
    "theta_est" = theta_est,
    "x" = list_input$x,
    "y" = list_input$y
  ) # return useful for information matrices
}

second_step_nlm_func_parallel_multiple <- function(list_input_all, # will contain lists of theta, x, x_sd, y
                                                   iterlim_nlm,
                                                   stepmax_nlm,
                                                   tau,
                                                   list_general,
                                                   method_c_tilde,
                                                   kappa_omega,
                                                   less_than_value) {
  parLapply(
    cl,
    list_input_all,
    second_step_nlm_func_parallel_single,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    tau = tau,
    list_general = list_general,
    method_c_tilde = method_c_tilde,
    kappa_omega = kappa_omega,
    less_than_value = less_than_value
  )
}

# MoM: method of moments initial values ----
kurt_func <- function(kappa, y) {
  (((gamma(5 / kappa) * gamma(1 / kappa)) / ((gamma(3 / kappa))^2)) - moments::kurtosis(y))^2
}

initial_values_MoM <- function(y) {
  kappa_est <- optimize(kurt_func, c(0.1, 10), y = y)$minimum

  phi_est <- (var(as.numeric(y)) * gamma(1 / kappa_est)) / gamma(3 / kappa_est)

  mu_est <- mean(y)

  list(
    "mu_est" = mu_est,
    "phi_est" = phi_est,
    "kappa_est" = kappa_est
  )
}

# Extract theta & plike_val -----------------------------------------------
extract_theta_plike_val <- function(fit_res,
                                    optimizer) {
  extract_type <- ifelse(optimizer == "profile_nu", "profile", "tele") # what sort of output do we have

  switch(extract_type,
         "tele" = {
           plike_val <- fit_res %>%
             select(plike_val) %>%
             slice_tail() %>% # take last row
             pull()
           theta <- fit_res %>%
             select(contains(c("beta", "alpha", "nu"))) %>%
             slice_tail() %>% # take last row
             unlist()
         },
         "profile" = {
           extract_profile <- nu_profile_over_extract(fit_res)
           tele_profile <- extract_profile$df_opt # df containing optimum estimates for each nu
           theta <- extract_profile$theta # optimium estimates based on nu that give maximum penalized likelihood value

           plike_val <- tele_profile %>%
             slice_max(plike_val) %>%
             pull(plike_val)
         }
  )
  list(
    "theta" = theta,
    "plike_val" = plike_val
  )
}

# Initial values function ------------------------------------------------
initial_values_func <- function(x1, # data should be scaled already
                                x2,
                                x3,
                                y,
                                method_initial_values,
                                other_initial_values_vec, # can manually input vector of initial values
                                kappa_omega,
                                fix_nu_value_init,
                                tau_1,
                                epsilon_1,
                                iterlim_nlm,
                                stepmax_nlm,
                                list_family,
                                method_c_tilde,
                                method_c_tilde_deriv,
                                algorithm,
                                h_kappa,
                                lambda_beta,
                                lambda_alpha,
                                lambda_nu,
                                initial_step,
                                max_step_it,
                                tol,
                                max_it) {
  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1
  p3 <- ncol(x3) - 1

  switch(method_initial_values,
         "lm" = { # kappa init = 2 here for normal values
           lm_fit <- lm(y ~ x1[, -1]) # remove intercept column
           lm_coef_sig <- c(
             unname(lm_fit$coefficients),
             log((summary(lm_fit)$sigma)^2)
           )
           theta_init <- c(lm_coef_sig, rep(0, p2), kappa_to_nu(
             kappa = 2,
             kappa_omega = kappa_omega
           )) # kappa = exp(nu_0) + omega, so then nu_0 = log(kappa - omega)
           theta_init
         },
         "MoM" = { # method of moments
           est_MoM <- initial_values_MoM(y)

           theta_init <- c(
             est_MoM$mu_est, rep(0, p1),
             log(est_MoM$phi_est), rep(0, p2),
             kappa_to_nu(est_MoM$kappa_est, kappa_omega)
           )
           theta_init
         },
         "lm_MoM" = {
           lm_fit <- lm(y ~ x1[, -1]) # remove intercept column
           lm_coef_sig <- c(
             unname(lm_fit$coefficients),
             log((summary(lm_fit)$sigma)^2)
           )
           est_MoM <- initial_values_MoM(y)

           theta_init <- c(lm_coef_sig, rep(0, p2), kappa_to_nu(est_MoM$kappa_est, kappa_omega)) # kappa = exp(nu_0) + omega, so then nu_0 = log(kappa - omega)
           theta_init
         },
         "fixed_nu_fullopt" = { # get initial values for a fixed kappa/nu
           # get lm values to kick off
           lm_fit <- lm(y ~ x1[, -1]) # remove intercept column
           lm_coef_sig <- c(
             unname(lm_fit$coefficients),
             log((summary(lm_fit)$sigma)^2)
           )
           theta_init_lm <- c(lm_coef_sig, rep(0, p2), kappa_to_nu(2, kappa_omega)) # kappa = exp(nu_0) + omega, so then nu_0 = log(kappa - omega)

           # set nu value to desired value
           theta_init_fixed <- theta_init_lm
           theta_init_fixed[length(theta_init_lm)] <- fix_nu_value_init

           it_loop(
             theta = theta_init_fixed,
             FUN = nr_penalty,
             tol = tol,
             max_it = max_it,
             fix_index = length(theta_init_fixed),
             x1 = x1,
             x2 = x2,
             x3 = x3,
             y = y,
             tau = tau_1,
             epsilon = epsilon_1,
             list_family = list_family,
             method_c_tilde = method_c_tilde,
             method_c_tilde_deriv = method_c_tilde_deriv,
             kappa_omega = kappa_omega,
             lambda_beta = lambda_beta,
             lambda_alpha = lambda_alpha,
             lambda_nu = lambda_nu,
             initial_step = initial_step,
             max_step_it = max_step_it,
             algorithm = algorithm,
             h_kappa = h_kappa
           )$estimate
         },
         "nlm" = {
           # get lm values to kick off
           lm_fit <- lm(y ~ x1[, -1]) # remove intercept column
           lm_coef_sig <- c(
             unname(lm_fit$coefficients),
             log((summary(lm_fit)$sigma)^2)
           )
           theta_init_lm <- c(lm_coef_sig, rep(0, p2), kappa_to_nu(2, kappa_omega)) # kappa = exp(nu_0) + omega, so then nu_0 = log(kappa - omega)

           nlm(likelihood_penalized_neg,
               p = theta_init_lm,
               x1 = x1,
               x2 = x2,
               x3 = x3,
               y = y,
               tau = tau_1,
               epsilon = epsilon_1,
               list_general = list_family$general,
               method_c_tilde = method_c_tilde,
               kappa_omega = kappa_omega,
               lambda_beta = lambda_beta,
               lambda_alpha = lambda_alpha,
               lambda_nu = lambda_nu,
               iterlim = iterlim_nlm,
               stepmax = stepmax_nlm
           )$estimate
         },
         "fixed_nu_nlm" = {
           # get lm values to kick off
           lm_fit <- lm(y ~ x1[, -1]) # remove intercept column
           lm_coef_sig <- c(
             unname(lm_fit$coefficients),
             log((summary(lm_fit)$sigma)^2)
           )
           theta_init_lm <- c(lm_coef_sig, rep(0, p2), kappa_to_nu(2, kappa_omega)) # kappa = exp(nu_0) + omega, so then nu_0 = log(kappa - omega)

           # set nu value to desired value
           theta_init_fixed <- theta_init_lm
           theta_init_fixed[length(theta_init_lm)] <- fix_nu_value_init

           nlm_fixed(likelihood_penalized_neg_fixed, # nlm minimizes -> likelihood_penalized_neg
                     p = theta_init_fixed,
                     fix_indices = length(theta_init_fixed),
                     x1 = x1,
                     x2 = x2,
                     x3 = x3,
                     y = y,
                     tau = tau_1,
                     epsilon = epsilon_1,
                     list_general = list_family$general,
                     method_c_tilde = method_c_tilde,
                     kappa_omega = kappa_omega,
                     lambda_beta = lambda_beta,
                     lambda_alpha = lambda_alpha,
                     lambda_nu = lambda_nu,
                     iterlim = iterlim_nlm,
                     stepmax = stepmax_nlm
           )$estimate
         },
         "other" = other_initial_values_vec # manually input vector
  )
}

# wrapper fitting function ------------------------------------------------
fitting_func_complete <- function(x1_raw, # unscaled data
                                  x2_raw,
                                  x3_raw,
                                  y,
                                  method_initial_values,
                                  other_initial_values_vec,
                                  fix_nu_value_init,
                                  optimizer,
                                  iterlim_nlm,
                                  stepmax_nlm,
                                  nu_profile_vec,
                                  theta_init,
                                  tau_1,
                                  tau_T,
                                  epsilon_1,
                                  epsilon_T,
                                  steps_T,
                                  list_family,
                                  method_c_tilde,
                                  method_c_tilde_deriv,
                                  algorithm,
                                  h_kappa,
                                  kappa_omega,
                                  lambda_beta,
                                  lambda_alpha,
                                  lambda_nu,
                                  initial_step,
                                  max_step_it,
                                  tol,
                                  max_it) {
  # Scale data depending on setting -----------------------
  # mpr_equal if x1 and x2 identical ----
  # mpr_diff if x1 not identical to x2 and ncol(x2) > 1 ----
  # spr if ncol(x2 == 1) ----
  stopifnot(ncol(x3_raw) == 1)

  n <- nrow(x1_raw)

  identical_lgl <- identical(x1_raw, x2_raw)
  ncol_x2_raw <- ncol(x2_raw)

  if (identical_lgl == TRUE) {
    type <- "mpr_equal"
  } else if (identical_lgl == FALSE & ncol_x2_raw > 1) {
    type <- "mpr_diff"
  } else if (ncol_x2_raw == 1) {
    type <- "spr"
  }

  switch(type,
         "mpr_equal" = {
           ## Scale X ----
           x_scale <- scale(x1_raw[, -1],
                            center = FALSE,
                            scale = apply(x1_raw[, -1], 2, sd)
           )
           x_scale <- cbind(rep(1, n), x_scale) # column of 1's for intercept
           x_sd <- apply(x1_raw[, -1], 2, sd) # save sd for later to transform back

           x_sd_theta <- c(1, x_sd, 1, x_sd, 1)

           x1 <- x_scale
           x2 <- x_scale
           x3 <- x3_raw
         },
         "mpr_diff" = {
           ## Scale both X matrices ---
           x1_scale <- scale(x1_raw[, -1],
                             center = FALSE,
                             scale = apply(x1_raw[, -1], 2, sd)
           )
           x1_scale <- cbind(rep(1, n), x1_scale) # column of 1's for intercept
           x1_sd <- apply(x1_raw[, -1], 2, sd) # save sd for later to transform back

           x2_scale <- scale(x2_raw[, -1],
                             center = FALSE,
                             scale = apply(x2_raw[, -1], 2, sd)
           )
           x2_scale <- cbind(rep(1, n), x2_scale) # column of 1's for intercept
           x2_sd <- apply(x2_raw[, -1], 2, sd) # save sd for later to transform back

           x_sd_theta <- c(1, x1_sd, 1, x2_sd, 1)
           x1 <- x1_scale
           x2 <- x2_scale
           x3 <- x3_raw
         },
         "spr" = {
           ## Scale X ----
           x_scale <- scale(x1_raw[, -1],
                            center = FALSE,
                            scale = apply(x1_raw[, -1], 2, sd)
           )
           x_scale <- cbind(rep(1, n), x_scale) # column of 1's for intercept
           x_sd <- apply(x1_raw[, -1], 2, sd) # save sd for later to transform back

           x_sd_theta <- c(1, x_sd, 1, 1)

           x1 <- x_scale
           x2 <- x2_raw
           x3 <- x3_raw
         }
  )

  # Initial Values ----
  theta_init <- initial_values_func(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    method_initial_values = method_initial_values,
    other_initial_values_vec = other_initial_values_vec,
    kappa_omega = kappa_omega,
    fix_nu_value_init = fix_nu_value_init,
    tau = tau_1,
    epsilon_1 = epsilon_1,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  )
  ## Fit ----
  fit_res <- fitting_func_base(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    optimizer = optimizer,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    nu_profile_vec = nu_profile_vec,
    theta_init = theta_init,
    tau_1 = tau_1,
    tau_T = tau_T,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  )
  ## Extract estimate ----
  extract_res <- extract_theta_plike_val(
    fit_res = fit_res,
    optimizer = optimizer
  )
  theta_scale <- extract_res$theta
  theta <- as.vector(theta_scale / x_sd_theta)

  list(
    "theta" = theta,
    "plike_val" = extract_res$plike_val,
    "call" = c(
      method_initial_values,
      optimizer,
      method_c_tilde
    )
  )
}

# Single Sim --------------------------------------------------------------
single_sim <- function(i,
                       n,
                       beta_true,
                       alpha_true,
                       kappa_true, # kappa = exp(nu_0) + omega
                       binary_positions,
                       kappa_omega,
                       method_initial_values,
                       other_initial_values_vec,
                       fix_nu_value_init,
                       tau_1,
                       tau_T,
                       epsilon_1,
                       epsilon_T,
                       steps_T,
                       list_family,
                       method_c_tilde,
                       method_c_tilde_deriv,
                       algorithm,
                       h_kappa,
                       iterlim_nlm,
                       stepmax_nlm,
                       lambda_beta,
                       lambda_alpha,
                       lambda_nu,
                       initial_step,
                       max_step_it,
                       tol,
                       max_it) {
  stopifnot(length(kappa_true) == 1) # just single value
  stopifnot(kappa_true > kappa_omega)

  # Data with binary covariates? ----
  run_binary <- ifelse(is.na(binary_positions[1]), "no_binary", "yes_binary")

  switch(run_binary,
         "no_binary" = {
           raw_data <- data_sim(
             n = n,
             beta_true = beta_true,
             alpha_true = alpha_true,
             kappa_true = kappa_true,
             tau = tau_T # run for smallest tau
           )
         },
         "yes_binary" = {
           raw_data <- data_sim_binary(
             n = n,
             beta_true = beta_true,
             alpha_true = alpha_true,
             kappa_true = kappa_true,
             tau = tau_T, # run for smallest tau
             binary_positions = binary_positions
           )
         }
  )

  nu_true <- kappa_to_nu(kappa_true, kappa_omega)

  x <- raw_data$x
  y <- raw_data$y
  p <- ncol(x) - 1 # not including intercept


  names_x <- paste0("X_", 1:p)
  ## Scale X ----
  x_scale <- scale(x[, -1],
                   center = FALSE,
                   scale = apply(x[, -1], 2, sd)
  )
  x_scale <- cbind(rep(1, n), x_scale) # column of 1's for intercept
  x_sd <- apply(x[, -1], 2, sd) # save sd for later to transform back
  x_sd_theta <- c(1, x_sd, 1, x_sd, 1)

  x1 <- x_scale
  x2 <- x_scale
  x3 <- as.matrix(rep(1, n))

  ## Initial values ----
  theta_init <- initial_values_func(
    x1 = x1, # data should be scaled already
    x2 = x2,
    x3 = x3,
    y = y,
    method_initial_values = method_initial_values,
    other_initial_values_vec = other_initial_values_vec, # can manually input vector of initial values
    kappa_omega = kappa_omega,
    fix_nu_value_init = fix_nu_value_init,
    tau_1 = tau_1,
    epsilon_1 = epsilon_1,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  )

  # Fits --------------------------------------------------------------------
  time_nlm <- system.time(tele_nlm <- tryCatch(fitting_func_base(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    optimizer = "nlm",
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    nu_profile_vec = NA,
    theta_init = theta_init, # initial values
    tau_1 = tau_1,
    tau_T = tau_T,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = NA,
    algorithm = NA,
    h_kappa = NA,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  ), error = function(err) NA))

  extract_nlm <- tryCatch(extract_theta_plike_val(
    tele_nlm,
    "nlm"
  ),
  error = function(err) NA
  )
  theta_nlm <- tryCatch(as.vector(extract_nlm$theta / x_sd_theta),
                        error = function(err) NA
  ) # unscale
  plike_val_nlm <- tryCatch(extract_nlm$plike_val, error = function(err) NA)

  ## Output ----
  list(
    "data_info" = list(
      "x" = x,
      "y" = y,
      "x_sd" = x_sd,
      "x_sd_theta" = x_sd_theta,
      "theta_true" = c(beta_true, alpha_true, kappa_to_nu(kappa_true, kappa_omega)), # nu_true
      "epsilon_info" = c(epsilon_1, epsilon_T, steps_T),
      "binary_positions" = binary_positions
    ),
    "time_fits" = list(
      "time_nlm" = time_nlm
    ),
    "fits" = list(
      "tele_nlm" = tele_nlm
    ),
    "thetas" = list(
      "theta_nlm" = theta_nlm
    ),
    "plike_vals" = list(
      "plike_val_nlm" = plike_val_nlm
    )
  )
}

# Multiple Sims -----------------------------------------------------------
multiple_sims <- function(nsims,
                          n,
                          beta_true,
                          alpha_true,
                          kappa_true, # kappa = exp(nu_0) + omega
                          binary_positions,
                          kappa_omega,
                          method_initial_values,
                          other_initial_values_vec,
                          fix_nu_value_init,
                          tau_1,
                          tau_T,
                          epsilon_1,
                          epsilon_T,
                          steps_T,
                          list_family,
                          method_c_tilde,
                          method_c_tilde_deriv,
                          algorithm,
                          h_kappa,
                          iterlim_nlm,
                          stepmax_nlm,
                          lambda_beta,
                          lambda_alpha,
                          lambda_nu,
                          initial_step,
                          max_step_it,
                          tol,
                          max_it) {
  simres <- parLapply(
    cl,
    1:nsims,
    single_sim,
    n = n,
    beta_true = beta_true,
    alpha_true = alpha_true,
    kappa_true = kappa_true, # kappa = exp(nu_0) + omega
    binary_positions = binary_positions,
    kappa_omega = kappa_omega,
    method_initial_values = method_initial_values,
    other_initial_values_vec = other_initial_values_vec,
    fix_nu_value_init = fix_nu_value_init,
    tau_1 = tau_1,
    tau_T = tau_T,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  )
  simres # output length(nsims)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Slice -------------------------------------------------------------------
# * fixed coef with warm starts throughout fixed values -------------------
fix_coef_warmstarts_telefits <- function(fix_index,
                                         fix_values_oneside, # either lhs or rhs vector
                                         theta_largeeps, # optimum values for large epsilon (scaled)
                                         x1, # data should be scaled
                                         x2,
                                         x3,
                                         y,
                                         optimizer,
                                         iterlim_nlm,
                                         stepmax_nlm,
                                         tau_tele_vec,
                                         eps_tele_vec,
                                         list_family,
                                         method_c_tilde,
                                         method_c_tilde_deriv,
                                         kappa_omega,
                                         lambda_beta,
                                         lambda_alpha,
                                         lambda_nu,
                                         initial_step,
                                         max_step_it,
                                         algorithm,
                                         h_kappa,
                                         tol,
                                         max_it) {
  theta_warmstart <- theta_largeeps # kick off warmstarts
  length_theta <- length(theta_largeeps)

  telefit_list <- list()
  for (i in fix_values_oneside) {
    theta_warmstart[fix_index] <- i # replace with new fixed value
    pos <- which(fix_values_oneside == i)


    switch(optimizer,
           "fullopt_fixed" = {
             telefit <- tryCatch(telescope(
               theta_init = theta_warmstart,
               FUN_tele = nr_penalty,
               fix_index = fix_index,
               x1 = x1,
               x2 = x2,
               x3 = x3,
               y = y,
               tau_tele_vec = tau_tele_vec,
               eps_tele_vec = eps_tele_vec,
               list_family = list_family,
               method_c_tilde = method_c_tilde,
               method_c_tilde_deriv = method_c_tilde_deriv,
               kappa_omega = kappa_omega,
               lambda_beta = lambda_beta,
               lambda_alpha = lambda_alpha,
               lambda_nu = lambda_nu,
               initial_step = initial_step,
               max_step_it = max_step_it,
               algorithm = algorithm,
               h_kappa = h_kappa,
               tol = tol,
               max_it = max_it
             ),
             error = function(err) NA
             ) # return NA if error
           },
           "nlm_fixed" = {
             telefit <- tryCatch(telescope_nlm_fixed(
               theta_init = theta_warmstart, # updated estimates
               fix_indices = fix_index, # index for fixed coef
               x1 = x1,
               x2 = x2,
               x3 = x3,
               y = y,
               tau_tele_vec = tau_tele_vec,
               eps_tele_vec = eps_tele_vec,
               list_family = list_family,
               method_c_tilde = method_c_tilde,
               kappa_omega = kappa_omega,
               lambda_beta = lambda_beta,
               lambda_alpha = lambda_alpha,
               lambda_nu = lambda_nu,
               iterlim = iterlim_nlm,
               stepm = stepmax_nlm
             ),
             error = function(err) NA
             ) # return NA if error
           }
    )

    if (class(telefit)[1] != "logical") {
      # if not NA, then update initial estimates
      # update initial values
      # careful with index - depends on shape of telescope output
      theta_warmstart <- unlist(telefit[1, 2:(length_theta + 1)]) # large epsilon value
      # can use grep to detect as an alternative
      # pos_extract <- c(grep("beta_", colnames(telefit_init)),
      #                  grep("alpha_", colnames(telefit_init)))
      # theta_warmstart <- unlist(telefit_init[1, pos_extract])
    }
    telefit_list[[pos]] <- telefit
  }
  telefit_list
}

# Select range want to see values ----
fix_coef_optimizer_range_mpr <- function(fix_index,
                                         run_range,
                                         seq_length, # either side (+ 1)
                                         telemat_opt, # scaled
                                         x1, # data should be scaled
                                         x2,
                                         x3,
                                         y,
                                         kappa_omega,
                                         optimizer,
                                         iterlim_nlm,
                                         stepmax_nlm,
                                         tau_1,
                                         tau_T,
                                         epsilon_1,
                                         epsilon_T,
                                         steps_T,
                                         list_family,
                                         method_c_tilde,
                                         method_c_tilde_deriv,
                                         algorithm,
                                         h_kappa,
                                         lambda_beta,
                                         lambda_alpha,
                                         lambda_nu,
                                         initial_step,
                                         max_step_it,
                                         tol,
                                         max_it) {
  eps_tele_vec <- rev(exp(seq(log(epsilon_T), log(epsilon_1), length = steps_T)))
  tau_tele_vec <- rev(exp(seq(log(tau_T), log(tau_1), length = steps_T)))

  theta_smalleps <- telemat_opt %>% # optimum values at smallest epsilon
    select(contains(c("beta", "alpha", "nu"))) %>%
    slice_tail() %>% # last row = small epsilon
    unlist()

  theta_largeeps <- telemat_opt %>% # optimium values at largest epsilon
    select(contains(c("beta", "alpha", "nu"))) %>%
    slice(1) %>% # first row = large epsilon
    unlist()

  fix_opt <- theta_smalleps[fix_index] # optimium value for fixed coef
  stopifnot(fix_opt > run_range[1] & fix_opt < run_range[2]) # range contains optimium value?

  # * Fix values sequence ----
  seq_values <- c(
    seq(
      from = run_range[1],
      to = run_range[2],
      length = seq_length
    ),
    unname(fix_opt)
  )
  fix_values_lhs <- seq_values[which(seq_values <= fix_opt)]
  fix_values_rhs <- seq_values[which(seq_values > fix_opt)] # don't include opt value

  fix_values_list_init <- list(
    fix_values_lhs,
    fix_values_rhs
  )

  # ** Contain zero? ----
  changesign <- map_dbl(fix_values_list_init, ~ {
    vec <- unique(sign(.x))
    changesign <- ifelse(length(vec) == 2, 1, 0) # true if change sign (length vec == 2), false if not
  })
  fix_values_list <- map2(
    .x = fix_values_list_init,
    .y = changesign,
    ~ {
      if (.y == 1) { # if changesign == true, then include 0
        vec <- c(.x, 0)
        # vec[-which(vec == fix_opt)] # if include 0, remove optimum value (will fit on the other side)
      } else {
        .x # == false, leave vector as is
      }
    }
  ) %>%
    map2(.x = ., .y = c(
      TRUE, # sort lhs by decreasing, rhs increasing
      FALSE
    ), ~ {
      sort(.x, decreasing = .y)
    })

  # * Optimization lhs ----
  telefit_lhs <- fix_values_list[[1]] %>%
    fix_coef_warmstarts_telefits(
      fix_index = fix_index,
      fix_values_oneside = ., # fix_values
      theta_largeeps = theta_largeeps, # large epsilon
      x1 = x1, # data should be scaled
      x2 = x2,
      x3 = x3,
      y = y,
      optimizer = optimizer,
      iterlim_nlm = iterlim_nlm,
      stepmax_nlm = stepmax_nlm,
      tau_tele_vec = tau_tele_vec,
      eps_tele_vec = eps_tele_vec,
      list_family = list_family,
      method_c_tilde = method_c_tilde,
      method_c_tilde_deriv = method_c_tilde_deriv,
      kappa_omega = kappa_omega,
      lambda_beta = lambda_beta,
      lambda_alpha = lambda_alpha,
      lambda_nu = lambda_nu,
      initial_step = initial_step,
      max_step_it = max_step_it,
      algorithm = algorithm,
      h_kappa = h_kappa,
      tol = tol,
      max_it = max_it
    ) # list of length(fix_values) of telemat results

  # * Get theta_large eps from lhs to use as initial for rhs
  theta_largeeps_lhs <- telefit_lhs[[1]] %>% # will be 1st in list
    select(contains(c("beta", "alpha", "nu"))) %>%
    slice(1) %>% # first row = large epsilon
    unlist()

  # Optimization rhs ----
  telefit_rhs <- fix_values_list[[2]] %>%
    fix_coef_warmstarts_telefits(
      fix_index = fix_index,
      fix_values_oneside = ., # fix_values
      theta_largeeps = theta_largeeps_lhs, # large epsilon from left hand side
      x1 = x1, # data should be scaled
      x2 = x2,
      x3 = x3,
      y = y,
      optimizer = optimizer,
      iterlim_nlm = iterlim_nlm,
      stepmax_nlm = stepmax_nlm,
      tau_tele_vec = tau_tele_vec,
      eps_tele_vec = eps_tele_vec,
      list_family = list_family,
      method_c_tilde = method_c_tilde,
      method_c_tilde_deriv = method_c_tilde_deriv,
      kappa_omega = kappa_omega,
      lambda_beta = lambda_beta,
      lambda_alpha = lambda_alpha,
      lambda_nu = lambda_nu,
      initial_step = initial_step,
      max_step_it = max_step_it,
      algorithm = algorithm,
      h_kappa = h_kappa,
      tol = tol,
      max_it = max_it
    ) # list of length(fix_values) of telemat results

  telefit_list <- list(
    telefit_lhs,
    telefit_rhs
  )

  # ** Number of NAs in LHS, RHS ----
  any_na <- telefit_list %>%
    map(~ map_lgl(.x, ~ {
      class(.x)[1] == "logical" # any na?
    })) %>%
    map_dbl(sum)

  seq_length_total <- sum(
    length(telefit_list[[1]]),
    length(telefit_list[[2]])
  )
  # ** Combine to large df ----
  telefit_df <- telefit_list %>%
    flatten() %>%
    list_na_remove() %>%
    # remove any NAs
    map2(
      .x = ., .y = seq(1:((seq_length_total) - sum(any_na))), # seq_length - any_na
      ~ {
        .x %>%
          add_column(
            id = !!.y,
            .before = 1
          )
      }
    ) %>%
    data.table::rbindlist()

  # * Return ----
  list(
    "telefit_df" = telefit_df,
    "any_na_sides" = any_na
  )
}

# ** Remove NAs from list -------------------------------------------------
list_na_remove <- function(list1) { # remove NAs from list
  pos <- which(is.na(list1))
  list1[pos] <- NULL
  list1
}

# Slice but keep nu fixed at a value as well ------------------------------
fix_coef_warmstarts_telefits_fix_nu <- function(fix_index,
                                                fix_values_oneside, # either lhs or rhs vector
                                                nu_index,
                                                nu_fix_value,
                                                theta_largeeps, # optimum values for large epsilon (scaled)
                                                x1, # data should be scaled
                                                x2,
                                                x3,
                                                y,
                                                optimizer,
                                                iterlim_nlm,
                                                stepmax_nlm,
                                                tau_tele_vec,
                                                eps_tele_vec,
                                                list_family,
                                                method_c_tilde,
                                                method_c_tilde_deriv,
                                                kappa_omega,
                                                lambda_beta,
                                                lambda_alpha,
                                                lambda_nu,
                                                initial_step,
                                                max_step_it,
                                                algorithm,
                                                h_kappa,
                                                tol,
                                                max_it) {
  theta_warmstart <- theta_largeeps # kick off warmstarts
  length_theta <- length(theta_largeeps)

  # Fix nu to a certain value ----
  theta_warmstart[nu_index] <- nu_fix_value

  telefit_list <- list()
  for (i in fix_values_oneside) {
    theta_warmstart[fix_index] <- i # replace with new fixed value
    pos <- which(fix_values_oneside == i)


    switch(optimizer,
           "fullopt_fixed" = {
             telefit <- tryCatch(telescope(
               theta_init = theta_warmstart,
               FUN_tele = nr_penalty,
               fix_index = c(fix_index, nu_index), # include nu_index
               x1 = x1,
               x2 = x2,
               x3 = x3,
               y = y,
               tau_tele_vec = tau_tele_vec,
               eps_tele_vec = eps_tele_vec,
               list_family = list_family,
               method_c_tilde = method_c_tilde,
               method_c_tilde_deriv = method_c_tilde_deriv,
               kappa_omega = kappa_omega,
               lambda_beta = lambda_beta,
               lambda_alpha = lambda_alpha,
               lambda_nu = lambda_nu,
               initial_step = initial_step,
               max_step_it = max_step_it,
               algorithm = algorithm,
               h_kappa = h_kappa,
               tol = tol,
               max_it = max_it
             ),
             error = function(err) NA
             ) # return NA if error
           },
           "nlm_fixed" = {
             telefit <- tryCatch(telescope_nlm_fixed(
               theta_init = theta_warmstart, # updated estimates
               fix_indices = c(fix_index, nu_index), # index for fixed coef
               x1 = x1,
               x2 = x2,
               x3 = x3,
               y = y,
               tau_tele_vec = tau_tele_vec,
               eps_tele_vec = eps_tele_vec,
               list_family = list_family,
               method_c_tilde = method_c_tilde,
               kappa_omega = kappa_omega,
               lambda_beta = lambda_beta,
               lambda_alpha = lambda_alpha,
               lambda_nu = lambda_nu,
               iterlim = iterlim_nlm,
               stepmax = stepmax_nlm
             ),
             error = function(err) NA
             ) # return NA if error
           }
    )

    if (class(telefit)[1] != "logical") {
      # if not NA, then update initial estimates
      # update initial values
      # careful with index - depends on shape of telescope output
      theta_warmstart <- unlist(telefit[1, 2:(length_theta + 1)]) # large epsilon value
      # can use grep to detect as an alternative
      # pos_extract <- c(grep("beta_", colnames(telefit_init)),
      #                  grep("alpha_", colnames(telefit_init)))
      # theta_warmstart <- unlist(telefit_init[1, pos_extract])
    }
    telefit_list[[pos]] <- telefit
  }
  telefit_list
}

# Select range want to see values ----
fix_coef_optimizer_range_mpr_fix_nu <- function(fix_index,
                                                run_range,
                                                seq_length, # either side (+ 1)
                                                nu_index, # index of nu
                                                nu_fix_value, # fix nu to this value
                                                telemat_opt, # scaled
                                                x1, # data should be scaled
                                                x2,
                                                x3,
                                                y,
                                                kappa_omega,
                                                optimizer,
                                                iterlim_nlm,
                                                stepmax_nlm,
                                                tau_1,
                                                tau_T,
                                                epsilon_1,
                                                epsilon_T,
                                                steps_T,
                                                list_family,
                                                method_c_tilde,
                                                method_c_tilde_deriv,
                                                algorithm,
                                                h_kappa,
                                                lambda_beta,
                                                lambda_alpha,
                                                lambda_nu,
                                                initial_step,
                                                max_step_it,
                                                tol,
                                                max_it) {
  eps_tele_vec <- rev(exp(seq(log(epsilon_T), log(epsilon_1), length = steps_T)))
  tau_tele_vec <- rev(exp(seq(log(tau_T), log(tau_1), length = steps_T)))

  theta_smalleps <- telemat_opt %>% # optimum values at smallest epsilon
    select(contains(c("beta", "alpha", "nu"))) %>%
    slice_tail() %>% # last row = small epsilon
    unlist()

  theta_largeeps <- telemat_opt %>% # optimium values at largest epsilon
    select(contains(c("beta", "alpha", "nu"))) %>%
    slice(1) %>% # first row = large epsilon
    unlist()

  fix_opt <- theta_smalleps[fix_index] # optimium value for fixed coef
  stopifnot(fix_opt > run_range[1] & fix_opt < run_range[2]) # range contains optimium value?

  # * Fix values sequence ----
  seq_values <- c(
    seq(
      from = run_range[1],
      to = run_range[2],
      length = seq_length
    ),
    unname(fix_opt)
  )
  fix_values_lhs <- seq_values[which(seq_values <= fix_opt)]
  fix_values_rhs <- seq_values[which(seq_values > fix_opt)] # don't include opt value

  fix_values_list_init <- list(
    fix_values_lhs,
    fix_values_rhs
  )

  # ** Contain zero? ----
  changesign <- map_dbl(fix_values_list_init, ~ {
    vec <- unique(sign(.x))
    changesign <- ifelse(length(vec) == 2, 1, 0) # true if change sign (length vec == 2), false if not
  })
  fix_values_list <- map2(
    .x = fix_values_list_init,
    .y = changesign,
    ~ {
      if (.y == 1) { # if changesign == true, then include 0
        vec <- c(.x, 0)
        # vec[-which(vec == fix_opt)] # if include 0, remove optimum value (will fit on the other side)
      } else {
        .x # == false, leave vector as is
      }
    }
  ) %>%
    map2(.x = ., .y = c(
      TRUE, # sort lhs by decreasing, rhs increasing
      FALSE
    ), ~ {
      sort(.x, decreasing = .y)
    })

  # * Optimization lhs ----
  telefit_lhs <- fix_values_list[[1]] %>%
    fix_coef_warmstarts_telefits_fix_nu(
      fix_index = fix_index,
      fix_values_oneside = ., # fix_values
      nu_index = nu_index,
      nu_fix_value = nu_fix_value,
      theta_largeeps = theta_largeeps, # large epsilon
      x1 = x1, # data should be scaled
      x2 = x2,
      x3 = x3,
      y = y,
      optimizer = optimizer,
      iterlim_nlm = iterlim_nlm,
      stepmax_nlm = stepmax_nlm,
      tau_tele_vec = tau_tele_vec,
      eps_tele_vec = eps_tele_vec,
      list_family = list_family,
      method_c_tilde = method_c_tilde,
      method_c_tilde_deriv = method_c_tilde_deriv,
      kappa_omega = kappa_omega,
      lambda_beta = lambda_beta,
      lambda_alpha = lambda_alpha,
      lambda_nu = lambda_nu,
      initial_step = initial_step,
      max_step_it = max_step_it,
      algorithm = algorithm,
      h_kappa = h_kappa,
      tol = tol,
      max_it = max_it
    ) # list of length(fix_values) of telemat results

  # * Get theta_large eps from lhs to use as initial for rhs
  theta_largeeps_lhs <- telefit_lhs[[1]] %>% # will be 1st in list
    select(contains(c("beta", "alpha", "nu"))) %>%
    slice(1) %>% # first row = large epsilon
    unlist()

  # Optimization rhs ----
  telefit_rhs <- fix_values_list[[2]] %>%
    fix_coef_warmstarts_telefits_fix_nu(
      fix_index = fix_index,
      fix_values_oneside = ., # fix_values
      nu_index = nu_index,
      nu_fix_value = nu_fix_value,
      theta_largeeps = theta_largeeps_lhs, # large epsilon from left hand side
      x1 = x1, # data should be scaled
      x2 = x2,
      x3 = x3,
      y = y,
      optimizer = optimizer,
      iterlim_nlm = iterlim_nlm,
      stepmax_nlm = stepmax_nlm,
      tau_tele_vec = tau_tele_vec,
      eps_tele_vec = eps_tele_vec,
      list_family = list_family,
      method_c_tilde = method_c_tilde,
      method_c_tilde_deriv = method_c_tilde_deriv,
      kappa_omega = kappa_omega,
      lambda_beta = lambda_beta,
      lambda_alpha = lambda_alpha,
      lambda_nu = lambda_nu,
      initial_step = initial_step,
      max_step_it = max_step_it,
      algorithm = algorithm,
      h_kappa = h_kappa,
      tol = tol,
      max_it = max_it
    ) # list of length(fix_values) of telemat results

  telefit_list <- list(
    telefit_lhs,
    telefit_rhs
  )

  # ** Number of NAs in LHS, RHS ----
  any_na <- telefit_list %>%
    map(~ map_lgl(.x, ~ {
      class(.x)[1] == "logical" # any na?
    })) %>%
    map_dbl(sum)

  seq_length_total <- sum(
    length(telefit_list[[1]]),
    length(telefit_list[[2]])
  )
  # ** Combine to large df ----
  telefit_df <- telefit_list %>%
    flatten() %>%
    list_na_remove() %>%
    # remove any NAs
    map2(
      .x = ., .y = seq(1:((seq_length_total) - sum(any_na))), # seq_length - any_na
      ~ {
        .x %>%
          add_column(
            id = !!.y,
            .before = 1
          )
      }
    ) %>%
    data.table::rbindlist()

  # * Return ----
  list(
    "telefit_df" = telefit_df,
    "any_na_sides" = any_na
  )
}

# get thetas from dataf ---------------------------------------------------
get_thetas_list_df <- function(fit,
                               df) { # all_thetas_dataf
  df %>%
    filter(fit_type == fit) %>%
    select(-sim_id, -fit_type, -sample_kappa_index, -kappa_est) %>%
    group_by(sample_size) %>%
    group_split(.keep = FALSE) %>%
    map(~ {
      .x %>%
        group_by(kappa_index) %>%
        group_split(.keep = FALSE) %>%
        map(~ as_tibble(t(.x)))
    })
}


# Data Specific -----------------------------------------------------------
# information criterion ---------------------------------------------------
likelihood_laplace <- function(theta, # incl s
                               x, # incl int
                               y) {
  n <- length(y)
  p <- ncol(x) - 1

  beta <- theta[1:(p + 1)]
  s <- sqrt(exp(theta[p + 2])) # s^2 = exp(alpha_0)

  mu <- x %*% beta

  mad_vec <- abs(y - mu)

  (-n * log(2 * s)) - ((1 / s) * (sum(mad_vec)))
}

likelihood_laplace_neg <- function(...) -likelihood_laplace(...)
# for BIC impact ---- need to harshly set some coefs to zero
likelihood_laplace_neg_fixed <- function(theta_estimate, # just coefs to be estimaed
                                         fix_indices, # positions of coefs that have been removed already
                                         fix_values, # match fix_indices
                                         x,
                                         y) {
  theta <- theta_estimate
  for (i in fix_indices) {
    pos <- which(fix_indices == i)
    theta <- append(
      x = theta,
      values = fix_values[pos],
      after = i - 1
    ) # insert value after fix_index - 1 position
  }

  likelihood_laplace_neg(theta = theta,
                         x = x,
                         y = y)
}

nlm_fixed_laplace <- function(p, # "initial values should have hte fixed values included here
                              fix_indices,
                              x,
                              y,
                              iterlim_nlm) {
  theta_estimate <- p[-fix_indices] # remove fix_indices value from the vector

  fix_values <- p[fix_indices] # the values that are fixed

  res <- nlm(likelihood_laplace_neg_fixed,
             p = theta_estimate,
             fix_indices = fix_indices,
             fix_values = fix_values,
             x = x,
             y = y,
             iterlim = iterlim_nlm)

  # fill in the fixed values so have complete vector returned
  theta <- res$estimate
  for (i in fix_indices) {
    pos <- which(fix_indices == i)
    theta <- append(
      x = theta,
      values = fix_values[pos],
      after = i - 1
    ) # insert value after fix_index - 1 position
  }

  res$estimate <- theta # replace nlm estimates with complete vector (includes fixed)
  res
}

info_criterion_laplace <- function(theta, # including s
                                   x,
                                   y,
                                   df,
                                   lambda) {
  n <- nrow(x)
  lambda <- eval(parse(text = lambda))

  (-2 * likelihood_laplace(
    theta = theta,
    x = x,
    y = y
  )) + (lambda * (df + 2)) # + 2 for estimate of intercept and b parameter
}

likelihood_normal <- function(theta,
                              x,
                              y) {
  n <- nrow(x)
  p <- ncol(x) - 1

  beta <- theta[1:(p + 1)]
  sigma_sq <- exp(theta[(p + 2)]) # convert to sigma (alpha_0 = log(sigma_sq))

  mu <- x %*% beta

  like <- -(n / 2) * log(2 * pi) - ((n / 2) * log(sigma_sq)) - ((1 / (2 * (sigma_sq))) * sum((y - mu)^2))
  like # return likelihood value
}

likelihood_normal_neg <- function(...) -likelihood_normal(...)

likelihood_normal_neg_fixed <- function(theta_estimate, # just coefs to be estimaed
                                        fix_indices, # positions of coefs that have been removed already
                                        fix_values, # match fix_indices
                                        x,
                                        y) {
  theta <- theta_estimate
  for (i in fix_indices) {
    pos <- which(fix_indices == i)
    theta <- append(
      x = theta,
      values = fix_values[pos],
      after = i - 1
    ) # insert value after fix_index - 1 position
  }

  likelihood_normal_neg(theta = theta,
                        x = x,
                        y = y)
}

nlm_fixed_normal <- function(p, # "initial values should have hte fixed values included here
                             fix_indices,
                             x,
                             y,
                             iterlim_nlm) {
  theta_estimate <- p[-fix_indices] # remove fix_indices value from the vector

  fix_values <- p[fix_indices] # the values that are fixed

  res <- nlm(likelihood_normal_neg_fixed,
             p = theta_estimate,
             fix_indices = fix_indices,
             fix_values = fix_values,
             x = x,
             y = y,
             iterlim = iterlim_nlm)

  # fill in the fixed values so have complete vector returned
  theta <- res$estimate
  for (i in fix_indices) {
    pos <- which(fix_indices == i)
    theta <- append(
      x = theta,
      values = fix_values[pos],
      after = i - 1
    ) # insert value after fix_index - 1 position
  }

  res$estimate <- theta # replace nlm estimates with complete vector (includes fixed)
  res
}

info_criterion_normal <- function(theta,
                                  x,
                                  y,
                                  df,
                                  lambda) {
  n <- nrow(x)
  p <- ncol(x) - 1
  lambda <- eval(parse(text = lambda))

  beta <- theta[1:(p + 1)]
  sigma_sq <- exp(theta[(p + 2)]) # convert to sigma (alpha_0 = log(sigma_sq))

  mu <- x %*% beta

  like <- -(n / 2) * log(2 * pi) - ((n / 2) * log(sigma_sq)) - ((1 / (2 * (sigma_sq))) * sum((y - mu)^2))

  (-2 * like) + (lambda * (df + 2)) # info criterion
}

info_criterion_sgnd <- function(theta, # smooth gnd
                                x1,
                                x2,
                                x3,
                                y,
                                tau,
                                list_general,
                                method_c_tilde,
                                kappa_omega,
                                df,
                                lambda) {
  n <- nrow(x1)
  lambda <- eval(parse(text = lambda))

  (-2 * likelihood(theta = theta,
                   x1 = x1,
                   x2 = x2,
                   x3 = x3,
                   y = y,
                   tau = tau,
                   list_general = list_general,
                   method_c_tilde = method_c_tilde,
                   kappa_omega = kappa_omega)) + (lambda * (df + 3)) # + 3 for estimate of intercepts (beta, alpha, nu)
}

# boostrap kurtosis -------------------------------------------------------
bootstrap_kurtosis <- function(x_raw_incl_int, # error distrib after fitting lm
                               y_raw,
                               n_boots) {
  n_all <- nrow(x_raw_incl_int)
  p <- ncol(x_raw_incl_int) - 1

  indices_list <- replicate(
    n = n_boots,
    expr = sample(1:n_all, size = n_all, rep = TRUE),
    simplify = FALSE
  )

  excess_kurt_vec <- NULL
  for (i in 1:n_boots) {
    x_now <- x_raw_incl_int[indices_list[[i]], ]
    y_now <- y_raw[indices_list[[i]]]

    ## Scale X ----
    x_scale <- scale(x_now[, -1],
                     center = FALSE,
                     scale = apply(x_now[, -1], 2, sd)
    )
    x_scale <- cbind(rep(1, n_all), x_scale) # column of 1's for intercept # n_all same as n_now
    x_sd <- apply(x_now[, -1], 2, sd) # save sd for later to transform back

    lm_fit <- lm(y_now ~ x_scale[, -1])

    coef_fit <- coef(lm_fit) / c(1, x_sd) # rescale

    pred <- x_now %*% coef_fit

    resid <- y_now - pred

    excess_kurt <- moments::kurtosis(resid) - 3 # subtract 3 for normal dist
    excess_kurt_vec[i] <- excess_kurt
  }
  excess_kurt_vec
}

# ladlasso ic -------------------------------------------------------------
flare_ladlasso_ic_selection <- function(x_raw_int,
                                        y_raw,
                                        lambda_ic,
                                        lambda_vec) {
  n <- length(y_raw)
  lambda_ic <- eval(parse(text = lambda_ic))

  fit_ic <- flare::slim(
    X = x_raw_int[, -1],
    Y = y_raw,
    lambda = lambda_vec,
    method = "lq",
    q = 1
  )

  # extract coef
  lambda_vec <- fit_ic$lambda
  beta_0_vec <- fit_ic$intercept

  names_coef_labels <- colnames(x_raw_int)[-1]

  coef_mat_prep <- rbind(
    beta_0_vec,
    fit_ic$beta
  ) %>%
    t() %>%
    as_tibble() %>%
    rename_all(~ c("intercept", names_coef_labels))

  # df_vec <- coef_mat_prep[, -1] %>%
  #   apply(., 1, function(x) sum(x != 0)) # count number of nonzeros
  df_vec <- fit_ic$df # extract df from function

  dataf_ic <- coef_mat_prep %>%
    add_column(
      lambda = lambda_vec,
      df = df_vec,
      .before = 1
    )

  # b estimate (scale parameter)
  list_all_coef <- coef_mat_prep %>%
    split(seq(nrow(coef_mat_prep))) %>%
    map(~ unlist(.x)) # list of coefficients

  dataf_pred <- list_all_coef %>%
    map_dfc(~ (x_raw_int %*% .x))
  colnames(dataf_pred) <- paste0("lambda", 1:ncol(dataf_pred))

  # sigma est = s est
  s_est <- dataf_pred %>%
    map_dbl(~ 1.4826 * (median(abs(y_raw - .x)))) # MAD

  ## add to dataf
  dataf_ic <- dataf_ic %>%
    add_column(s = s_est) %>%
    mutate(log_ssquare = log(s^2), # alpha_0 = log(s^2)
           .keep = "unused")

  ## Calculate IC ----
  list_df_theta <- dataf_ic %>%
    select(-lambda) %>%
    split(seq(nrow(dataf_ic))) %>%
    map(~ unlist(.x))

  vec_ic <- list_df_theta %>%
    map_dbl(~ (info_criterion_laplace(
      theta = .x[-1], # includes alpha_0
      x = x_raw_int,
      y = y_raw,
      df = .x[1],
      lambda = lambda_ic
    )))
  dataf_ic <- dataf_ic %>%
    add_column(ic = vec_ic) %>%
    rowid_to_column()
  dataf_ic
}

flare_ladlasso_cv_base <- function(x_raw_int,
                                   y_raw,
                                   indices,
                                   lambda_val) { # for single lambda
  test_error <- numeric()
  for (i in 1:max(indices)) {
    # Split test and train
    x_raw_train <- x_raw_int[indices != i, ]
    x_raw_test <- x_raw_int[indices == i, ]

    y_raw_train <- y_raw[indices != i]
    y_raw_test <- y_raw[indices == i]

    # Fit ladlasso
    fit_flare <- flare::slim(
      X = x_raw_train[, -1],
      Y = y_raw_train,
      lambda = lambda_val, # single lambda value
      method = "lq",
      q = 1
    )

    coef_flare <- c(fit_flare$intercept, fit_flare$beta)

    pred <- x_raw_test %*% coef_flare

    test_error[i] <- mean(abs(y_raw_test - pred)) # Mean absolute deviation
  }
  mean(test_error)
}

flare_ladlasso_cv_selection <- function(x_raw_int,
                                        y_raw,
                                        indices,
                                        lambda_vec) {
  length_lambda_vec <- length(lambda_vec)
  test_error <- numeric()

  for (i in lambda_vec) {
    pos <- which(lambda_vec == i)

    cv_res <- flare_ladlasso_cv_base(
      x_raw_int = x_raw_int,
      y_raw = y_raw,
      indices = indices,
      lambda_val = i
    ) # different lambda values
    test_error[pos] <- cv_res # store mean error over folds
    print((pos / length_lambda_vec) * 100)
  }
  tibble(
    lambda = lambda_vec,
    cv_error = test_error
  )
}

s_est_func <- function(x_raw_int,
                       y_raw,
                       theta) {
  pred <- x_raw_int %*% theta
  s_est <- 1.4826 * median(abs(y_raw - pred))
  s_est
}

get_see_now <- function(info_list) {
  inv_pen <- solve(info_list$observed_information_penalized)

  sqrt(diag(inv_pen %*% info_list$observed_information_unpenalized %*% inv_pen))
}

# normal adaptive lasso ----
lasso_alasso_scale_ic_selection <- function(x_scale, # include col of 1s
                                            y,
                                            lambda_ic,
                                            ...) {
  n <- length(y)
  lambda_ic <- eval(parse(text = lambda_ic))

  # glmnet fit ----------
  fit_ic <- glmnet::glmnet(
    x = x_scale[, -1],
    y = y,
    alpha = 1,
    standardize = FALSE,
    ...
  ) # ... for penalty.factor

  # Extract
  coef_mat_scale <- as_tibble(t(as.matrix(coef(fit_ic)))) %>%
    rename("X0" = "(Intercept)")
  # Dataframe
  dataf_ic <- coef_mat_scale %>%
    add_column(
      lambda = fit_ic$lambda,
      df = fit_ic$df,
      .before = "X0"
    )
  # Sigma Estimate --------------------
  ## Predictions
  list_all_coef <- coef_mat_scale %>%
    split(seq(nrow(coef_mat_scale))) %>%
    map(~ unlist(.x)) # list of coefficients

  dataf_pred <- list_all_coef %>%
    map_dfc(~ (x_scale %*% .x))
  colnames(dataf_pred) <- paste0("lambda", 1:ncol(dataf_pred))

  ## (y - yhat)^2
  dataf_y_yhatsq <- dataf_pred %>%
    map_dfc(~ ((y - .x)^2))

  ## rss
  rss <- dataf_y_yhatsq %>%
    map_dbl(~ sum(.x))
  ## sigma
  sigma_est <- tibble(rss,
                      df = dataf_ic$df
  ) %>%
    mutate(sigma = (sqrt(rss / (n - df - 1)))) %>%
    pull(sigma)
  ## add to dataf
  dataf_ic <- dataf_ic %>%
    add_column(
      sigma = sigma_est,
      logsigmasq = log(sigma^2)
    )

  # Calculate IC ---------------------
  ## rearrange df and coef into list
  list_df_theta <- dataf_ic %>% # df, and coef
    select(-lambda, -sigma) %>%
    split(seq(nrow(dataf_ic))) %>%
    map(~ unlist(.x))

  # calculate IC
  vec_ic <- list_df_theta %>%
    map_dbl(~ (info_criterion_normal(
      theta = .x[-1],
      x = x_scale,
      y = y,
      df = .x[1],
      lambda = lambda_ic # BIC or AIC
    )))
  dataf_ic_scale <- dataf_ic %>%
    add_column(ic = vec_ic) %>%
    rowid_to_column() # scaled coef

  # Output: IC Dataframe -------------
  return(dataf_ic_scale)
}

# ** LASSO ALASSO COEF FINDER CV ------------------------------------------
coef_sigma_scale_cvglmnet <- function(fit_scale, # extract standardized coef and calculate sigma
                                      x_scale,
                                      y,
                                      lambda_choice) { # "lambda.1se", "lambda.min"
  n <- length(y)
  beta_scale <- as.vector(coef(fit_scale, s = lambda_choice)) # standardized coef
  df <- sum(beta_scale[-1] != 0)

  # Sigma Estimate
  yhat <- as.vector(x_scale %*% beta_scale)
  rss <- sum((y - yhat)^2)
  sigmasq <- (rss / (n - df - 1))
  # Coef
  return(theta_scale <- c(beta_scale, log(sigmasq)))
}

# ** fitting functions for data comparison --------------------------------
fit_models <- function(type, # either "normal" for lm and adaptive lasso, or "robust" for rlm and ladlasso
                       x_raw_incl_int,
                       y_raw,
                       optimizer,
                       iterlim_nlm,
                       stepmax_nlm,
                       tau_1_mpr,
                       tau_T_mpr,
                       tau_1_spr,
                       tau_T_spr,
                       epsilon_1,
                       epsilon_T,
                       steps_T,
                       list_family,
                       method_c_tilde,
                       method_c_tilde_deriv,
                       algorithm,
                       h_kappa,
                       kappa_omega,
                       lambda_beta,
                       lambda_alpha,
                       lambda_nu,
                       initial_step,
                       max_step_it,
                       tol,
                       max_it,
                       lambda_vec,
                       indices, # for ladlasso
                       formula_bamlss) {
  x_raw_int <- x_raw_incl_int
  stopifnot(all(x_raw_int[, 1] == 1)) # ensure column for intercept included

  n_raw <- length(y_raw)

  p <- ncol(x_raw_int) - 1 # number of variables, exclude intercept
  names_coef_labels <- colnames(x_raw_int)[-1]

  ## Scale X --
  x_scale <- x_raw_int[, -1] %>%
    as_tibble() %>%
    rename_with(~ paste0("X", 1:p)) %>%
    scale(
      center = FALSE,
      scale = apply(., 2, sd)
    ) %>%
    as_tibble() %>%
    add_column(
      X0 = 1,
      .before = "X1"
    ) %>%
    as.matrix()

  x_sd <- apply(x_raw_int[, -1], 2, sd) # save sd for later
  x_sd_theta <- c(1, x_sd, 1, x_sd, 1)
  x_sd_theta_spr <- c(1, x_sd, 1, 1)

  x3 <- as.matrix(rep(1, n_raw))

  dataset_scale <- as.data.frame(x_scale[, -1]) %>%
    add_column(y_raw) %>%
    rename_all(~ c(names_coef_labels, "y"))

  ## Fits --
  theta_init <- initial_values_func(
    x1 = x_scale,
    x2 = x_scale,
    x3 = x3, # col of 1s
    y = y_raw,
    method_initial_values = "lm",
    other_initial_values_vec = NA,
    kappa_omega = kappa_omega,
    fix_nu_value_init = NA,
    tau = tau_1_mpr,
    epsilon_1 = epsilon_1,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  )
  theta_init_spr <- c(theta_init[1:(p + 2)], theta_init[(2 * p) + 3]) # incl alpha_0 and nu_0

  # mpr -
  time_mpr <- system.time(fit_mpr <- fitting_func_base(
    x1 = x_scale,
    x2 = x_scale,
    x3 = x3,
    y = y_raw,
    optimizer = optimizer,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    nu_profile_vec = NA,
    theta_init = theta_init,
    tau_1 = tau_1_mpr,
    tau_T = tau_T_mpr,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  ))
  extract_mpr <- extract_theta_plike_val(
    fit_mpr,
    optimizer
  )
  theta_mpr_scale <- extract_mpr$theta
  theta_mpr <- as.vector(theta_mpr_scale / x_sd_theta)

  ## spr -
  time_spr <- system.time(fit_spr <- fitting_func_base(
    x1 = x_scale,
    x2 = x3, # col of 1s
    x3 = x3,
    y = y_raw,
    optimizer = optimizer,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    nu_profile_vec = NA,
    theta_init = theta_init_spr,
    tau_1 = tau_1_spr,
    tau_T = tau_T_spr,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it
  ))
  extract_spr <- extract_theta_plike_val(
    fit_spr,
    optimizer
  )
  theta_spr_scale_prep <- extract_spr$theta # include 0s for alphas
  theta_spr_scale <- c(theta_spr_scale_prep[1:(p + 2)], rep(0, p), theta_spr_scale_prep[length(theta_spr_scale_prep)])
  theta_spr <- as.vector(theta_spr_scale / x_sd_theta)

  ## bamlss -
  time_bamlss <- system.time(fit_bamlss <- bamlss(formula_bamlss,
                                                  data = dataset_scale, # scaled dataset
                                                  family = PE,
                                                  optimizer = opt_lasso_MON, # no penalty on the shape
                                                  criterion = "BIC",
                                                  multiple = TRUE))

  theta_bamlss_scale <- extract_coef_bamlss(fit_bamlss,
                                            names_coef_labels)$value

  theta_bamlss <- as.vector(theta_bamlss_scale / x_sd_theta)

  # Fit either normal or robust functions ----
  switch (type,
          "robust" = {
            ## ladlassoic -
            time_ladlassoic <- system.time(fit_ladlassoic <- flare_ladlasso_ic_selection(
              x_raw_int = x_scale, # scaled data
              y_raw = y_raw,
              lambda_ic = "log(n)",
              lambda_vec = lambda_vec
            ))
            opt_fit_ladlassoic <- fit_ladlassoic %>%
              slice_min(ic)

            theta_ladlassoic_scale_prep <- opt_fit_ladlassoic %>%
              select(intercept, contains("X"), log_ssquare) %>%
              unlist()

            theta_ladlassoic_scale <- c(theta_ladlassoic_scale_prep, rep(0, p), 0) # include 0's for alphas and kappa
            theta_ladlassoic <- as.vector(theta_ladlassoic_scale / x_sd_theta)

            ## rlm -
            time_rlm <- system.time(fit_rlm <- MASS::rlm(y_raw ~ x_scale[, -1]))
            theta_rlm_scale_prep <- c(coef(fit_rlm), log(summary(fit_rlm)$sigma^2))
            theta_rlm_scale <- c(theta_rlm_scale_prep, rep(0, p), 0)
            theta_rlm <- as.vector(theta_rlm_scale / x_sd_theta)

            out <- list(
              "data_info" = list(
                "x_raw_incl_int" = x_raw_incl_int,
                "y_raw" = y_raw,
                "x_sd" = x_sd
              ),
              "fits" = list(
                "tele_mpr" = fit_mpr,
                "tele_spr" = fit_spr,
                "fit_ladlassoic" = fit_ladlassoic,
                "fit_rlm" = fit_rlm,
                "fit_bamlss" = fit_bamlss
              ),
              "thetas" = list(
                "theta_mpr" = theta_mpr,
                "theta_spr" = theta_spr,
                "theta_ladlassoic" = unname(theta_ladlassoic),
                "theta_rlm" = unname(theta_rlm),
                "theta_bamlss" = theta_bamlss
              ),
              "thetas_scale" = list(
                "theta_mpr_scale" = theta_mpr_scale,
                "theta_spr_scale" = theta_spr_scale,
                "theta_ladlassoic_scale" = unname(theta_ladlassoic_scale),
                "theta_rlm_scale" = unname(theta_rlm_scale),
                "theta_bamlss_scale" = theta_bamlss_scale
              ),
              "time_fits" = list(
                "time_mpr" = time_mpr,
                "time_spr" = time_spr,
                "time_ladlassoic" = time_ladlassoic,
                "time_rlm" = time_rlm,
                "time_bamlss" = time_bamlss
              )
            )
          },
          # normal functions ----
          "normal" = {

            time_lm <- system.time(lm_fit <- lm(y_raw ~ x_scale[, -1]))
            lm_coef <- coef(lm_fit)

            theta_lm_scale_prep <- c(coef(lm_fit), log((summary(lm_fit)$sigma)^2))
            theta_lm_scale <- c(theta_lm_scale_prep, rep(0, p), 0)
            theta_lm <- as.vector(theta_lm_scale / x_sd_theta)

            # adaptive lasso
            param_weight_alasso <- 1 / abs(lm_coef[-1])

            time_alassoic <- system.time(fit_alassoic_scale <- lasso_alasso_scale_ic_selection(
              x_scale = x_scale,
              y = y_raw,
              lambda_ic = "log(n)", # BIC
              penalty.factor = param_weight_alasso
            ))
            theta_alassoic_scale_prep <- fit_alassoic_scale %>%
              filter(ic == min(ic)) %>%
              select(-rowid, -lambda, -df, -sigma, -ic) %>%
              unlist() %>%
              as.numeric()
            theta_alassoic_scale <- c(theta_alassoic_scale_prep, rep(0, p), 0)
            theta_alassoic <- as.vector(theta_alassoic_scale / x_sd_theta)

            out <- list(
              "data_info" = list(
                "x_raw_incl_int" = x_raw_incl_int,
                "y_raw" = y_raw,
                "x_sd" = x_sd
              ),
              "fits" = list(
                "tele_mpr" = fit_mpr,
                "tele_spr" = fit_spr,
                "fit_alassoic" = fit_alassoic_scale,
                "fit_lm" = lm_fit,
                "fit_bamlss" = fit_bamlss
              ),
              "thetas" = list(
                "theta_mpr" = theta_mpr,
                "theta_spr" = theta_spr,
                "theta_alassoic" = unname(theta_alassoic),
                "theta_lm" = unname(theta_lm),
                "theta_bamlss" = theta_bamlss
              ),
              "thetas_scale" = list(
                "theta_mpr_scale" = theta_mpr_scale,
                "theta_spr_scale" = theta_spr_scale,
                "theta_alassoic_scale" = unname(theta_alassoic_scale),
                "theta_lm_scale" = unname(theta_lm_scale),
                "theta_bamlss_scale" = theta_bamlss_scale
              ),
              "time_fits" = list(
                "time_mpr" = time_mpr,
                "time_spr" = time_spr,
                "time_alassoic" = time_alassoic,
                "time_lm" = time_lm,
                "time_bamlss" = time_bamlss
              )
            )
          }
  )
  out
}

get_thetas_vec <- function(fit,
                           df,
                           model_type) { # for single dataframe

  switch (model_type,
          "mpr" = {
            df %>%
              filter(fit_type == fit) %>%
              select(-fit_type) %>%
              unlist()
          },
          "spr" = {
            df %>%
              filter(fit_type == fit) %>%
              select(-fit_type) %>%
              select(contains(c("beta", "alpha_0", "nu_0"))) %>%
              unlist()
          },
          "no_nu" = {
            df %>%
              filter(fit_type == fit) %>%
              select(-fit_type) %>%
              select(contains(c("beta", "alpha_0"))) %>%
              unlist()
          }
  )
}


# ** out of sample cv -----------------------------------------------------
out_of_sample_metrics <- function(index_oos,
                                  type,
                                  x_raw_incl_int,
                                  y_raw,
                                  indices_oos,
                                  less_than_value,
                                  less_than_value_bamlss,
                                  optimizer,
                                  iterlim_nlm,
                                  stepmax_nlm,
                                  tau_1_mpr,
                                  tau_T_mpr,
                                  tau_1_spr,
                                  tau_T_spr,
                                  epsilon_1,
                                  epsilon_T,
                                  steps_T,
                                  list_family,
                                  method_c_tilde,
                                  method_c_tilde_deriv,
                                  algorithm,
                                  h_kappa,
                                  kappa_omega,
                                  lambda_beta,
                                  lambda_alpha,
                                  lambda_nu,
                                  initial_step,
                                  max_step_it,
                                  tol,
                                  max_it,
                                  lower_cutoff,
                                  upper_cutoff,
                                  formula_bamlss) {
  x_raw_int <- x_raw_incl_int
  # Split test and train --
  x_raw_train <- x_raw_int[indices_oos != index_oos, ]
  x_raw_test <- x_raw_int[indices_oos == index_oos, ]

  y_raw_train <- y_raw[indices_oos != index_oos]
  y_raw_test <- y_raw[indices_oos == index_oos]

  # Fit models --
  indices <- sample(1:10, size = nrow(x_raw_train), rep = TRUE) # indices for training
  p <- ncol(x_raw_train) - 1
  lambda_vec <- NA # will overwrite if type == "robust"
  # Scale training data to get lambda vec from glmnet
  if (type == "robust") { # not required if type == "normal"
    x_scale_t <- x_raw_train[, -1] %>%
      as_tibble() %>%
      rename_with(~ paste0("X", 1:p)) %>%
      scale(
        center = FALSE,
        scale = apply(., 2, sd)
      ) %>%
      as_tibble() %>%
      add_column(
        X0 = 1,
        .before = "X1"
      ) %>%
      as.matrix()
    lambda_vec <- glmnet::glmnet(
      x = x_scale_t[, -1],
      y = y_raw_train,
      alpha = 1,
      standardize = FALSE
    )$lambda # get lambda vec for ladlasso
  }

  fitres <- fit_models(
    type = type,
    x_raw_incl_int = x_raw_train,
    y_raw = y_raw_train,
    optimizer = optimizer,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    tau_1_mpr = tau_1_mpr,
    tau_T_mpr = tau_T_mpr,
    tau_1_spr = tau_1_spr,
    tau_T_spr = tau_T_spr,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it,
    lambda_vec = lambda_vec,
    indices = indices,
    formula_bamlss = formula_bamlss
  )

  x_sd_train <- fitres$data_info$x_sd
  x_sd_train_spr <- c(1, x_sd_train, 1, 1) # no alphas
  x_sd_train_theta <- c(1, x_sd_train, 1, x_sd_train, 1)

  x_scale_train <- as_tibble(x_raw_train) %>%
    map2_dfc(.x = .,
             .y = c(1, x_sd_train), ~ {
               .x / .y # scale the training data (required for some BIC calculations later)
             }) %>%
    as.matrix()

  x3_train <- as.matrix(rep(1, nrow(x_raw_train)))

  names_coef_no_nu <- c(
    paste0("beta_", 0:p),
    paste0("alpha_", 0:p)
  )
  names_coef <- c(names_coef_no_nu, "nu_0")

  names_fits <- names(fitres$thetas) %>%
    str_remove("theta_")

  # get estimates --
  thetas_dataf <- fitres$thetas %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    rename_all(~names_coef) %>%
    add_column(
      fit_type = factor(names_fits, levels = names_fits),
      .before = "beta_0"
    )

  thetas_dataf_t <- thetas_dataf %>%
    select(-fit_type) %>%
    t() %>%
    as_tibble()

  betas_dataf_t <- thetas_dataf %>% # transpose
    select(
      -contains(c("alpha", "nu")), # remove alpha
      -fit_type
    ) %>%
    t() %>%
    as_tibble() # (p + 1) * n_fits (each column = 1 fit)

  # for mpr, spr and bamlss, extract scaled coef and then see what are zero coef
  thetas_dataf_scale <- fitres$thetas_scale %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    rename_all(~names_coef) %>%
    add_column(
      fit_type = factor(names_fits, levels = names_fits),
      .before = "beta_0"
    )

  mpr_spr_dataf_round <- thetas_dataf_scale %>%
    filter(str_detect(fit_type, "mpr|spr")) %>%
    mutate_at(vars(contains(c("beta", "alpha"))), ~ case_when(
      abs(.) <= less_than_value ~ 0,
      TRUE ~ .
    )) %>%
    group_by(fit_type) %>%
    group_split(.keep = FALSE) %>%
    map(~ {
      (unlist(.x) / x_sd_train_theta) %>% # unscale
        as_tibble_row()
    }) %>%
    map2(.x = ., .y = c("mpr", "spr"), ~ {
      .x %>%
        add_column(fit_type = .y,
                   .before = 1)
    }) %>%
    data.table::rbindlist() %>%
    as_tibble()

  mpr_spr_dataf_round <- thetas_dataf_scale %>%
    filter(str_detect(fit_type, "mpr|spr")) %>%
    mutate_at(vars(contains(c("beta", "alpha"))), ~ case_when(
      abs(.) <= less_than_value ~ 0,
      TRUE ~ .
    )) %>%
    group_by(fit_type) %>%
    group_split(.keep = FALSE) %>%
    map(~ {
      (unlist(.x) / x_sd_train_theta) %>%
        as_tibble_row()
    }) %>%
    map2(.x = ., .y = c("mpr", "spr"), ~ {
      .x %>%
        add_column(fit_type = .y,
                   .before = 1)
    }) %>%
    data.table::rbindlist() %>%
    as_tibble()

  thetas_dataf_round <- bind_rows(
    mpr_spr_dataf_round,
    (thetas_dataf %>%
       filter(!str_detect(fit_type, "mpr|spr")))
  ) %>%
    mutate(fit_type = factor(fit_type, levels = names_fits))

  nu_est_bamlss <- exp(fitres$thetas$theta_bamlss[(2*p) + 3]) # loglink func
  alpha_int <- p + 2
  int_remove <- c(1, alpha_int, (2 * p) + 3) # positions of intercepts to remove

  # Unseen test data: predictions --
  predictions_list <- betas_dataf_t %>%
    map(~ {
      as.vector((x_raw_test %*% .x)) # X %*% beta
    }) # list of n_fit predictions as vector

  # * MSE --
  y_yhatsq_list <- predictions_list %>%
    map(~ {
      ((y_raw_test - .x)^2) # (true - pred)^2
    }) # list of n_fit (y - y_hat)^2

  mse_vec <- y_yhatsq_list %>%
    map_dbl(mean) # vector of mse

  # * MAD --
  y_yhatabs_list <- predictions_list %>%
    map(~ {
      (abs(y_raw_test - .x)) # abs(true - pred)
    }) # list of n_fit abs(y - y_hat)

  mad_vec <- y_yhatabs_list %>%
    map_dbl(mean) # vector of mad

  # * MAPE --
  y_yhatmape_list <- predictions_list %>%
    map(~ {
      (abs((y_raw_test - .x) / y_raw_test)) # mean<U+FFFD>(true - predict) / true<U+FFFD>
    }) # list of n_fit

  mape_vec <- y_yhatmape_list %>%
    map_dbl(mean) # vector of mape

  # * Prediction Intervals ----
  nominal <- 0.95
  sig <- 1 - nominal
  z_alpha_value <- nominal + (0.5 * sig) # 0.975

  zk_vals <- thetas_dataf %>%
    filter(fit_type %in% c("mpr", "spr")) %>%
    pull(nu_0) %>%
    nu_to_kappa(., kappa_omega) %>%
    map_dbl(~ {
      qgndnorm(p = z_alpha_value,
               mu = 0,
               s = 1,
               kappa = .x, # different kappa values
               tau = tau_T_mpr)
    })
  stopifnot(tau_T_mpr == tau_T_spr)
  switch (type,
          "robust" = {
            zother_vals <- rmutil::qlaplace(z_alpha_value, m = 0, s = 1) # for laplace dist
          },
          "normal" = {
            zother_vals <- qnorm(z_alpha_value) # normal distribution
          }
  )

  zbamlss <- qPE(z_alpha_value, # this actually doesn't matter - will replace it with upper and lower intervals from the function
                 mu = 0,
                 sigma = 1,
                 nu = nu_est_bamlss) # take nu estimate from training
  zk_vals_all <- c(zk_vals, rep(zother_vals, 2), # for lasso type and lm type
                   zbamlss)

  # sigma estimate on test data --
  # (only mpr will depend on test data, other fits have constant (alpha0) * 1, remaining alpha = 0)
  mpr_spr_s_list <- thetas_dataf %>%
    filter(str_detect(fit_type, "mpr|spr")) %>%
    group_by(fit_type) %>%
    select(
      fit_type,
      contains("alpha")
    ) %>%
    group_split(.keep = FALSE) %>%
    map(~ {
      as.numeric(sqrt(exp(x_raw_test %*% unlist(.x)))) # s estimate = sqrt(phi)
    })

  # now sigma = s (not sqrt(variance))
  other_sigma_list <- thetas_dataf %>%
    filter(!str_detect(fit_type, "mpr|spr|bamlss")) %>%
    pull(alpha_0) %>%
    map(~ {

      switch (type,
              "robust" = {
                s2 <- exp(.x) # s^2 = exp(alpha_0)
                s <- sqrt(s2) # s = scale parameter

                out <- rep(s, nrow(x_raw_test)) # return vector of s
              },
              "normal" = {
                s2 <- exp(.x)
                sigma <- sqrt(s2) # variance normal = s^2 = sigma^2
                out <- rep(sigma, nrow(x_raw_test)) # return vector of sigma
              }
      )
      out

    })

  bamlss_sigma_list <- list(as.numeric(exp(x_raw_test %*% (thetas_dataf %>%
                                                             filter(fit_type == "bamlss") %>%
                                                             select(contains("alpha")) %>%
                                                             unlist()))))

  s_list <- append(
    mpr_spr_s_list,
    other_sigma_list
  ) %>%
    append(.,
           bamlss_sigma_list)

  s_width_dataf <- s_list %>%
    map(~ {
      .x %>%
        as_tibble() %>%
        t() %>%
        as_tibble()
    }) %>%
    do.call(rbind, .) %>%
    rename_all(~ c(paste0("s_n", 1:nrow(x_raw_test)))) %>%
    add_column(
      fit_type = factor(names_fits, levels = names_fits),
      .before = 1
    )

  # combine to list of tibbles
  ytest_yhat_s_list <- map2(
    .x = predictions_list,
    .y = s_list, ~ {
      tibble(
        ytest = y_raw_test,
        yhat = .x,
        s = .y
      )
    }
  ) # list of nfit, each with tibble

  # perform operations
  prediction_intervals_list <- map2(.x = ytest_yhat_s_list,
                                    .y = zk_vals_all, ~ {

                                      z_val <- .y # different z_val for each fit_type
                                      .x %>%
                                        mutate(
                                          P_low = (yhat - (z_val * s)), # lower limit
                                          P_upp = (yhat + (z_val * s)), # upper limit
                                          P = ((ytest > P_low) & (ytest < P_upp)) # lie in limit?
                                        )
                                    })

  # prop true
  mtrue_p_vec <- prediction_intervals_list %>%
    map_dbl(~ {
      mean(.x$P)
    })

  names_prediction_metrics <- c(
    "mse",
    "mad",
    "mape",
    "mtrue_p"
  )

  names_prediction_metrics_vec <- paste0(names_prediction_metrics, "_vec")

  prediction_metrics_dataf <- map(names_prediction_metrics_vec, ~ {
    get(.x) %>%
      setNames(names_fits)
  }) %>%
    setNames(names_prediction_metrics) %>%
    purrr::transpose() %>%
    map(unlist) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    add_column(
      fit_type = names_fits,
      .before = "mse"
    )

  # * BIC -- on the training data
  mpr_bic_val <- -2 * (fitres$fits$tele_mpr %>%
                         extract_theta_plike_val(., "nlm") %>%
                         .$"plike_val") # -2 * plike_val
  spr_bic_val <- -2 * (fitres$fits$tele_spr %>%
                         extract_theta_plike_val(., "nlm") %>%
                         .$"plike_val") # -2 * plike_val

  # * bamlss --
  bamlss_bic_val <- BIC(fitres$fits$fit_bamlss)

  switch (type,
          "robust" = {
            # ladlassoic -- estimate alpha_0 and kappa
            ladlassoic_sgnd_est_scale <- estimate_scale_kappa_func(theta_spr_scale = get_thetas_vec(fit = "ladlassoic",
                                                                                                    df = thetas_dataf_scale,
                                                                                                    model_type = "no_nu"),
                                                                   nu_0 = kappa_to_nu(1, kappa_omega), # initial value for opt
                                                                   x_scale = x_scale_train,
                                                                   y = y_raw_train,
                                                                   tau = tau_T_mpr,
                                                                   list_general = list_family$general,
                                                                   method_c_tilde = method_c_tilde,
                                                                   kappa_omega = kappa_omega,
                                                                   iterlim_nlm = iterlim_nlm,
                                                                   stepmax = NA)
            ladlassoic_bic_val <- info_criterion_sgnd(theta = ladlassoic_sgnd_est_scale / x_sd_train_spr, # unscale
                                                      x1 = x_raw_train,
                                                      x2 = x3_train,
                                                      x3 = x3_train,
                                                      y = y_raw_train,
                                                      tau = tau_T_mpr,
                                                      list_general = list_family$general,
                                                      method_c_tilde = method_c_tilde,
                                                      kappa_omega = kappa_omega,
                                                      df = (fitres$fits$fit_ladlassoic %>% slice_min(ic) %>% pull(df)),
                                                      lambda = "log(n)")
            # rlm -- estimate alpha_0 and kappa
            rlm_sgnd_est_scale <- estimate_scale_kappa_func(theta_spr_scale = get_thetas_vec(fit = "rlm",
                                                                                             df = thetas_dataf_scale,
                                                                                             model_type = "no_nu"),
                                                            nu_0 = kappa_to_nu(1, kappa_omega), # initial value for opt
                                                            x_scale = x_scale_train,
                                                            y = y_raw_train,
                                                            tau = tau_T_mpr,
                                                            list_general = list_family$general,
                                                            method_c_tilde = method_c_tilde,
                                                            kappa_omega = kappa_omega,
                                                            iterlim_nlm = iterlim_nlm,
                                                            stepmax = NA)

            rlm_bic_val <- info_criterion_sgnd(theta = rlm_sgnd_est_scale / x_sd_train_spr, # unscale
                                               x1 = x_raw_train,
                                               x2 = x3_train,
                                               x3 = x3_train,
                                               y = y_raw_train,
                                               tau = tau_T_mpr,
                                               list_general = list_family$general,
                                               method_c_tilde = method_c_tilde,
                                               kappa_omega = kappa_omega,
                                               df = p, # no penalization
                                               lambda = "log(n)")

            train_bic <- c(mpr_bic_val,
                           spr_bic_val,
                           ladlassoic_bic_val,
                           rlm_bic_val,
                           bamlss_bic_val)
          },
          "normal" = {
            # alassoic -- estimate alpha_0 and kappa
            alassoic_sgnd_est_scale <- estimate_scale_kappa_func(theta_spr_scale = get_thetas_vec(fit = "alassoic",
                                                                                                  df = thetas_dataf_scale,
                                                                                                  model_type = "no_nu"),
                                                                 nu_0 = kappa_to_nu(2, kappa_omega), # initial value for opt
                                                                 x_scale = x_scale_train,
                                                                 y = y_raw_train,
                                                                 tau = tau_T_mpr,
                                                                 list_general = list_family$general,
                                                                 method_c_tilde = method_c_tilde,
                                                                 kappa_omega = kappa_omega,
                                                                 iterlim_nlm = iterlim_nlm,
                                                                 stepmax = NA)
            alassoic_bic_val <- info_criterion_sgnd(theta = alassoic_sgnd_est_scale / x_sd_train_spr, # unscale
                                                    x1 = x_raw_train,
                                                    x2 = x3_train,
                                                    x3 = x3_train,
                                                    y = y_raw_train,
                                                    tau = tau_T_mpr,
                                                    list_general = list_family$general,
                                                    method_c_tilde = method_c_tilde,
                                                    kappa_omega = kappa_omega,
                                                    df = (fitres$fits$fit_alassoic %>% slice_min(ic) %>% pull(df)),
                                                    lambda = "log(n)")
            # lm -- estimate alpha_0 and kappa
            lm_sgnd_est_scale <- estimate_scale_kappa_func(theta_spr_scale = get_thetas_vec(fit = "lm",
                                                                                            df = thetas_dataf_scale,
                                                                                            model_type = "no_nu"),
                                                           nu_0 = kappa_to_nu(2, kappa_omega), # initial value for opt
                                                           x_scale = x_scale_train,
                                                           y = y_raw_train,
                                                           tau = tau_T_mpr,
                                                           list_general = list_family$general,
                                                           method_c_tilde = method_c_tilde,
                                                           kappa_omega = kappa_omega,
                                                           iterlim_nlm = iterlim_nlm,
                                                           stepmax = NA)

            lm_bic_val <- info_criterion_sgnd(theta = lm_sgnd_est_scale / x_sd_train_spr, # unscale
                                              x1 = x_raw_train,
                                              x2 = x3_train,
                                              x3 = x3_train,
                                              y = y_raw_train,
                                              tau = tau_T_mpr,
                                              list_general = list_family$general,
                                              method_c_tilde = method_c_tilde,
                                              kappa_omega = kappa_omega,
                                              df = p, # no penalization
                                              lambda = "log(n)")

            train_bic <- c(mpr_bic_val,
                           spr_bic_val,
                           alassoic_bic_val,
                           lm_bic_val,
                           bamlss_bic_val)
          }
  )

  # include train bic in table
  prediction_metrics_dataf <- prediction_metrics_dataf %>%
    add_column(bic = train_bic,
               .before = "mse")

  # s width categories
  s_vec_mpr <- s_list[[which(names_fits == "mpr")]] # based on test set

  index_cat <- case_when(s_vec_mpr <= lower_cutoff ~ 1,
                         (s_vec_mpr > lower_cutoff) & (s_vec_mpr <= upper_cutoff) ~ 2,
                         s_vec_mpr > upper_cutoff ~ 3
  )

  prediction_intervals_list_s_cat <-  prediction_intervals_list %>%
    map(~ {
      .x %>%
        select(s,
               P) %>%
        add_column(index_cat = index_cat) # include indices in prediction intervals dataf
    })

  # get mean for each category
  mtrue_p_s_cat_dataf <- prediction_intervals_list_s_cat %>%
    map(~ {
      .x %>%
        group_by(index_cat) %>%
        summarise(mean_P = mean(P), .groups = "keep")
    }) %>%
    map2(.x = ., .y = names_fits, ~ {
      .x %>%
        add_column(
          fit_type = .y,
          .before = 1
        )
    }) %>%
    rbindlist()

  # return list --
  list(
    "thetas_dataf" = thetas_dataf,
    "prediction_metrics_dataf" = prediction_metrics_dataf,
    "s_width_dataf" = s_width_dataf,
    "s_width_categories_dataf" = mtrue_p_s_cat_dataf
  )
}

out_of_sample_metrics_parallel <- function(index_oos_vec, # 1:10
                                           type,
                                           x_raw_incl_int,
                                           y_raw,
                                           indices_oos,
                                           less_than_value,
                                           less_than_value_bamlss,
                                           optimizer,
                                           iterlim_nlm,
                                           stepmax_nlm,
                                           tau_1_mpr,
                                           tau_T_mpr,
                                           tau_1_spr,
                                           tau_T_spr,
                                           epsilon_1,
                                           epsilon_T,
                                           steps_T,
                                           list_family,
                                           method_c_tilde,
                                           method_c_tilde_deriv,
                                           algorithm,
                                           h_kappa,
                                           kappa_omega,
                                           lambda_beta,
                                           lambda_alpha,
                                           lambda_nu,
                                           initial_step,
                                           max_step_it,
                                           tol,
                                           max_it,
                                           lower_cutoff,
                                           upper_cutoff,
                                           formula_bamlss) {
  res <- parLapply(
    cl,
    index_oos_vec,
    out_of_sample_metrics,
    type = type,
    x_raw_incl_int = x_raw_incl_int,
    y_raw = y_raw,
    indices_oos = indices_oos,
    less_than_value = less_than_value,
    less_than_value_bamlss = less_than_value_bamlss,
    optimizer = optimizer,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    tau_1_mpr = tau_1_mpr,
    tau_T_mpr = tau_T_mpr,
    tau_1_spr = tau_1_spr,
    tau_T_spr = tau_T_spr,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it = max_it,
    lower_cutoff = lower_cutoff,
    upper_cutoff = upper_cutoff,
    formula_bamlss = formula_bamlss
  )
  res
}

# Bootstrap SEE -----------------------------------------------------------
bootstrap_single_data <- function(indices,
                                  fit_type, # mpr or spr
                                  x_raw_int,
                                  y_raw,
                                  kappa_omega,
                                  method_initial_values,
                                  tau_1,
                                  tau_T,
                                  epsilon_1,
                                  epsilon_T,
                                  steps_T,
                                  list_family,
                                  method_c_tilde,
                                  iterlim_nlm,
                                  stepmax_nlm,
                                  lambda_beta,
                                  lambda_alpha,
                                  lambda_nu) {

  x_now <- x_raw_int[indices, ]
  y_now <- y_raw[indices]

  n_all <- nrow(x_raw_int)
  p <- ncol(x_raw_int) - 1

  ## Scale X ----
  names_x <- paste0("X_", 1:p)
  x_scale <- scale(x_now[, -1],
                   center = FALSE,
                   scale = apply(x_now[, -1], 2, sd)
  )
  x_scale <- cbind(rep(1, n_all), x_scale) # column of 1's for intercept
  x_sd <- apply(x_now[, -1], 2, sd) # save sd for later to transform back

  x1 <- x_scale

  switch (fit_type,
          "mpr" = {
            x2 <- x_scale
            x3 <- as.matrix(rep(1, n_all))
            x_sd_theta <- c(1, x_sd, 1, x_sd, 1)
          },
          "spr" = {
            x2 <- x3 <- as.matrix(rep(1, n_all))
            x_sd_theta <- c(1, x_sd, 1, 1) # spr version
          }
  )

  theta_init <- initial_values_func(
    x1 = x1, # data should be scaled already
    x2 = x2,
    x3 = x3,
    y = y_now,
    method_initial_values = method_initial_values,
    other_initial_values_vec = NA, # can manually input vector of initial values
    kappa_omega = kappa_omega,
    fix_nu_value_init = NA,
    tau_1 = tau_1,
    epsilon_1 = epsilon_1,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = NA,
    algorithm = NA,
    h_kappa = NA,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = NA,
    max_step_it = NA,
    tol = NA,
    max_it = NA
  )

  tele_nlm <- tryCatch(fitting_func_base(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y_now,
    optimizer = "nlm",
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    nu_profile_vec = NA,
    theta_init = theta_init, # initial values
    tau_1 = tau_1,
    tau_T = tau_T,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = NA,
    algorithm = NA,
    h_kappa = NA,
    kappa_omega = kappa_omega,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu,
    initial_step = NA,
    max_step_it = NA,
    tol = NA,
    max_it = NA
  ), error = function(err) NA)

  extract_nlm <- tryCatch(extract_theta_plike_val(
    tele_nlm,
    "nlm"
  ),
  error = function(err) NA
  )
  theta_nlm <- tryCatch(as.vector(extract_nlm$theta / x_sd_theta),
                        error = function(err) NA
  ) # unscale

  list("theta_nlm" = theta_nlm,
       "x_sd" = x_sd)
}

bootstrap_multiple_data <- function(indices_list,
                                    fit_type, # mpr or spr
                                    x_raw_int,
                                    y_raw,
                                    kappa_omega,
                                    method_initial_values,
                                    tau_1,
                                    tau_T,
                                    epsilon_1,
                                    epsilon_T,
                                    steps_T,
                                    list_family,
                                    method_c_tilde,
                                    iterlim_nlm,
                                    stepmax_nlm,
                                    lambda_beta,
                                    lambda_alpha,
                                    lambda_nu) {
  parLapply(
    cl,
    indices_list, # list of n_boots indices
    bootstrap_single_data,
    fit_type = fit_type,
    x_raw_int = x_raw_int,
    y_raw = y_raw,
    kappa_omega = kappa_omega,
    method_initial_values = method_initial_values,
    tau_1 = tau_1,
    tau_T = tau_T,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    iterlim_nlm = iterlim_nlm,
    stepmax_nlm = stepmax_nlm,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    lambda_nu = lambda_nu
  )
}


# bamlss ------------------------------------------------------------------
opt_lasso_MON <- function (x, y, start = NULL, adaptive = TRUE, lower = 0.001,
                           upper = 1000, nlambda = 100, lambda = NULL, multiple = FALSE,
                           verbose = TRUE, digits = 4, flush = TRUE, nu = NULL, stop.nu = NULL,
                           ridge = .Machine$double.eps^0.5, zeromodel = NULL, ...)
{
  method <- list(...)$method
  if (is.null(method))
    method <- 1
  if (is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, update = bfit_iwls, ...)
  start2 <- start
  if (lower < 0.00000000000000000001)
    lower <- 0.00000000000000000001
  lambdas <- if (is.null(lambda)) {
    exp(seq(log(upper), log(lower), length = nlambda))
  }
  else lambda
  lambdas <- rep(list(lambdas), length = length(x)-1) # -1
  lambdas <- append(lambdas, list(0)) # zero for nu
  names(lambdas) <- names(x)
  lambdas <- as.matrix(do.call(if (multiple)
    "expand.grid"
    else "cbind", lambdas))
  if (length(verbose) < 2)
    verbose <- c(verbose, FALSE)
  ia <- if (flush)
    interactive()
  else FALSE
  par <- list()
  ic <- NULL
  ptm <- proc.time()
  fuse <- NULL
  for (i in names(x)) {
    for (j in names(x[[i]]$smooth.construct)) {
      if (inherits(x[[i]]$smooth.construct[[j]], "lasso.smooth")) {
        x[[i]]$smooth.construct[[j]]$state$do.optim <- FALSE
        x[[i]]$smooth.construct[[j]]$fxsp <- TRUE
        fuse <- c(fuse, x[[i]]$smooth.construct[[j]]$fuse)
        if (adaptive) {
          tau2 <- get.par(x[[i]]$smooth.construct[[j]]$state$parameters,
                          "tau2")
          tau2 <- rep(1/ridge, length.out = length(tau2))
          x[[i]]$smooth.construct[[j]]$state$parameters <- set.par(x[[i]]$smooth.construct[[j]]$state$parameters,
                                                                   tau2, "tau2")
          x[[i]]$smooth.construct[[j]]$LAPEN <- x[[i]]$smooth.construct[[j]]$S
          x[[i]]$smooth.construct[[j]]$S <- list(diag(length(get.par(x[[i]]$smooth.construct[[j]]$state$parameters,
                                                                     "b"))))
        }
      }
    }
  }
  fuse <- if (is.null(fuse))
    FALSE
  else any(fuse)
  if (!is.null(nu))
    nu <- rep(nu, length.out = 2)
  if (!is.null(stop.nu))
    stop.nu <- rep(stop.nu, length.out = 2)
  if (adaptive & fuse) {
    if (verbose[1] & is.null(zeromodel))
      cat("Estimating adaptive weights\n---\n")
    if (is.null(zeromodel)) {
      if (method == 1) {
        zeromodel <- opt_bfit(x = x, y = y, start = start,
                              verbose = verbose[1], nu = nu[2], stop.nu = stop.nu[2],
                              ...)
      }
      else {
        zeromodel <- opt_optim(x = x, y = y, start = start,
                               verbose = verbose[1], ...)
      }
    }
    x <- lasso_transform(x, zeromodel, nobs = nrow(y))
  }
  else {
    if (!is.null(zeromodel)) {
      x <- lasso_transform(x, zeromodel, nobs = nrow(y))
    }
  }
  for (l in 1:nrow(lambdas)) {
    if (l > 1)
      start <- unlist(par[[l - 1]])
    tau2 <- NULL
    for (i in names(x)) {
      for (j in names(x[[i]]$smooth.construct)) {
        if (inherits(x[[i]]$smooth.construct[[j]], "lasso.smooth")) {
          tau2 <- get.par(x[[i]]$smooth.construct[[j]]$state$parameters,
                          "tau2")
          nt <- names(tau2)
          tau2 <- rep(1/lambdas[l, i], length.out = length(tau2))
          names(tau2) <- paste(i, "s", x[[i]]$smooth.construct[[j]]$label,
                               nt, sep = ".")
          if (!is.null(start) & (l > 1)) {
            if (all(names(tau2) %in% names(start))) {
              start[names(tau2)] <- tau2
            }
            else {
              start <- c(start, tau2)
            }
          }
          else {
            start <- c(start, tau2)
          }
        }
      }
    }
    if ((l < 2) & !is.null(start2)) {
      start <- c(start, start2)
      start <- start[!duplicated(names(start))]
    }
    if (method == 1) {
      b <- opt_bfit(x = x, y = y, start = start, verbose = verbose[2],
                    nu = nu[2], stop.nu = stop.nu[2], ...)
    }
    else {
      b <- opt_optim(x = x, y = y, start = start, verbose = verbose[2],
                     ...)
    }
    nic <- grep("ic", names(b), value = TRUE, ignore.case = TRUE)
    if (!length(nic)) {
      b$edf <- sum(abs(unlist(b$parameters)) > .Machine$double.eps^0.25)
      b$BIC <- -2 * b$logLik + b$edf * log(nrow(y))
    }
    nic <- grep("ic", names(b), value = TRUE, ignore.case = TRUE)
    par[[l]] <- unlist(b$parameters)
    mstats <- c(b$logLik, b$logPost, b[[nic]], b[["edf"]])
    names(mstats) <- c("logLik", "logPost", nic, "edf")
    ic <- rbind(ic, mstats)
    if (!is.null(list(...)$track)) {
      plot(ic[, nic] ~ c(1:l), type = "l", xlab = "Iteration",
           ylab = nic)
    }
    if (!is.null(stop.nu)) {
      if (l > stop.nu)
        nu <- NULL
    }
    if (verbose[1]) {
      cat(if (ia)
        "\r"
        else if (l > 1)
          "\n"
        else NULL)
      vtxt <- paste(nic, " ", bamlss:::fmt(b[[nic]], width = 8,
                                           digits = digits), " edf ", bamlss:::fmt(mstats["edf"],
                                                                                   width = 6, digits = digits), " lambda ", paste(bamlss:::fmt(if (!multiple)
                                                                                     lambdas[l, 1]
                                                                                     else lambdas[l, ], width = 6, digits = digits),
                                                                                     collapse = ","), " iteration ", formatC(l, width = nchar(nlambda)),
                    sep = "")
      cat(vtxt)
      if (.Platform$OS.type != "unix" & ia)
        flush.console()
    }
  }
  elapsed <- c(proc.time() - ptm)[3]
  if (verbose[1]) {
    et <- if (elapsed > 60) {
      paste(formatC(format(round(elapsed/60, 2), nsmall = 2),
                    width = 5), "min", sep = "")
    }
    else paste(formatC(format(round(elapsed, 2), nsmall = 2),
                       width = 5), "sec", sep = "")
    cat("\nelapsed time: ", et, "\n", sep = "")
  }
  colnames(lambdas) <- paste("lambda", names(x), sep = ".")
  ic <- cbind(ic, lambda = lambdas)
  rownames(ic) <- NULL
  attr(ic, "multiple") <- multiple
  class(ic) <- c("lasso.stats", "matrix")
  list(parameters = do.call("rbind", par), lasso.stats = ic,
       nobs = nrow(y))
}


extract_coef_bamlss <- function(fit_bamlss,
                                names_coef_labels) {

  levels_coef <- c(paste0("beta", "_", c("0", names_coef_labels)),
                   paste0("alpha", "_", c("0", names_coef_labels)),
                   "nu_0")

  coef_mat <- coef(fit_bamlss, mstop = lasso_stop(fit_bamlss))
  coef_mat %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_tibble() %>%
    rename(type = rowname) %>%
    filter(!(str_detect(type, c("tau")) | str_detect(type, c("alpha")))) %>%
    select(type, Mean) %>%
    mutate(coef_type = sub("\\..*", "", type),
           coef_name = sub('.*\\.', '', type),
           .keep = "unused",
           .before = 1) %>%
    mutate(coef_type = case_when(
      coef_type == "mu" ~ "beta",
      coef_type == "sigma" ~ "alpha",
      TRUE ~ coef_type
    )) %>%
    mutate(coef = paste0(coef_type, "_", coef_name),
           coef = case_when(
             str_detect(coef, "(Intercept)") ~ str_replace(coef, "\\(Intercept\\)", "0"),
             TRUE ~ coef
           ),
           .before = "Mean") %>%
    mutate(coef = factor(coef, levels = levels_coef)) %>%
    arrange(coef) %>%
    select(-coef_type, -coef_name) %>%
    rename(value = Mean)

}

extract_credible_intervals_bamlss <- function(fit_bamlss,
                                              names_coef_labels) {
  levels_coef <- c(paste0("beta", "_", c("0", names_coef_labels)),
                   paste0("alpha", "_", c("0", names_coef_labels)),
                   "nu_0")

  coef_mat <- coef(fit_bamlss, mstop = lasso_stop(fit_bamlss))
  ci_mat <- coef_mat %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_tibble() %>%
    rename(type = rowname) %>%
    filter(!(str_detect(type, c("tau")) | str_detect(type, c("alpha")))) %>%
    select(-4) %>% # remove 50% column
    mutate(
      coef_type = sub("\\..*", "", type),
      coef_name = sub(".*\\.", "", type),
      .keep = "unused",
      .before = 1
    ) %>%
    mutate(coef_type = case_when(
      coef_type == "mu" ~ "beta",
      coef_type == "sigma" ~ "alpha",
      TRUE ~ coef_type
    )) %>%
    mutate(
      coef = paste0(coef_type, "_", coef_name),
      coef = case_when(
        str_detect(coef, "(Intercept)") ~ str_replace(coef, "\\(Intercept\\)", "0"),
        TRUE ~ coef
      ),
      .before = "Mean"
    ) %>%
    mutate(coef = factor(coef, levels = levels_coef)) %>%
    arrange(coef) %>%
    select(-coef_type, -coef_name) %>%
    rename_all(~ c("coef", "estimate", "lower", "upper"))


  # significant?
  ci_sig_mat_scale <- ci_mat %>%
    mutate(signif = case_when(
      sign(lower) == sign(upper) ~ 1, # doesn't cross zero -> significant
      sign(lower) != sign(upper) ~ 0 # if different signs, then crosses zero
    )) # scaled estimates
  ci_sig_mat_scale
}

# Estimate kappa and scale for spr fits like lasso ----
estimate_scale_kappa_func <- function(theta_spr_scale,
                                      nu_0, # initial value for nu
                                      x_scale,
                                      y,
                                      tau,
                                      list_general,
                                      method_c_tilde,
                                      kappa_omega,
                                      iterlim_nlm,
                                      stepmax) {
  p <- ncol(x_scale) - 1
  x3 <- as.matrix(rep(1, nrow(x_scale)))
  res <- nlm_fixed_unpenalized(p = c(theta_spr_scale, nu_0),
                               fix_indices = 1:(p + 1), # estimate alpha_0 and nu_0
                               x1 = x_scale,
                               x2 = x3,
                               x3 = x3,
                               y = y,
                               tau = tau,
                               list_general = list_general,
                               method_c_tilde = method_c_tilde,
                               kappa_omega = kappa_omega,
                               iterlim_nlm = iterlim_nlm,
                               stepmax = stepmax)
  res$estimate

}
