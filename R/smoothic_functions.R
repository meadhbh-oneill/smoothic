#' @title Variable Selection Using a Smooth Information Criterion (SIC)
#'
#' @description Implements the SIC \eqn{\epsilon}-telescope method, either using
#' single or multi-parameter regression. Returns estimated coefficients, estimated
#' standard errors (SEE) and the value of the penalized likelihood function.
#' Note that the function will scale the predictors to have unit variance, however,
#' the final estimates are converted back to their original scale.
#'
#' @param formula An object of class \code{"\link{formula}"}: a two-sided object
#' with response on the left hand side and the model variables on the right hand side.
#'
#' @param data A data frame containing the variables in the model; the data frame
#' should be unstandardized.
#'
#' @param model The type of regression to be implemented, either \code{model = "mpr"}
#' for multi-parameter regression, or \code{model = "spr"} for single parameter
#' regression (i.e., classical normal linear regression). Defaults to \code{model="mpr"}.
#' @param lambda Value of penalty tuning parameter. Suggested values are
#' \code{"log(n)"} and \code{"2"} for the BIC and AIC respectively. Defaults to
#' \code{lambda ="log(n)"} for the BIC case.
#' @param epsilon_1 Starting value for \eqn{\epsilon}-telescope. Defaults to 10.
#' @param epsilon_T Final value for \eqn{\epsilon}-telescope. Defaults to
#' \code{1e-05}.
#' @param steps_T Number of steps in \eqn{\epsilon}-telescope. Defaults to 100.
#' @param zero_tol Coefficients below this value are treated as being zero.
#' Defaults to \code{1e-08}.
#' @param tol Convergence tolerance for the optimization. Defaults to
#' \code{1e-08}.
#' @param max_it Maximum number of iterations to performed before the
#' optimization is terminated. Defaults to \code{1e+04}.
#' @param initial_step Initial step length for step halving in Newton-Raphson
#' algorithm. Defaults to 10.
#' @param max_step_it Maximum allowable number of steps to take for step
#' halving in Newton-Raphson algorithm. Defaults to \code{1e+03}.
#'
#' @return A list with estimates and estimated standard errors.
#' \itemize{
#'   \item \code{coefficients} - vector of coefficients.
#'   \item \code{see} - vector of estimated standard errors.
#'   \item \code{model} - the matched type of model which is called.
#'   \item \code{plike} - value of the penalized likelihood function.
#'   }
#'
#' @author Meadhbh O'Neill
#'
#' @references O'Neill, M. and Burke, K. (2021) Variable Selection Using a Smooth
#' Information Criterion for Multi-Parameter Regression Models. arXiv:2110.02643 [stat.ME]
#'
#' @examples
#' # Sniffer Data --------------------
#' # MPR Model ----
#' results <- smoothic(
#'   formula = y ~ .,
#'   data = sniffer
#'   model = "mpr"
#' )
#' summary(results)
#'
#' @importFrom stats sd lm
#'
#' @export

smoothic <- function(formula,
                     data,
                     model = "mpr", # either "mpr" or "spr"
                     lambda = "log(n)", # lambda_beta = lambda_alpha
                     epsilon_1 = 10,
                     epsilon_T = 1e-05,
                     steps_T = 100,
                     zero_tol = 1e-08, # values less than this treated as zero
                     tol = 1e-08,
                     max_it = 1e+04,
                     initial_step = 10,
                     max_step_it = 1e+03) {
  # Formula object ----
  # x
  x <- model.matrix(formula, data = data)[, -1] # remove column of 1s

  # y
  y <- model.extract(model.frame(formula, data = data), "response")


  if (all(x[, 1] == 1)) {
    stop("Error: Make sure column of 1's for intercept is not included in x")
  }

  if (any(is.na(x) == TRUE) | any(is.na(y) == TRUE)) {
    stop("Error: Make sure data does not contain any NAs.")
  }

  n <- length(y)
  p <- ncol(x) # not including intercept

  if (is.null(colnames(x))) {
    colnames_x <- paste0("X_", 0:p)
  } else {
    colnames_x <- paste0(c("intercept", colnames(x)), "_", 0:p)
  }

  # Scale x ----
  x_scale <- scale(x,
    center = FALSE,
    scale = apply(x, 2, sd)
  )
  x_scale <- cbind(rep(1, n), x_scale) # column of 1's for intercept
  x_sd <- apply(x, 2, sd) # save sd for later to transform back

  # Initial values ----
  lm_fit <- lm(y ~ x_scale[, -1]) # remove intercept column
  lm_coef_sig <- c(
    unname(lm_fit$coefficients),
    log((summary(lm_fit)$sigma)^2)
  )
  theta_init <- c(lm_coef_sig, rep(0, p)) # vector of zeros as initial values for alpha parameters

  # Epsilon telescope vector ----
  eps_tele_vec <- rev(exp(seq(log(epsilon_T), log(epsilon_1), length = steps_T)))

  # MPR ---------------------------------
  # model: "mpr" ----
  if (model == "mpr") {
    # Telescope ----
    tele_mat_scale <- telescope_mpr_base(
      x = x_scale,
      y = y,
      theta_init = theta_init,
      eps_tele_vec = eps_tele_vec,
      lambda = lambda,
      initial_step = initial_step,
      max_step_it = max_step_it,
      tol = tol,
      max_it = max_it
    )
    theta_scale <- tele_mat_scale[steps_T, (2:((2 * p) + 3))] # extract estimates
    theta <- as.vector(theta_scale / c(1, x_sd, 1, x_sd)) # unscale (convert back)
    theta_final <- theta # treat some as zero

    zero_pos <- which(abs(theta_final) < zero_tol)
    theta_final[zero_pos] <- 0 # treat values less than zero_tol to zero

    names_coef <- c(
      paste0(colnames_x, "_beta"),
      paste0(colnames_x, "_alpha")
    )
    names(theta_final) <- names_coef
    names(theta) <- names_coef

    # Get standard errors ----
    info_mat_list <- information_matrices_mpr(
      theta = theta,
      x = as.matrix(cbind(rep(1, n), x)), # include col of 1's in raw data
      y = y,
      lambda = lambda,
      epsilon = epsilon_T # final epsilon
    )

    see_vec <- get_see_func(info_mat_list) # calculate standard errors
    see_vec[zero_pos] <- 0 # if coef treated as zero then change SEE to zero
    names(see_vec) <- names_coef

    # Penalized likelihood value ----
    plike <- tele_mat_scale[steps_T, "criterion_val"]
  } else if (model == "spr") {
    # MPR ---------------------------------
    # model: "spr" ----
    # Telescope ----
    tele_mat_scale <- telescope_spr_base(
      x = x_scale,
      y = y,
      theta_init = lm_coef_sig,
      eps_tele_vec = eps_tele_vec,
      lambda = lambda,
      initial_step = initial_step,
      max_step_it = max_step_it,
      tol = tol,
      max_it = max_it
    )
    theta_scale <- tele_mat_scale[steps_T, (2:(p + 3))] # extract estimates
    theta <- as.vector(theta_scale / c(1, x_sd, 1)) # unscale (convert back)
    theta_final <- theta

    zero_pos <- which(abs(theta_final) < zero_tol)
    theta_final[zero_pos] <- 0 # set values less than zero_tol to zero

    names_coef <- c(paste0(colnames_x, "_beta"), "0_alpha")

    names(theta_final) <- names_coef
    names(theta) <- names_coef

    # Get standard errors ----
    info_mat_list <- information_matrices_spr(
      theta_incl0 = c(theta, rep(0, p)), # include 0's for alpha vector for this func
      x = as.matrix(cbind(rep(1, n), x)), # include col of 1's in raw data
      y = y,
      lambda = lambda,
      epsilon = epsilon_T # final epsilon
    )

    see_vec <- get_see_func(info_mat_list) # calculate standard errors
    see_vec[zero_pos] <- 0 # if coef treated as zero then change SEE to zero
    names(see_vec) <- names_coef

    # Penalized likelihood value ----
    plike <- tele_mat_scale[steps_T, "criterion_val"]
  }

  # Output ----
  out <- list(
    "coefficients" = theta_final,
    "see" = see_vec,
    "model" = model,
    "plike" = plike
  )
  class(out) <- "smoothic"
  out
}

#' @title Summarising Smooth Information Criterion (SIC) Fits
#'
#' @aliases summary.smoothic print.summary.smoothic
#'
#' @description \code{summary} method class \dQuote{\code{smoothic}}
#'
#' @param object an object of class \dQuote{\code{smoothic}} which is the result
#' of a call to \code{\link{smoothic}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{model} - the matched model from the \code{smoothic} object.
#'   \item \code{coefmat} - a typical coefficient matrix whose columns are the
#'   estimated regression coefficients, estimated standard errors (SEE) and p-values.
#'   \item \code{plike} - value of the penalized likelihood function.
#'   }
#'
#' @author Meadhbh O'Neill
#'
#' @examples
#' # Sniffer Data --------------------
#' # MPR Model ----
#' results <- smoothic(
#'   formula = y ~ .,
#'   data = sniffer
#'   model = "mpr"
#' )
#' summary(results)
#'
#' @importFrom stats pnorm
#'
#' @export

summary.smoothic <- function(object, ...) {
  coefficients <- object$coefficients
  see <- object$see

  zeropos <- which(coefficients == 0)

  coefficients[zeropos] <- NA
  see[zeropos] <- NA

  zval <- coefficients / see
  pval <- 1 * pnorm(abs(zval), lower.tail = FALSE)

  coefmat <- cbind(
    Estimate = coefficients,
    SEE = see,
    Z = zval,
    Pvalue = pval
  )
  out <- list(
    model = object$model,
    coefmat = coefmat,
    plike = unname(object$plike)
  )
  class(out) <- "summary.smoothic"
  out
}



# print.smoothic ----------------------------------------------------------
#' @export
print.smoothic <- function(object, ...) {
  cat("Model:\n")
  print(object$model)
  cat("\nCoefficients:\n")
  print(object$coefficients)
}

# print.summary.smoothic --------------------------------------------------
#' @importFrom stats printCoefmat
#' @export
print.summary.smoothic <- function(object, ...) {
  cat("Model:\n")
  print(object$model)
  cat("\nCoefficients:\n")
  printCoefmat(object$coefmat,
    cs.ind = 1:2,
    tst.ind = 3,
    P.values = TRUE,
    has.Pvalue = TRUE,
    signif.legend = TRUE,
    na.print = "0" # change NA to 0 for printing
  )
  cat("Penalized Likelihood:\n")
  print(object$plike) # BIC or AIC = -2*plike
}

# Iterative Loop ----------------------------------------------------------
# * For both MPR & SPR ----------------------------------------------------
it_loop <- function(theta,
                    x,
                    y,
                    tol,
                    max_it,
                    FUN,
                    ...) {
  theta1 <- theta + 1
  it <- 0
  max_step_reached <- NULL
  steps <- NULL
  while (max(abs(theta1 - theta)) > tol & it < max_it) {
    theta1 <- theta
    fit <- FUN(theta1, x, y, ...)
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

# * Standard Errors -------------------------------------------------------
get_see_func <- function(list_of_2_info_matrices) {
  observed_information_penalized <- list_of_2_info_matrices$observed_information_penalized
  observed_information_unpenalized <- list_of_2_info_matrices$observed_information_unpenalized

  inverse_observed_information_penalized <- solve(observed_information_penalized)

  # Square root of the diagonal of the variance-covariance matrix
  see <- sqrt(diag(inverse_observed_information_penalized %*% observed_information_unpenalized %*% inverse_observed_information_penalized)) # sandwich formula

  see
}


# MPR Functions -----------------------------------------------------------
# * MPR Normal Likelihood -------------------------------------------------
normallike <- function(theta,
                       x,
                       y) {
  n <- length(y)

  beta <- theta[1:(0.5 * length(theta))] # first half of vector = beta
  alpha <- theta[-(1:(0.5 * length(theta)))] # second half of vector = alpha

  p <- length(beta) - 1 # p same for alpha and beta

  mu <- x %*% beta
  phi <- exp(x %*% alpha)

  like <- -(n / 2) * log(2 * pi) - (1 / 2) * (sum(x %*% alpha)) - (1 / 2) * (sum((1 / phi) * ((y - mu)^2)))
  like # maximize likelihood
}


# * MPR Negative Normal Likelihood ----------------------------------------
neg_normallike <- function(...) {
  -normallike(...)
}

# * MPR Penalized Likelihood ----------------------------------------------
penalizedlike <- function(theta,
                          x,
                          y,
                          lambda,
                          epsilon) {
  # maximise this
  n <- length(y)

  beta <- theta[1:(0.5 * length(theta))] # first half of vector = beta
  alpha <- theta[-(1:(0.5 * length(theta)))] # second half of vector = alpha

  p <- length(beta) - 1 # p same for alpha and beta

  lambda <- (eval(parse(text = lambda)))

  # remove intercepts as won't be penalized
  beta1_p <- beta[-1]
  alpha1_p <- alpha[-1]

  # both sets of coefficients penalized
  p_like <- normallike(theta, x, y) - ((lambda / 2) * (1 + (sum((beta1_p^2) / (beta1_p^2 + epsilon^2))))) - ((lambda / 2) * (1 + (sum((alpha1_p^2) / (alpha1_p^2 + epsilon^2)))))
  p_like # divided by two for penalized setup, +1 for estimation of intercepts
}


# * MPR Negative Penalized Likelihood -------------------------------------
neg_penalizedlike <- function(...) {
  -penalizedlike(...)
}

# * MPR Newton Raphson No Penalty -----------------------------------------
mpr_normal <- function(theta,
                       x,
                       y,
                       initial_step,
                       max_step_it) {
  n <- length(y)

  beta <- theta[1:(0.5 * length(theta))] # first half of vector = beta
  alpha <- theta[-(1:(0.5 * length(theta)))] # second half of vector = alpha

  p <- length(beta) - 1 # p same for alpha and beta

  mu <- x %*% beta
  phi <- exp(x %*% alpha)

  tx <- t(x) # save transpose of x as an object

  ## Score
  ## Beta
  z_beta <- ((1 / phi) * (y - mu))
  score_beta <- tx %*% z_beta

  ## Alpha
  z_alpha <- (-1 / 2) + (1 / (2 * phi)) * (((y - mu)^2))
  score_alpha <- tx %*% z_alpha

  score <- c(score_beta, score_alpha)

  ## Information Matrix
  W_beta <- (1 / phi)
  I_beta <- (t(x * c(W_beta))) %*% x

  W_alpha <- (1 / (2 * phi)) * ((y - mu)^2)
  I_alpha <- (t(x * c(W_alpha))) %*% x

  W_beta_alpha <- rep(0, n) # cross-derivatives = 0 for optimization
  # W_beta_alpha <- (1 / phi) * (y - mu) # true derivative
  I_beta_alpha <- t(x * c(W_beta_alpha)) %*% x
  I_alpha_beta <- t(I_beta_alpha)

  information <- rbind(
    cbind(I_beta, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha)
  )

  delta <- solve(information, score)

  ## Step halving
  like_current <- normallike(theta, x, y)
  like_new <- like_current - 1
  j <- 0
  while (like_new < like_current & j < max_step_it) {
    theta_new <- theta + ((initial_step * delta) / (2^j))
    like_new <- normallike(theta_new, x, y)
    j <- j + 1
  }
  max_step_reached <- ifelse(j == max_step_it, TRUE, FALSE)
  return(list(
    "estimate" = theta_new, "maximum" = like_new,
    "steps" = c(max_step_reached, j)
  ))
}

# * MPR Newton Raphson Penalty --------------------------------------------
mpr_penalty <- function(theta,
                        x,
                        y,
                        lambda,
                        epsilon,
                        initial_step,
                        max_step_it) {
  n <- length(y)

  beta <- theta[1:(0.5 * length(theta))] # first half of vector = beta
  alpha <- theta[-(1:(0.5 * length(theta)))] # second half of vector = alpha

  p <- length(beta) - 1 # p same for alpha and beta

  beta1_p <- beta[-1] # remove intercepts as won't be penalized
  alpha1_p <- alpha[-1]

  mu <- x %*% beta
  phi <- exp(x %*% alpha)

  tx <- t(x) # save transpose of x as an object

  lambda <- (eval(parse(text = lambda)))

  ## Score
  ## Beta
  z_beta <- ((1 / phi) * (y - mu))
  v_beta <- c(0, (lambda / 2) * ((2 * beta1_p * (epsilon^2)) / ((beta1_p^2 + epsilon^2)^2)))
  score_beta <- (tx %*% z_beta) - v_beta

  ## Alpha
  z_alpha <- (-1 / 2) + (1 / (2 * phi)) * (((y - mu)^2))
  v_alpha <- c(0, (lambda / 2) * ((2 * alpha1_p * (epsilon^2)) / ((alpha1_p^2 + epsilon^2)^2)))
  score_alpha <- (tx %*% z_alpha) - v_alpha

  score <- c(score_beta, score_alpha)

  ## Information Matrix
  # Beta
  W_beta <- (1 / phi)
  Sigma_beta <- diag(c(
    0,
    (lambda / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
  ),
  nrow = p + 1,
  ncol = p + 1
  )
  I_beta <- ((t(x * c(W_beta))) %*% x) + Sigma_beta

  # Alpha
  W_alpha <- (1 / (2 * phi)) * ((y - mu)^2)
  Sigma_alpha <- diag(c(
    0,
    (lambda / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (alpha1_p^2))) / ((alpha1_p^2 + epsilon^2)^3))
  ),
  nrow = p + 1,
  ncol = p + 1
  )
  I_alpha <- ((t(x * c(W_alpha))) %*% x) + Sigma_alpha

  # Cross deriv
  W_beta_alpha <- rep(0, n) # Can change depending on orthogonality of parameters
  # W_beta_alpha <- (1 / phi) * (y - mu) # true derivative
  I_beta_alpha <- t(x * c(W_beta_alpha)) %*% x
  I_alpha_beta <- t(I_beta_alpha)

  information <- rbind(
    cbind(I_beta, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha)
  )

  delta <- solve(information, score)

  ## Step halving
  like_current <- penalizedlike(theta, x, y, lambda, epsilon)
  like_new <- like_current - 1
  j <- 0
  while (like_new < like_current & j < max_step_it) {
    theta_new <- theta + ((initial_step * delta) / (2^j))
    like_new <- penalizedlike(theta_new, x, y, lambda, epsilon)
    j <- j + 1
  }
  max_step_reached <- ifelse(j == max_step_it, TRUE, FALSE)
  return(list(
    "estimate" = theta_new, "maximum" = like_new,
    "steps" = c(max_step_reached, j)
  ))
}


# * MPR Information Matrices ----------------------------------------------
information_matrices_mpr <- function(theta,
                                     x,
                                     y,
                                     lambda,
                                     epsilon) {
  n <- length(y)

  beta <- theta[1:(0.5 * length(theta))] # first half of vector = beta
  alpha <- theta[-(1:(0.5 * length(theta)))] # second half of vector = alpha

  p <- length(beta) - 1 # p same for alpha and beta

  beta1_p <- beta[-1] # remove intercepts as won't be penalized
  alpha1_p <- alpha[-1]

  mu <- x %*% beta
  phi <- exp(x %*% alpha)

  tx <- t(x) # save transpose of x as an object

  lambda <- (eval(parse(text = lambda)))

  ## Penalized Observed Information Matrix
  # Beta
  W_beta <- (1 / phi)
  Sigma_beta <- diag(c(
    0,
    (lambda / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
  ),
  nrow = p + 1,
  ncol = p + 1
  )
  I_beta <- ((t(x * c(W_beta))) %*% x) + Sigma_beta

  # Alpha
  W_alpha <- (1 / (2 * phi)) * ((y - mu)^2)
  Sigma_alpha <- diag(c(
    0,
    (lambda / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (alpha1_p^2))) / ((alpha1_p^2 + epsilon^2)^3))
  ),
  nrow = p + 1,
  ncol = p + 1
  )
  I_alpha <- ((t(x * c(W_alpha))) %*% x) + Sigma_alpha

  # Cross deriv
  # W_beta_alpha <- rep(0, n) # Can change depending on orthogonality of parameters
  W_beta_alpha <- (1 / phi) * (y - mu) #  negative derivative wrt to beta and alpha
  I_beta_alpha <- t(x * c(W_beta_alpha)) %*% x
  I_alpha_beta <- t(I_beta_alpha)

  observed_information_penalized <- rbind(
    cbind(I_beta, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha)
  )


  ## Unpenalized Information evaluated at penalized estimates
  I_beta_unpen <- ((t(x * c(W_beta))) %*% x) # no sigma
  I_alpha_unpen <- ((t(x * c(W_alpha))) %*% x) # no sigma

  observed_information_unpenalized <- rbind(
    cbind(I_beta_unpen, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha_unpen)
  )


  list(
    "observed_information_penalized" = observed_information_penalized,
    "observed_information_unpenalized" = observed_information_unpenalized
  )
}

# * MPR Telescope ---------------------------------------------------------
telescope_mpr_base <- function(x,
                               y,
                               theta_init,
                               eps_tele_vec,
                               lambda,
                               initial_step,
                               max_step_it,
                               tol,
                               max_it) {
  p <- ncol(x) - 1
  t_initial_guess <- theta_init
  t_res_mat <- matrix(NA,
    nrow = length(eps_tele_vec),
    ncol = (length(t_initial_guess) + 5)
  )
  colnames(t_res_mat) <- c(
    "epsilon",
    paste0("beta_", 0:p),
    paste0("alpha_", 0:p),
    "criterion_val",
    "step_full",
    "it_full",
    "it"
  )
  ## Telescope
  for (i in eps_tele_vec) {
    pos <- which(eps_tele_vec == i)
    t_fit <- it_loop(
      theta = t_initial_guess,
      x = x,
      y = y,
      tol = tol,
      max_it = max_it,
      FUN = mpr_penalty,
      lambda = lambda,
      epsilon = i,
      initial_step = initial_step,
      max_step_it = max_step_it
    )
    t_initial_guess <- t_fit$estimate
    steps_full <- ifelse(any(t_fit$max_step_reached == 1), 1, 0) # 1 if true, 0 if false

    t_res_mat[pos, ] <- c(i, t_fit$estimate, t_fit$maximum, steps_full, t_fit$iterations)
  }
  t_res_mat
}



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# SPR Functions -----------------------------------------------------------
# * SPR Normal Likelihood -------------------------------------------------
normallike_spr <- function(theta,
                           x,
                           y) { # theta = c(beta_0, ..., beta_p, log(sigma_sq))
  n <- length(y)
  p <- ncol(x) - 1

  beta <- theta[1:(p + 1)]
  sigma_sq <- exp(theta[p + 2]) # last element of theta = alpha (phi = exp(alpha))

  mu <- x %*% beta
  like <- -(n / 2) * log(2 * pi) - ((n / 2) * log(sigma_sq)) - ((1 / (2 * (sigma_sq))) * sum((y - mu)^2))
  like # maximise likelihood
}

# * SPR Negative Normal Likelihood ----------------------------------------
neg_normallike_spr <- function(...) {
  -normallike_spr(...)
}


# * SPR Penalized Likelihood ----------------------------------------------
# No alphas so only betas penalized
penalizedlike_spr <- function(theta,
                              x,
                              y,
                              lambda,
                              epsilon) { # maximise this
  n <- length(y)
  p <- ncol(x) - 1

  lambda <- (eval(parse(text = lambda)))
  beta1_p <- theta[2:(p + 1)] # no intercept

  normallike_spr(theta, x, y) - ((lambda / 2) * (2 + (sum((beta1_p^2) / (beta1_p^2 + epsilon^2))))) # normallike - (lambda/2)*penalty, +2 for estimation of beta0 and sigma sq, # divide by 2 as penalized likelihood setup
}


# ** SPR Negative Penalized Likelihood ------------------------------------
neg_penalizedlike_spr <- function(...) {
  -penalizedlike_spr(...)
}


# ** SPR Newton Raphson No Penalty ----------------------------------------
spr_normal <- function(theta,
                       x,
                       y,
                       initial_step,
                       max_step_it) {
  n <- length(y)
  p <- ncol(x) - 1
  x2 <- rep(1, n) # column of 1s for alpha_0

  beta <- theta[1:(p + 1)]
  alpha_0 <- theta[p + 2]

  mu <- x %*% beta
  phi <- exp(alpha_0)

  tx <- t(x)
  tx2 <- t(x2)

  ## Score
  ## Beta
  z_beta <- ((1 / phi) * (y - mu))
  score_beta <- tx %*% z_beta

  ## Alpha_0
  z_alpha <- (-1 / 2) + (1 / (2 * phi)) * (((y - mu)^2))
  score_alpha <- tx2 %*% z_alpha

  score <- c(score_beta, score_alpha)

  ## Information Matrix (vector operations faster than matrix ops)
  W_beta <- rep((1 / phi), n)
  I_beta <- t(x * c(W_beta)) %*% x

  W_alpha <- (1 / (2 * phi)) * ((y - mu)^2)
  I_alpha <- t(x2 * c(W_alpha)) %*% x2

  W_beta_alpha <- rep(0, n) # cross derivative = 0 for optimization
  # W_beta_alpha <- (1 / phi) * (y - mu) # true derivative
  I_beta_alpha <- t(x * c(W_beta_alpha)) %*% x2
  I_alpha_beta <- t(I_beta_alpha)

  information <- rbind(
    cbind(I_beta, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha)
  )

  delta <- solve(information, score)

  ## Step halving
  like_current <- normallike_spr(theta, x, y)
  like_new <- like_current - 1 # initialize loop
  j <- 0
  while (like_new < like_current & j < max_step_it) {
    theta_new <- theta + ((initial_step * delta) / (2^j))
    like_new <- normallike_spr(theta_new, x, y)
    j <- j + 1
  }
  max_step_reached <- ifelse(j == max_step_it, TRUE, FALSE)
  return(list(
    "estimate" = theta_new, "maximum" = like_new,
    "steps" = c(max_step_reached, j)
  ))
}


# ** SPR Newton Raphson Penalty -------------------------------------------
spr_penalty <- function(theta,
                        x,
                        y,
                        lambda,
                        epsilon,
                        initial_step,
                        max_step_it) {
  n <- length(y)
  p <- ncol(x) - 1
  x2 <- rep(1, n) # column of 1s for alpha_0

  beta <- theta[1:(p + 1)]
  alpha_0 <- theta[p + 2]

  beta1_p <- beta[-1] # remove intercept as won't penalize

  mu <- x %*% beta
  phi <- exp(alpha_0)

  tx <- t(x)
  tx2 <- t(x2)

  lambda <- (eval(parse(text = lambda)))
  ## Score
  ## Beta
  z_beta <- ((1 / phi) * (y - mu))
  v_beta <- c(0, (lambda / 2) * ((2 * beta1_p * (epsilon^2)) / ((beta1_p^2 + epsilon^2)^2)))
  score_beta <- (tx %*% z_beta) - v_beta

  ## Alpha_0
  z_alpha <- (-1 / 2) + (1 / (2 * phi)) * (((y - mu)^2))
  score_alpha <- tx2 %*% z_alpha # no penalty on constant variance

  score <- c(score_beta, score_alpha)

  ## Information Matrix (vector operations faster than matrix ops)
  W_beta <- rep((1 / phi), n)
  Sigma_beta <- diag(c(
    0,
    (lambda / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
  ),
  nrow = p + 1,
  ncol = p + 1
  )
  I_beta <- ((t(x * c(W_beta))) %*% x) + Sigma_beta

  W_alpha <- (1 / (2 * phi)) * ((y - mu)^2)
  I_alpha <- t(x2 * c(W_alpha)) %*% x2 # no penalty on constant variance

  W_beta_alpha <- rep(0, n) # cross derivative = 0 for optimization
  # W_beta_alpha <- (1 / phi) * (y - mu) # true derivative
  I_beta_alpha <- t(x * c(W_beta_alpha)) %*% x2
  I_alpha_beta <- t(I_beta_alpha)

  information <- rbind(
    cbind(I_beta, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha)
  )

  delta <- solve(information, score)

  ## Step halving
  like_current <- penalizedlike_spr(theta, x, y, lambda, epsilon) # x = x1 here as alpha doesn't depend on covariates
  like_new <- like_current - 1 # initialize loop
  j <- 0
  while (like_new < like_current & j < max_step_it) {
    theta_new <- theta + ((initial_step * delta) / (2^j))
    like_new <- penalizedlike_spr(theta_new, x, y, lambda, epsilon)
    j <- j + 1
  }
  max_step_reached <- ifelse(j == max_step_it, TRUE, FALSE)
  return(list(
    "estimate" = theta_new, "maximum" = like_new,
    "steps" = c(max_step_reached, j)
  ))
}


# ** SPR Information Matrices ---------------------------------------------
information_matrices_spr <- function(theta_incl0,
                                     x,
                                     y,
                                     lambda,
                                     epsilon) {
  n <- length(y)
  p <- ncol(x) - 1
  x2 <- rep(1, n) # column of 1s for alpha_0

  ## Get beta and alpha_0 (remove remainder of zeros - just in this simulation case to make iteration easier)
  theta <- theta_incl0[1:(p + 2)] # beta and alpha 0

  beta <- theta[1:(p + 1)]
  alpha_0 <- theta[p + 2]

  beta1_p <- beta[-1] # remove intercept as won't penalize

  mu <- x %*% beta
  phi <- exp(alpha_0)

  tx <- t(x)
  tx2 <- t(x2)

  lambda <- (eval(parse(text = lambda)))

  ## Penalized Observed Information Matrix
  W_beta <- rep((1 / phi), n)
  Sigma_beta <- diag(c(
    0,
    (lambda / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
  ),
  nrow = p + 1,
  ncol = p + 1
  )
  I_beta <- ((t(x * c(W_beta))) %*% x) + Sigma_beta

  W_alpha <- (1 / (2 * phi)) * ((y - mu)^2)
  I_alpha <- t(x2 * c(W_alpha)) %*% x2 # no penalty on constant variance

  W_beta_alpha <- (1 / phi) * (y - mu) # true derivative
  I_beta_alpha <- t(x * c(W_beta_alpha)) %*% x2
  I_alpha_beta <- t(I_beta_alpha)

  observed_information_penalized <- rbind(
    cbind(I_beta, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha)
  )

  ## Unpenalized Information at Penalized Estimates
  I_beta_unpen <- ((t(x * c(W_beta))) %*% x) # no sigma
  # I_alpha already unpenalized for constant variance
  observed_information_unpenalized <- rbind(
    cbind(I_beta_unpen, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha)
  )

  list(
    "observed_information_penalized" = observed_information_penalized,
    "observed_information_unpenalized" = observed_information_unpenalized
  )
}


# ** SPR Telescope --------------------------------------------------------
telescope_spr_base <- function(x,
                               y,
                               theta_init,
                               eps_tele_vec,
                               lambda,
                               initial_step,
                               max_step_it,
                               tol,
                               max_it) {
  p <- ncol(x) - 1
  t_initial_guess <- theta_init
  t_res_mat <- matrix(NA,
    nrow = length(eps_tele_vec),
    ncol = (length(t_initial_guess) + 5)
  )
  colnames(t_res_mat) <- c(
    "epsilon",
    paste0("beta_", 0:p),
    "logsigmasq",
    "criterion_val",
    "step_full",
    "it_full",
    "it"
  )

  ## Telescope
  for (i in eps_tele_vec) {
    pos <- which(eps_tele_vec == i)
    t_fit <- it_loop(
      theta = t_initial_guess,
      x = x,
      y = y,
      tol = tol,
      max_it = max_it,
      FUN = spr_penalty,
      lambda = lambda,
      epsilon = i,
      initial_step = initial_step,
      max_step_it = max_step_it
    )
    t_initial_guess <- t_fit$estimate
    steps_full <- ifelse(any(t_fit$max_step_reached == 1), 1, 0)
    t_res_mat[pos, ] <- c(
      i,
      t_fit$estimate,
      t_fit$maximum,
      steps_full,
      t_fit$iterations
    )
  }
  t_res_mat
}
