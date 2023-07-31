#' @title Variable Selection Using a Smooth Information Criterion (SIC)
#'
#' @description Implements the SIC \eqn{\epsilon}-telescope method, either using
#' single or multiparameter regression. Returns estimated coefficients, estimated
#' standard errors (SEE) and the value of the penalized likelihood function.
#' Note that the function will scale the predictors to have unit variance, however,
#' the final estimates are converted back to their original scale.
#'
#' @param formula An object of class \code{"\link{formula}"}: a two-sided object
#' with response on the left hand side and the model variables on the right hand side.
#'
#' @param data A data frame containing the variables in the model; the data frame
#' should be unstandardized.
#' @param family The family of the model, default is \code{family = "sgnd"} for the
#' "Smooth Generalized Normal Distribution" where the shape parameter kappa is also
#' estimated. Classical regression with normally distributed errors is performed
#' when \code{family = "normal"}. If \code{family = "laplace"}, this corresponds to
#' a robust regression with errors from a Laplace-like distribution. If \code{family = "laplace"},
#' then the default value of \code{tau = 0.15}, which is used to approximate the absolute value
#' in the Laplace density function.
#' @param model The type of regression to be implemented, either \code{model = "mpr"}
#' for multiparameter regression (i.e., location and scale), or \code{model = "spr"} for single parameter
#' regression (i.e., location only). Defaults to \code{model="mpr"}.
#' @param lambda Value of penalty tuning parameter. Suggested values are
#' \code{"log(n)"} and \code{"2"} for the BIC and AIC respectively. Defaults to
#' \code{lambda ="log(n)"} for the BIC case.
#' @param epsilon_1 Starting value for \eqn{\epsilon}-telescope. Defaults to 10.
#' @param epsilon_T Final value for \eqn{\epsilon}-telescope. Defaults to
#' \code{1e-04}.
#' @param steps_T Number of steps in \eqn{\epsilon}-telescope. Defaults to 100, must be
#' greater than or equal to 10.
#' @param zero_tol Coefficients below this value are treated as being zero.
#' Defaults to \code{1e-05}.
#' @param max_it Maximum number of iterations to be performed before the
#' optimization is terminated. Defaults to \code{1e+04}.
#' @param optimizer The optimization procedure to be used. Defaults to
#' \code{optimizer = "nlm"}, where the \code{\link{nlm}} function from the
#' \bold{stats} package is used. This tends to be more stable than the manually
#'  coded Newton-Raphson procedure that is used when \code{optimizer = "manual"}.
#' @param kappa Optional user-supplied positive kappa value (> 0.2 to avoid
#' computational issues) if \code{family = "sgnd"}. If supplied, the shape parameter
#' kappa will be fixed to this value in the optimization. If not supplied, kappa is
#' estimated from the data.
#' @param tau Optional user-supplied positive smoothing parameter value in the
#' "Smooth Generalized Normal Distribution" if \code{family = "sgnd"} or
#' \code{family = "laplace"}. If not supplied then \code{tau = 0.15}. If \code{family = "normal"}
#' then \code{tau = 0} is used. Smaller values of \code{tau} bring the approximation closer to the
#' absolute value function, but this can cause the optimization to become unstable. Some issues with
#' standard error calculation with smaller values of \code{tau} when using the Laplace distribution in
#' the robust regression setting.
#' @param max_it_vec Optional vector of length \code{steps_T} that contains the maximum number of
#' iterations to be performed in each \eqn{\epsilon}-telescope step. If not supplied, \code{max_it}
#' is the maximum number of iterations performed for 10 steps and then the maximum number of iterations
#' to be performed reduces to 10 for the remainder of the telescope.
#' @param stepmax_nlm Optional maximum allowable scaled step length (positive scalar) to be passed to
#' \code{\link{nlm}} if \code{optimizer = "nlm"}. If not supplied, default values in
#' \code{\link{nlm}} are used.
#'
#'
#' @return A list with estimates and estimated standard errors.
#' \itemize{
#'   \item \code{coefficients} - vector of coefficients.
#'   \item \code{see} - vector of estimated standard errors.
#'   \item \code{model} - the matched type of model which is called.
#'   \item \code{plike} - value of the penalized likelihood function.
#'   \item \code{kappa} - value of the estimated/fixed shape parameter kappa if \code{family = "sgnd"}.
#'   }
#'
#' @author Meadhbh O'Neill
#'
#' @references O'Neill, M. and Burke, K. (2023) Variable selection using a smooth information
#' criterion for distributional regression models. <doi:10.1007/s11222-023-10204-8>
#'
#' O'Neill, M. and Burke, K. (2022) Robust Distributional Regression with
#' Automatic Variable Selection. <arXiv:2212.07317>
#'
#' @examples
#' # Sniffer Data --------------------
#' # MPR Model ----
#' results <- smoothic(
#'   formula = y ~ .,
#'   data = sniffer,
#'   family = "normal",
#'   model = "mpr"
#' )
#' summary(results)
#' @importFrom stats sd lm model.matrix model.frame model.extract nlm integrate
#' @importFrom MASS rlm
#' @importFrom numDeriv grad hessian
#'
#' @export

smoothic <- function(formula,
                     data,
                     family = "sgnd",
                     model = "mpr",
                     lambda = "log(n)",
                     epsilon_1 = 10,
                     epsilon_T = 1e-4,
                     steps_T = 100,
                     zero_tol = 1e-5,
                     max_it = 1e4,
                     optimizer = "nlm",
                     kappa, # if missing then it is estimated
                     tau, # if missing and sgnd then set to 0.15
                     max_it_vec, # if missing then uses max_it for each step in the telescope
                     stepmax_nlm # if missing then uses nlm defaults
) {
  cl <- match.call()

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

  # mpr or spr ----
  x1 <- x_scale
  x3 <- as.matrix(rep(1, n))

  switch(model,
    "mpr" = {
      x2 <- x_scale
      x_sd_theta <- c(1, x_sd, 1, x_sd, 1) # mpr version
      names_coef <- c(
        paste0(colnames_x, "_beta"),
        paste0(colnames_x, "_alpha"),
        "nu_0"
      )
    },
    "spr" = {
      x2 <- x3
      x_sd_theta <- c(1, x_sd, 1, 1) # spr version
      names_coef <- c(
        paste0(colnames_x, "_beta"),
        "alpha_0",
        "nu_0"
      )
    }
  )

  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1
  p3 <- ncol(x3) - 1


  # Inputs ----
  if (!missing(stepmax_nlm)) {
    if (!(stepmax_nlm > 0)) {
      stop("Error: stepmax_nlm must be a positive scalar")
    }
  }

  if (missing(stepmax_nlm)) { # if not supplied, then set to NA so default values are used in nlm function
    stepmax_nlm <- NA
  }

  if (steps_T < 10) {
    stop("Error: steps_T must be greater than 10")
  }

  if (!missing(max_it_vec)) {
    if (!(length(max_it_vec) == steps_T)) {
      stop("Error: max_it_vec vector must be same length as steps_T")
    }
  }

  if (missing(max_it_vec)) {
    max_it_vec <- c(rep(max_it, times = 10), rep(10, times = steps_T - 10)) # less steps as go through telescope
  }


  if (!(family %in% c("old_normal", "normal", "laplace", "sgnd"))) {
    stop("Error: please choose a family from 'normal', 'laplace' or 'sgnd'")
  }


  if (family == "old_normal" & !missing(kappa)) {
    stop("Error: for 'old_normal' family, kappa cannot be fixed, please choose 'sgnd' family")
  }

  if (family == "normal" & !missing(kappa)) {
    stop("Error: for 'normal' family, kappa cannot be fixed, please choose 'sgnd' family")
  }

  if (family == "laplace" & !missing(kappa)) {
    stop("Error: for 'laplace' family, kappa cannot be fixed: please choose 'sgnd' family")
  }

  # If family == "normal" then use "basic" original coding method ----
  if (family == "old_normal") {
    fit_mat <- fitting_func_normal(
      x1 = x1,
      x2 = x2,
      y = y,
      optimizer = optimizer,
      epsilon_1 = epsilon_1,
      epsilon_T = epsilon_T,
      steps_T = steps_T,
      lambda_beta = lambda,
      lambda_alpha = lambda,
      max_it_vec = max_it_vec,
      stepmax_nlm = stepmax_nlm
    ) # returns matrix
    fit_out <- extract_theta_plike_val(fit_res = fit_mat)

    fit_out$theta <- c(fit_out$theta, "nu_0" = NA) # include NA for nu_0

    kappa <- NA # needed for output of function
    tau <- NA
    kappa_omega <- NA
  } else if (family != "old_normal") {
    # Inputs ----
    fix_kappa_lgl <- FALSE # assume not supplied so will be estimated (not fixed)

    if (!missing(kappa)) { # if kappa value supplied in input, then will be fixed
      if (!(kappa >= 0.2)) {
        stop("Error: kappa must be greater than 0.2 (must be positive parameter, but values less than 0.2 tend to be unstable")
      }
      fix_kappa_lgl <- TRUE
    }

    if (!fix_kappa_lgl) {
      kappa <- NA # if not fixing, then set to NA to be estimated
    }

    if (!missing(tau)) {
      if (!(tau >= 0)) {
        stop("Error: tau must be positive")
      }
    }

    if (missing(tau)) { # if tau not supplied, then make value
      tau <- 0.15
    }

    if (family == "laplace") { # fix kappa to 1 if family == "laplace"
      kappa <- 1
      fix_kappa_lgl <- TRUE
      # tau = 0.15 unless otherwise specified
    }

    if (family == "normal") { # fix kappa to 2 if family == "normal" # new normal
      kappa <- 2
      fix_kappa_lgl <- TRUE
      tau <- 0 # tau = 0 for normal distribution, tau only required for absolute value when kappa != 2
    }

    kappa_omega <- 0.2

    # Fit model ----
    fit_mat <- fitting_func_pkg(
      x1 = x1,
      x2 = x2,
      x3 = x3,
      y = y,
      family = family,
      optimizer = optimizer,
      lambda = lambda,
      epsilon_1 = epsilon_1,
      epsilon_T = epsilon_T,
      steps_T = steps_T,
      max_it_vec = max_it_vec,
      kappa = kappa,
      tau = tau,
      fix_kappa_lgl = fix_kappa_lgl,
      kappa_omega = kappa_omega,
      stepmax_nlm = stepmax_nlm
    )
    fit_out <- extract_theta_plike_val(fit_res = fit_mat)
  }

  # Extract values ----
  # Estimates ----
  theta_scale <- fit_out$theta
  int_pos <- c(1, (p1 + 2), (1 + p1 + 1 + p2 + 1)) # positions of intercepts
  zero_pos <- which(abs(theta_scale) < zero_tol)
  zero_pos_final <- zero_pos[!(zero_pos %in% int_pos)]

  theta <- theta_scale / x_sd_theta
  theta[zero_pos_final] <- 0 # set to zero
  names(theta) <- names_coef

  # Get standard errors ----
  if (family == "old_normal") {
    info_list_normal <- basic_information_matrices(
      theta = theta_scale[-length(theta_scale)], # remove nu_0
      x1 = x1,
      x2 = x2,
      y = y,
      lambda_beta = lambda,
      lambda_alpha = lambda,
      epsilon = epsilon_T
    )
    see_scale_normal <- c(get_see_normal(info_list_normal), NA) # NA for nu_0
    see <- see_scale_normal / x_sd_theta
    see[zero_pos_final] <- 0 # set to zero
    names(see) <- names_coef
  } else {
    see_scale <- suppressWarnings({
      get_see_func_user(
        theta = theta_scale,
        x1 = x1,
        x2 = x2,
        x3 = x3,
        y = y,
        tau = tau,
        epsilon = epsilon_T,
        kappa_omega = kappa_omega,
        list_family = list_families$smoothgnd,
        lambda = lambda
      )
    }) # scaled data

    see <- see_scale / x_sd_theta
    see[zero_pos_final] <- 0 # set to zero
    names(see) <- names_coef
  }

  kappa_val <- unname(nu_to_kappa(theta[length(theta)],
    kappa_omega = kappa_omega
  ))

  # Return ----
  out <- list(
    "coefficients" = theta,
    "see" = see,
    "family" = family,
    "model" = model,
    "plike" = fit_out$plike_val,
    "kappa" = kappa_val,
    "tau" = tau,
    "kappa_omega" = kappa_omega,
    "telescope_df" = fit_mat,
    "x" = x, # return supplied data -> matches the formula
    "y" = y,
    "call" = cl
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
#'   data = sniffer,
#'   family = "normal",
#'   model = "mpr"
#' )
#' summary(results)
#' @importFrom stats pnorm
#'
#' @export
summary.smoothic <- function(object, ...) {
  family <- object$family

  coefficients <- object$coefficients
  see <- object$see

  zeropos <- which(coefficients == 0)

  coefficients[zeropos] <- NA
  see[zeropos] <- NA

  zval <- coefficients / see

  if (family == "old_normal") {
    pval <- 1.96 * pnorm(abs(zval), lower.tail = FALSE)
  } else {
    kappa <- object$kappa
    kappa_omega <- object$kappa_omega
    tau <- object$tau

    # qgndnorm(p = 0.975, mu = 0, s = sqrt(2), kappa = 2, tau = 1e-6)
    times <- qgndnorm(
      p = 0.975, # like getting 1.96 from qnorm(0.975, lower.tail = FALSE)
      mu = 0,
      s = 1,
      kappa = kappa,
      tau = tau
    )

    pval <- times * pgndnorm_vec(
      q = abs(zval),
      mu = 0,
      s = 1,
      kappa = kappa,
      tau = tau
    )
  }

  coefmat <- cbind(
    Estimate = coefficients,
    SE = see,
    Z = zval,
    Pvalue = pval
  )

  out <- list(
    family = object$family,
    model = object$model,
    coefmat = coefmat,
    plike = unname(object$plike),
    kappa = object$kappa,
    call = object$call
  )
  class(out) <- "summary.smoothic"
  out
}

# print.smoothic ----------------------------------------------------------
#' @aliases print.smoothic
#' @export
print.smoothic <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Family:\n")
  print(x$family)
  cat("Model:\n")
  print(x$model)

  kappa_lgl <- ifelse(x$family == "sgnd", "yes_kappa", "no_kappa")

  switch(kappa_lgl,
    yes_kappa = {
      cat("\nCoefficients:\n")
      print(x$coefficients)
      cat("\nKappa:\n")
      print(x$kappa)
    },
    no_kappa = {
      coef_short <- x$coefficients[-length(x$coefficients)] # remove nu_0
      cat("\nCoefficients:\n")
      print(coef_short)
    }
  )
}

# print.summary.smoothic --------------------------------------------------
#' @aliases print.summary.smoothic
#' @importFrom stats symnum
#' @export
print.summary.smoothic <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Family:\n")
  print(x$family)
  cat("Model:\n")
  print(x$model)

  kappa_lgl <- ifelse(x$family == "sgnd", "yes_kappa", "no_kappa")

  # split coefmat into location, scale and shape
  names_coefmat <- row.names(x$coefmat)
  loc_ind <- grep("_beta", names_coefmat)
  scale_ind <- grep("alpha", names_coefmat)
  shape_ind <- grep("nu_0", names_coefmat)

  coefmat_loc <- x$coefmat[loc_ind, ]

  if (x$model == "mpr") {
    coefmat_scale <- x$coefmat[scale_ind, ]
  } else if (x$model == "spr") {
    coefmat_scale <- matrix(x$coefmat[shape_ind, ], nrow = 1) # if scale constant returns numeric, so need to convert to matrix
    rownames(coefmat_scale) <- "intercept_0_alpha"
    colnames(coefmat_scale) <- colnames(coefmat_loc)
  }

  coefmat_shape <- matrix(x$coefmat[shape_ind, ], nrow = 1)
  rownames(coefmat_shape) <- "intercept_0_nu   "
  colnames(coefmat_shape) <- colnames(coefmat_loc)


  switch(kappa_lgl, # if shape is estimated, then include nu, kappa in the output
    yes_kappa = {
      cat("\nCoefficients:\n")
      cat("\n")
      cat("Location:\n")
      printCoefmat_MON(coefmat_loc, # print the location
        cs.ind = 1:2,
        tst.ind = 3,
        P.values = TRUE,
        has.Pvalue = TRUE,
        signif.legend = FALSE,
        na.print = "0" # change NA to 0 for printing
      )

      cat("\n")
      cat("Scale:\n")
      printCoefmat_MON(coefmat_scale, # print the scale
        cs.ind = 1:2,
        tst.ind = 3,
        P.values = TRUE,
        has.Pvalue = TRUE,
        signif.legend = FALSE,
        na.print = "0" # change NA to 0 for printing
      )

      cat("\n")
      cat("Shape:\n")
      printCoefmat_MON(coefmat_shape, # print the shape
        cs.ind = 1:2,
        tst.ind = 3,
        P.values = TRUE,
        has.Pvalue = TRUE,
        signif.legend = FALSE,
        na.print = "0" # change NA to 0 for printing
      )

      # print signif stars
      cat("\n")
      printCoefmat_MON_stars(x$coefmat, # overall legend, just the stars
        cs.ind = 1:2,
        tst.ind = 3,
        P.values = TRUE,
        has.Pvalue = TRUE,
        signif.legend = TRUE,
        na.print = "0" # change NA to 0 for printing
      )
    },
    no_kappa = {
      cat("\nCoefficients:\n")
      cat("\n")
      cat("Location:\n")
      printCoefmat_MON(coefmat_loc,
        cs.ind = 1:2,
        tst.ind = 3,
        P.values = TRUE,
        has.Pvalue = TRUE,
        signif.legend = FALSE,
        na.print = "0" # change NA to 0 for printing
      )

      cat("\n")
      cat("Scale:\n")
      printCoefmat_MON(coefmat_scale,
        cs.ind = 1:2,
        tst.ind = 3,
        P.values = TRUE,
        has.Pvalue = TRUE,
        signif.legend = FALSE,
        na.print = "0" # change NA to 0 for printing
      )
      cat("\n")
      # print signif stars
      printCoefmat_MON_stars(x$coefmat[-nrow(x$coefmat), ], # overall legend, remove nu_0
        cs.ind = 1:2,
        tst.ind = 3,
        P.values = TRUE,
        has.Pvalue = TRUE,
        signif.legend = TRUE,
        na.print = "0" # change NA to 0 for printing
      )
    }
  )
  cat("\n")
  if (kappa_lgl == "yes_kappa") {
    cat("Kappa Estimate:\n")
    print(x$kappa)
  }

  cat("Penalized Likelihood:\n")
  print(x$plike) # BIC or AIC = -2*plike

  cat("IC Value:\n")
  print(-2 * x$plike) # BIC or AIC = -2*plike
}

# predict.smoothic --------------------------------------------------------
#' @title Predict smoothic
#'
#' @aliases predict.smoothic
#'
#' @description \code{predict} method class \dQuote{\code{smoothic}}
#'
#' @param object an object of class \dQuote{\code{smoothic}} which is the result
#' of a call to \code{\link{smoothic}}.
#' @param newdata new data object
#' @param ... further arguments passed to or from other methods.
#'
#' @return a matrix containing the predicted values for the location mu and scale s
#'
#' @author Meadhbh O'Neill
#'
#' @examples
#' # Sniffer Data --------------------
#' # MPR Model ----
#' results <- smoothic(
#'   formula = y ~ .,
#'   data = sniffer,
#'   family = "normal",
#'   model = "mpr"
#' )
#' predict(results)
#'
#' @importFrom stats as.formula model.matrix
#'
#' @export
predict.smoothic <- function(object,
                             newdata,
                             ...) {
  if (!inherits(object, "smoothic")) {
    stop("Error: object should have 'smoothic' class")
  }

  if (missing(newdata)) {
    new_x <- object$x # no intercept (matches the formula from the original fit)
  }

  colnames_x <- colnames(object$x)

  if (!missing(newdata)) {
    if (!(class(newdata)[1] %in% c("matrix", "data.frame"))) {
      stop("Error: newdata should be a matrix or data frame")
    }

    newdata_x <- model.matrix(as.formula(paste0("~ ", as.character(object$call$formula)[[3]])), data = newdata)[, -1] # remove column of 1s

    paste0("~ ", as.character(object$call$formula)[[3]])

    colnames_newdata_x <- colnames(newdata_x)

    if (!identical(colnames_x, colnames_newdata_x)) {
      stop("Error: newdata should have same column names as the original supplied data")
    }

    new_x <- newdata_x
  }

  new_x_mat <- cbind(1, as.matrix(new_x)) # include intercept

  # Estimates mpr or spr
  p <- ncol(object$x)

  if (object$model == "mpr") {
    theta_final <- object$coefficients
  } else if (object$model == "spr") {
    theta_final <- c(object$coefficients[1:(p + 2)], rep(0, p), object$coefficients[p + 3])
  }

  theta_beta <- theta_final[1:(p + 1)]
  theta_alpha <- theta_final[(p + 2):((2 * p) + 2)]

  # predict mu
  location_pred <- new_x_mat %*% theta_beta

  # predict scale
  scale_pred <- sqrt(exp(new_x_mat %*% theta_alpha))

  # combine
  out_mat <- cbind(
    location_pred,
    scale_pred
  )
  colnames(out_mat) <- c("mu", "s")

  out_mat
}

# plot_effects ------------------------------------------------------------
#' @title Plot conditional density curves
#'
#' @description This function plots the model-based conditional density curves for
#' different effect combinations. For example, take a particular covariate that is selected
#' in the final model. The other selected covariates are fixed at their median values by default
#' (see \code{covariate_fix} to fix at other values) and then the plotted red and blue densities
#' correspond to the modification of the chosen covariate as \dQuote{low} (25th quantile by default) and
#' \dQuote{high} (75th quantile by default).
#'
#' @param obj An object of class \dQuote{\code{smoothic}} which is the result
#' of a call to \code{\link{smoothic}}.
#' @param what The covariate effects to be plotted, default is \code{what = "all"}. The user
#' may supply a vector of covariate names to be plotted (only covariates selected in the final
#' model can be plotted).
#' @param show_average_indiv Should a \dQuote{baseline} or \dQuote{average} individual be shown,
#' default is \code{show_average_indiv = TRUE}. If \code{show_average_indiv = FALSE} then this
#' is not shown.
#' @param p The probabilities given to the \code{\link{quantile}} function. This corresponds to the plotted
#' red and blue density curves where the chosen covariate is modified as \dQuote{low} and \dQuote{high}.
#' The default is \code{p = c(0.25, 0.75)} to show the 25th and 75th quantiles.
#' @param covariate_fix Optional values to fix the covariates at that are chosen in the final model. When
#' not supplied, the covariates are fixed at their median values. See the example for more detail.
#' @param density_range Optional range for which the density curves should be plotted.
#'
#' @return A plot of the conditional density curves.
#'
#' @author Meadhbh O'Neill
#'
#' @examples
#' # Sniffer Data --------------------
#' # MPR Model ----
#' results <- smoothic(
#'   formula = y ~ .,
#'   data = sniffer,
#'   family = "normal",
#'   model = "mpr"
#' )
#' plot_effects(results)
#'
#' # Only plot gastemp and gaspres
#' # Do not show the average individual plot
#' # Plot the lower and upper density curves using 10th quantile (lower) and 90th quantile (upper)
#' # Fix violent to its violent to 820 and funding to 40
#'
#' plot_effects(results,
#'              what = c("gastemp", "gaspres"),
#'              show_average_indiv = FALSE,
#'              p = c(0.1, 0.9),
#'              covariate_fix = c("gastemp" = 70,
#'                                "gaspres" = 4))
#'
#' # The curves for the gastemp variable are computed by fixing gaspres = 4 (as is specified
#' # in the input). The remaining variables that are not specified in covariate_fix are fixed
#' # to their median values (i.e., tanktemp is fixed at its median). gastemp is then modified
#' # to be low (10th quantile) and high (90th quantile), as specified by p in the function.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import tibble
#' @importFrom rlang .data
#' @importFrom data.table rbindlist
#'
#' @export
plot_effects <- function(obj,
                         what = "all",
                         show_average_indiv = TRUE, # show the average individual?
                         p = c(0.25, 0.75), # quantiles to be plotted
                         covariate_fix, # specific values of the covariate densities (named vector)
                         density_range) {
  quantile_values <- p

  if (length(quantile_values) != 2) {
    stop("Error: p must be length 2 (lower and upper quantile values)")
  }

  if (any(quantile_values > 1) | any(quantile_values < 0)) {
    stop("Error: values of p must be between 0 and 1")
  }

  if (quantile_values[1] > quantile_values[2]) {
    stop("Error: the first element of p is the should be less than the second element of p")
  }


  fit_obj <- obj
  n <- length(fit_obj$y)
  p <- ncol(fit_obj$x) # not including intercept

  # Estimates mpr or spr
  if (fit_obj$model == "mpr") {
    theta_final <- fit_obj$coefficients
  } else if (fit_obj$model == "spr") {
    theta_final <- c(fit_obj$coefficients[1:(p + 2)], rep(0, p), fit_obj$coefficients[p + 3])
  }

  if (is.null(colnames(fit_obj$x))) {
    colnames_x <- paste0("X_", 1:p)
    colnames(fit_obj$x) <- colnames_x
  }

  names_coef_labels <- colnames(fit_obj$x)
  names_coef_labels_theta <- c(
    "inter", names_coef_labels,
    "inter", names_coef_labels,
    "inter"
  )

  # Estimates (theta)
  names_coef_keep <- theta_final %>%
    as_tibble(.name_repair = "minimal") %>%
    add_column(
      coeff = names_coef_labels_theta,
      .before = 1
    ) %>%
    filter(
      .data$coeff != "inter",
      .data$value != 0
    ) %>%
    pull(.data$coeff) %>%
    unique()

  # Response value
  response_val <- as.character(fit_obj$call$formula[[2]])

  # Build dataset with response included (no intercept columns)
  dataset_raw <- cbind(fit_obj$x, fit_obj$y) %>%
    as_tibble(.name_repair = "minimal")
  colnames(dataset_raw)[p + 1] <- response_val

  # Calculate
  levels_summary <- c("average", "lower", "upper")
  levels_quant <- c("lower", "upper")
  names_coef_labels <- colnames(dataset_raw)[-(p + 1)]

  names_coef_response_keep <- c("response", names_coef_keep)
  levels_response_quant <- c("response", levels_quant)

  # What is to be plotted
  if ("all" %in% what) {
    coef_plot_names <- names_coef_keep
  } else if (any(!(what %in% names_coef_keep))) {
    stop("Error: please choose a variable that is selected in the model")
  } else {
    coef_plot_names <- what
  }

  if (show_average_indiv == TRUE) {
    coef_plot_names <- c("response", coef_plot_names)
  }

  # Any binary
  bin_unique <- as.data.frame(fit_obj$x) %>%
    map_dbl(~ length(unique(.x)))

  bin_pos <- which(bin_unique == 2)

  # Values to fix the covariates to
  if (missing(covariate_fix)) {
    fix_row <- dataset_raw %>%
      dplyr::select(-all_of(response_val)) %>%
      reframe(across(everything(), ~ median(.))) %>% # use median values if vector of values not supplied
      add_column(typee = factor("fix_value"), .before = 1)
  } else if (!missing(covariate_fix)) {
    if (!IsNamedVector(covariate_fix)) { # vector should have covariate names
      stop("Error: covariate_fix should be a vector with variable names and numeric values")
    }

    covariate_fix_names <- names(covariate_fix)
    # make sure supplied names are actual variable names
    if (any(!(covariate_fix_names %in% names_coef_keep))) {
      stop("Error: only the variables selected in the final model should appear in covariate_fix")
    }

    # make sure fix values are within the range of each variable
    df_range <- dataset_raw %>%
      select(all_of(covariate_fix_names)) %>%
      reframe(across(everything(), ~ c(min(.), max(.))))

    range_lgl <- map2_lgl(.x = covariate_fix,
                      .y = df_range, ~ {
                      lgl_val <- (.x >= .y[1]) & (.x <= .y[2])
                      lgl_val
                      })
    if (any(range_lgl == FALSE)) {
      stop("Error: the values of covariate_fix should be within the range of each variable")
    }

    # make the row of fixed values
    fix_row <- dataset_raw %>%
      dplyr::select(-all_of(response_val)) %>%
      reframe(across(everything(), ~ median(.))) # use median values if vector of values not supplied

    fix_row_vec <- fix_row %>%
      unlist()

    # if values supplied, then insert them
    for (i in seq_along(covariate_fix)) {
      name_now <- names(covariate_fix[i])

      pos_now <- which(names(fix_row_vec) == name_now)

      fix_row_vec[pos_now] <- covariate_fix[i]
    }

    fix_row <- as.matrix(fix_row)
    fix_row[1,] <- fix_row_vec
    fix_row <- as_tibble(as.data.frame(fix_row)) %>%
      add_column(typee = factor("fix_value"), .before = 1) # tibble output
  }

  df_summary_1 <- dataset_raw %>%
    reframe(across(everything(), ~ c(
      median(.),
      quantile(., probs = quantile_values) # input p
    ))) %>%
    add_column(
      typee = factor(levels_summary, levels = levels_summary),
      .before = 1
    ) %>%
    dplyr::select(-all_of(response_val))

  df_summary <- bind_rows(df_summary_1,
                          fix_row)

  # Change binary
  for (i in bin_pos) {
    df_summary[, i + 1] <- c(
      Mode_calc(unlist(dataset_raw[, i])),
      min(unlist(dataset_raw[, i])),
      max(unlist(dataset_raw[, i])),
      Mode_calc(unlist(dataset_raw[, i])) # fixed value
    )
  }

  df_summary_t <- df_summary %>%
    pivot_longer(-.data$typee,
      names_to = "coeff"
    ) %>%
    pivot_wider(
      names_from = .data$typee,
      values_from = .data$value
    ) %>%
    mutate(coeff = factor(.data$coeff, levels = names_coef_labels))

  df_average_indiv <- df_summary %>%
    filter(.data$typee == "average") %>%
    add_column(
      inter = 1,
      .before = 2
    ) %>%
    add_column(
      coeff = "response",
      .before = 1
    ) %>%
    dplyr::select(-.data$typee) %>%
    add_column(
      typee = "response",
      .after = 1
    )

  df_dist_est_prep <- names_coef_keep %>%
    map(~ {
      coef_now <- .x

      baseline <- df_summary_t %>%
        filter(.data$coeff == coef_now) %>%
        dplyr::select(.data$lower, .data$upper) %>%
        unlist()

      fix_vec <- df_summary %>%
        filter(.data$typee == "fix_value") %>% # values that are specified, or else median if not specified
        dplyr::select(-.data$typee)

      new_df <- rbind(
        fix_vec,
        fix_vec
      )
      new_df[, which(colnames(new_df) == coef_now)] <- baseline
      new_df %>%
        add_column(
          coeff = coef_now,
          typee = levels_quant,
          .before = 1
        )
    }) %>%
    data.table::rbindlist()

  df_dist_est_1 <- df_dist_est_prep %>%
    add_column(
      "inter" = 1,
      .before = 3
    )

  df_dist_est <- bind_rows(
    df_average_indiv,
    df_dist_est_1
  ) %>%
    mutate(
      coeff = factor(.data$coeff, levels = names_coef_response_keep),
      typee = factor(.data$typee, levels = levels_response_quant)
    )
  # return dist_est in final object so can make sure everything ok
  df_dist_est_out <- df_dist_est %>%
    dplyr::select(.data$coeff, .data$typee, all_of(names_coef_keep)) %>%
    dplyr::filter(.data$coeff %in% coef_plot_names) %>%
    mutate(coeff = recode_factor(.data$coeff, "response" = "average"),
           typee = recode_factor(.data$typee, "response" = "average"))

  df_dist_est_mat <- df_dist_est %>%
    dplyr::select(-c(.data$coeff, .data$typee)) %>%
    as.matrix()

  # Beta and Alpha
  theta_final_beta <- theta_final[1:(p + 1)]
  theta_final_alpha <- theta_final[(p + 2):((2 * p) + 2)] # will be zeros if spr
  theta_final_nu <- theta_final[(2 * p) + 3]

  mu_i_vec <- as.numeric(df_dist_est_mat %*% theta_final_beta)
  s_i_vec <- as.numeric(sqrt(exp(df_dist_est_mat %*% theta_final_alpha)))
  kappa_const <- unname(nu_to_kappa(theta_final_nu, fit_obj$kappa_omega)) # constant nu

  names_dist_est_1 <- df_dist_est %>%
    dplyr::select(.data$coeff, .data$typee)

  names_dist_est <- split(names_dist_est_1, seq(1:nrow(names_dist_est_1))) %>%
    map(unlist)

  mu_s_list <- map2(
    .x = mu_i_vec,
    .y = s_i_vec,
    ~ c(.x, .y)
  )

  # Range of density plot
  if (!missing(density_range)) {
    x_seq <- seq(density_range[1], density_range[2], length.out = 200)
  } else {
    x_seq <- seq(range(fit_obj$y)[1], range(fit_obj$y)[2], length.out = 200)
  }

  df_sgnd_1 <- mu_s_list %>%
    map(~ {
      tibble(
        x = x_seq,
        y = dgndnorm(x = x_seq, mu = .x[1], s = .x[2], kappa = kappa_const, tau = fit_obj$tau)
      )
    })
  df_sgnd <- map2(.x = df_sgnd_1, .y = names_dist_est, ~ {
    .x %>%
      add_column(
        coeff = .y[1],
        typee = .y[2],
        .before = 1
      )
  }) %>%
    data.table::rbindlist() %>%
    mutate(
      coeff = factor(.data$coeff, levels = names_coef_response_keep),
      typee = factor(.data$typee, levels = levels_response_quant)
    )

  labels_facets <- c("average", names_coef_keep)
  names(labels_facets) <- names_coef_response_keep

  names_coef_input <- c("average", names_coef_keep)

  # label quantile values according to p input vector (lower, upper)
  levels_resquant_p <- c("response", quantile_values)
  names(levels_resquant_p) <- c("response", levels_quant)

  levels_quant_p <- levels_resquant_p[2:3]

  col_pal_typee <- c(
    "black",
    "#E41A1C", # red
    "#0072B2"
  ) # blue
  names(col_pal_typee) <- levels_resquant_p

  fill_pal_typee <- c(
    "NA",
    "#E41A1C", # red
    "#0072B2"
  ) # blue
  names(fill_pal_typee) <- levels_resquant_p


  # Remove tail values less that 1e-5
  df_sgnd_rough <- df_sgnd %>%
    mutate(y_rough = ifelse(.data$y < 1e-5, NA, .data$y)) %>%
    filter(.data$coeff %in% coef_plot_names) %>%
    mutate(coeff = factor(.data$coeff, levels = coef_plot_names)) %>%
    mutate(typee = recode_factor(.data$typee, !!!levels_resquant_p))


  fig_effects <- df_sgnd_rough %>%
    filter(!is.na(.data$y_rough)) %>%
    ggplot(aes(
      x = .data$x,
      y = .data$y_rough,
      colour = .data$typee,
      fill = .data$typee
    )) +
    facet_wrap(~ .data$coeff,
      ncol = 1,
      labeller = as_labeller(labels_facets)
    ) +
    geom_line() +
    geom_ribbon(aes(ymin = 0, ymax = .data$y_rough),
      alpha = 0.25
    ) +
    scale_colour_manual(values = col_pal_typee, breaks = levels_quant_p, name = "Quantile") +
    scale_fill_manual(values = fill_pal_typee, breaks = levels_quant_p, name = "Quantile") +
    theme_bw() +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    theme(
      # legend.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    labs(x = response_val)
  fig_effects$covariate_vals <- df_dist_est_out
  fig_effects
}


# plot_paths --------------------------------------------------------------
#' @title Plot the \eqn{\epsilon}-telescope coefficient paths
#'
#' @description This function plots the standardized coefficient values with respect
#' to the \eqn{\epsilon}-telescope for the location (and dispersion) components.
#'
#' @param obj An object of class \dQuote{\code{smoothic}} which is the result
#' of a call to \code{\link{smoothic}}.
#' @param log_scale_x Default is \code{log_scale_x = TRUE}, which uses a log scale
#' on the x-axis. If \code{log_scale_x = FALSE}, then the raw values of the \eqn{\epsilon}-telescope are plotted.
#' @param log_scale_x_pretty Default is \code{log_scale_x_pretty = TRUE}, where the x-axis labels are \dQuote{pretty}.
#' \code{epsilon_1} and \code{epsilon_T} must be a number to the power of 10 for this to apply.
#' @param facet_scales Default is \code{facet_scales = "fixed"}. This is supplied to \code{\link{facet_wrap}}.
#'
#' @return A plot of the standardized coefficient values through the \eqn{\epsilon}-telescope.
#'
#' @author Meadhbh O'Neill
#'
#' @examples
#' # Sniffer Data --------------------
#' # MPR Model ----
#' results <- smoothic(
#'   formula = y ~ .,
#'   data = sniffer,
#'   family = "normal",
#'   model = "mpr"
#' )
#' plot_paths(results)
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom purrr map_dbl
#' @importFrom rlang .data
#'
#' @export
plot_paths <- function(obj,
                       log_scale_x = TRUE,
                       log_scale_x_pretty = TRUE,
                       facet_scales = "fixed") {
  if (log_scale_x == FALSE & log_scale_x_pretty == TRUE) {
    stop("Error: if log_scale_x = FALSE, then set log_scale_x_pretty = FALSE")
  }

  fit_obj <- obj
  telescope_df <- fit_obj$telescope_df

  if (nrow(telescope_df) < 10) {
    warning("Ensure an adequate number of steps_T are used")
  }

  names_coef_fit <- names(fit_obj$coefficients)
  names_coef <- names_coef_fit[!names_coef_fit %in% c(
    "intercept_0_beta",
    "intercept_0_alpha",
    "alpha_0",
    "nu_0"
  )]
  plot_df_prep <- telescope_df %>%
    dplyr::select(
      .data$epsilon,
      contains(c("beta", "alpha")),
      -c(
        "beta_0",
        "alpha_0"
      )
    ) %>%
    rename_all(~ c("epsilon", names_coef)) %>%
    pivot_longer(-.data$epsilon) %>%
    mutate(type = case_when(
      grepl("_beta", .data$name) ~ "Location",
      grepl("_alpha", .data$name) ~ "Scale"
    )) %>%
    mutate(coeff = sub("_.*", "", .data$name)) # extract variable name

  if (log_scale_x == TRUE) {
    plot_df <- plot_df_prep %>%
      mutate(epsilon_plot = log(.data$epsilon))
    x_label <- "log(epsilon)"
  } else {
    plot_df <- plot_df_prep %>%
      mutate(epsilon_plot = .data$epsilon)
    x_label <- "epsilon"
  }

  if (log_scale_x == TRUE & log_scale_x_pretty == TRUE) {
    x_label <- "epsilon"
  }

  fig_paths <- plot_df %>%
    ggplot(aes(
      x = .data$epsilon_plot,
      y = .data$value,
      colour = .data$coeff
    )) +
    facet_wrap(~ .data$type,
      scales = facet_scales
    ) +
    geom_line() +
    labs(
      y = "Standardized Coefficient Value",
      x = x_label
    ) +
    guides(colour = guide_legend("Variable")) +
    theme_bw()

  if (log_scale_x == TRUE & log_scale_x_pretty == TRUE) {
    range_eps <- range(telescope_df$epsilon)
    range_eps_char <- changeSciNot(range_eps)

    suitable_vec <- sub("\\*.*", "", range_eps_char)

    if (any(suitable_vec != "1") == TRUE) {
      stop("Error: if log_scale_x_pretty = TRUE, then epsilon_1 and epsilon_T must be a power of 10")
    }

    range_eps_num <- sub(".*\\^", "", range_eps_char)
    eps_num <- seq(as.numeric(range_eps_num[1]), as.numeric(range_eps_num[2]))
    eps_choice <- as.numeric(paste0("1e", eps_num))
    eps_choice_char <- paste0("10^{", eps_num, "}")

    eps_tele_vec <- telescope_df$epsilon
    eps_vals <- eps_tele_vec[eps_choice %>%
      purrr::map_dbl(~ {
        which(abs(eps_tele_vec - .x) == min(abs(eps_tele_vec - .x)))
      })]
    log_eps_vals <- log(eps_vals)
    labels_eps <- paste0("c(", toString((paste0("'", log_eps_vals, "'=expression(", eps_choice_char, ")"))), ")")
    labels_eps <- eval(parse(text = labels_eps))

    fig_paths <- fig_paths +
      scale_x_continuous(
        labels = labels_eps,
        breaks = log_eps_vals,
        expand = expansion(mult = c(0.02, 0.02))
      )
  } else {
    fig_paths <- fig_paths +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)))
  }

  fig_paths
}

# Check that vector has names for plot_effects function
IsNamedVector <- function(VECTOR) {
  is.vector(VECTOR) & is.numeric(VECTOR) & !is.null(names(VECTOR)) & !any(is.na(names(VECTOR)))
  }

# Calculate mode
Mode_calc <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Scientific notation for epsilon plotting pretty
changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) # Transforms the number into scientific notation even if small
  output <- sub("e", "*10^", output) # Replace e with 10^
  output <- sub("\\+0?", "", output) # Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) # Leaves - symbol but removes leading zeros on expoent, if < 1
  output
}

# Fitting function for package --------------------------------------------
fitting_func_pkg <- function(x1,
                             x2,
                             x3,
                             y,
                             family,
                             optimizer,
                             lambda,
                             epsilon_1,
                             epsilon_T,
                             steps_T,
                             max_it_vec,
                             kappa,
                             tau,
                             fix_kappa_lgl,
                             kappa_omega,
                             # nlm
                             stepmax_nlm,
                             method_c_tilde = "integrate",
                             # manual
                             tol = 1e-8,
                             initial_step = 10,
                             max_step_it = 1e3,
                             method_c_tilde_deriv = "finite_diff",
                             algorithm = "cross_zero",
                             h_kappa = 1e-5,
                             # initial values
                             method_initial_values,
                             ...) {
  list_family <- list_families$smoothgnd

  if (missing(method_initial_values)) {
    switch(family,
      "sgnd" = {
        method_initial_values <- "lm"
      },
      "normal" = {
        method_initial_values <- "lm"
      },
      "laplace" = {
        method_initial_values <- "rlm" # robust linear model
      }
    )
  }

  # Initial Values ----
  theta_init <- initial_values_func(
    x1 = x1,
    x2 = x2,
    x3 = x3, # col of 1s
    y = y,
    method_initial_values = method_initial_values, # either "lm" or "rlm
    kappa_omega = kappa_omega
  )

  # Is kappa fixed? ----
  pos_nu <- length(theta_init)
  if (fix_kappa_lgl) { # if true, then need to replace nu in theta_init with fixed value
    theta_init[pos_nu] <- kappa_to_nu(kappa, kappa_omega)
  }

  # Optimizer choice ----
  if (optimizer == "manual") {
    fit_opt <- "fullopt"
    if (fix_kappa_lgl) { # if fixed nu
      fit_opt <- "fullopt_fixed_nu"
    }
  } else if (optimizer == "nlm") {
    fit_opt <- "nlm"
    if (fix_kappa_lgl) {
      fit_opt <- "nlm_fixed_nu"
    }
  }

  # Fitting function ----
  fit_out <- fitting_func_base(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    optimizer = fit_opt,
    stepmax_nlm = stepmax_nlm,
    theta_init = theta_init, # initial values for scaled data, should contain fixed value for nu if taking that option
    tau_1 = tau,
    tau_T = tau,
    epsilon_1 = epsilon_1,
    epsilon_T = epsilon_T,
    steps_T = steps_T,
    list_family = list_family,
    method_c_tilde = method_c_tilde,
    method_c_tilde_deriv = method_c_tilde_deriv,
    algorithm = algorithm,
    h_kappa = h_kappa,
    kappa_omega = kappa_omega,
    lambda_beta = lambda,
    lambda_alpha = lambda,
    lambda_nu = lambda,
    initial_step = initial_step,
    max_step_it = max_step_it,
    tol = tol,
    max_it_vec = max_it_vec
  )
  # extract_fit <- extract_theta_plike_val(fit_res = fit_out)
  # extract_fit
  fit_out
}


# SEE function package ----------------------------------------------------
get_see_func_user <- function(theta,
                              x1,
                              x2,
                              x3,
                              y,
                              tau,
                              epsilon,
                              kappa_omega,
                              list_family,
                              lambda,
                              method_c_tilde = "integrate",
                              ...) {
  info_mat_list <- information_matrices_numDeriv(
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
    lambda_beta = lambda,
    lambda_alpha = lambda,
    lambda_nu = lambda
  )
  see <- get_see_now(info_mat_list)
  see
}

get_see_normal <- function(list_of_2_info_matrices) {
  observed_information_penalized <- list_of_2_info_matrices$observed_information_penalized
  observed_information_unpenalized <- list_of_2_info_matrices$observed_information_unpenalized

  inverse_observed_information_penalized <- solve(observed_information_penalized)

  # Square root of the diagonal of the variance-covariance matrix
  see <- sqrt(diag(inverse_observed_information_penalized %*% observed_information_unpenalized %*% inverse_observed_information_penalized)) # sandwich formula

  see
}


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

  F2x(
    u = p,
    F = F
  )
}

dgndnorm <- function(x,
                     mu,
                     s,
                     kappa,
                     tau) {
  c_tilde_val <- c_tilde_integrate_smoothgnd(
    kappa = kappa,
    tau = tau
  )

  density_gnd_approx(
    x = x,
    mu_loc = mu,
    s_scale = s,
    kappa_shape = kappa,
    tau = tau,
    c_tilde = c_tilde_val
  )
}

pgndnorm <- function(q,
                     mu,
                     s,
                     kappa,
                     tau) {
  if (!is.na(q)) {
    cdf_df <- as.data.frame(F_cdf_full(
      from = -20, # cdf
      to = 20,
      by = 0.001,
      kappa = kappa,
      tau = tau,
      mu = mu,
      s = s
    ))

    # Find values of x in dataframe closest to q
    pos <- which(abs(q - cdf_df$x) == min(abs(q - cdf_df$x)))

    # lower.tail = FALSE (subtract from 1)
    1 - cdf_df[pos, ]$F
  } else {
    NA
  }
}

pgndnorm_vec <- Vectorize(pgndnorm,
  vectorize.args = "q"
)

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
    c_tilde_trapezoid = c_tilde_trapezoid_smoothgnd
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

  likelihood_neg(
    theta = theta,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y = y,
    tau = tau,
    list_general = list_general,
    method_c_tilde = method_c_tilde,
    kappa_omega = kappa_omega
  )
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
  Sigma_beta <- diag(
    c(
      0,
      (lambda_beta / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
    ),
    nrow = p1 + 1,
    ncol = p1 + 1
  )
  I_beta <- ((t(x1 * c(W_beta))) %*% x1) + Sigma_beta

  ## Alpha
  Sigma_alpha <- diag(
    c(
      0,
      (lambda_alpha / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (alpha1_p^2))) / ((alpha1_p^2 + epsilon^2)^3))
    ),
    nrow = p2 + 1,
    ncol = p2 + 1
  )
  I_alpha <- ((t(x2 * c(W_alpha))) %*% x2) + Sigma_alpha

  ## Nu
  Sigma_nu <- diag(
    c(
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
  Sigma_beta <- diag(
    c(
      0,
      (lambda_beta / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
    ),
    nrow = p1 + 1,
    ncol = p1 + 1
  )
  I_beta <- ((t(x1 * c(W_beta))) %*% x1) + Sigma_beta

  ## Alpha
  Sigma_alpha <- diag(
    c(
      0,
      (lambda_alpha / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (alpha1_p^2))) / ((alpha1_p^2 + epsilon^2)^3))
    ),
    nrow = p2 + 1,
    ncol = p2 + 1
  )
  I_alpha <- ((t(x2 * c(W_alpha))) %*% x2) + Sigma_alpha

  ## Nu
  Sigma_nu <- diag(
    c(
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

  observed_information_penalized <- tryCatch(
    -numDeriv::hessian(likelihood_penalized, # negative second deriv
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

  observed_information_unpenalized <- tryCatch(
    -numDeriv::hessian(likelihood, # negative second deriv, unpen
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
                      max_it_vec) {
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
    iter_val <- max_it_vec[i]

    pos <- i

    t_fit <- it_loop(
      theta = t_initial_guess,
      FUN = FUN_tele, # nr_penalty
      tol = tol,
      max_it = iter_val,
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
    iter_val <- iterlim[i]
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
            iterlim = iter_val,
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
            iterlim = iter_val
          )
        })
      }
    )

    t_initial_guess <- t_fit$estimate
    steps_full <- 0 # 1 if true, 0 if false (just set to zero for nlm as not dealing with step halving)

    it_full <- ifelse(t_fit$iterations == iter_val, 1, 0)

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

  switch(stepmax_lgl,
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
      res <- nlm(
        f = likelihood_neg_fixed, # likelihood_neg_fixed - nlm minimizes, unpenalized
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
      res <- nlm(
        f = likelihood_neg_fixed, # likelihood_neg_fixed - nlm minimizes, unpenalized
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
                                iterlim_nlm,
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
    iter_val <- iterlim_nlm[i]
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
        iterlim_nlm = iter_val,
        stepmax = stepmax
      )
    })

    t_initial_guess <- t_fit$estimate
    steps_full <- 0 # 1 if true, 0 if false (just set to zero for nlm as not dealing with step halving)

    it_full <- ifelse(t_fit$iterations == iter_val, 1, 0)

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
                              stepmax_nlm,
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
                              max_it_vec) {
  # tau_tele_vec <- rev(exp(seq(log(tau_T), log(tau_1), length = steps_T)))
  tau_tele_vec <- rep(tau_1, times = steps_T)
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
      max_it_vec = max_it_vec
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
      iterlim = max_it_vec,
      stepmax = stepmax_nlm
    ),
    "fullopt_fixed_nu" = telescope(
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
      max_it_vec = max_it_vec
    ),
    "nlm_fixed_nu" = telescope_nlm_fixed(
      theta_init = theta_init, # contains fixed value
      fix_indices = length(theta_init), # fix nu
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
      iterlim_nlm = max_it_vec,
      stepmax = stepmax_nlm
    )
  )
}

# Extract theta & plike_val -----------------------------------------------
extract_theta_plike_val <- function(fit_res) {
  # last row
  row_n <- fit_res[nrow(fit_res), ]
  plike_val <- row_n$plike_val

  colnames_now <- colnames(fit_res)

  # get positions of column names of coefficients
  pos_extract <- c(
    grep("beta_", colnames_now),
    grep("alpha_", colnames_now),
    grep("nu_0", colnames_now)
  )
  theta <- unlist(row_n[, pos_extract])

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
                                kappa_omega) {
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
    "rlm" = {
      rlm_fit <- rlm(y ~ x1[, -1])
      rlm_coef_sig <- c(
        unname(rlm_fit$coefficients),
        log((summary(rlm_fit)$sigma^2))
      )
      theta_init <- c(rlm_coef_sig, rep(0, p2), kappa_to_nu(
        kappa = 1,
        kappa_omega = kappa_omega
      )) # kappa = exp(nu_0) + omega, so then nu_0 = log(kappa - omega)
      theta_init
    }
  )
}

# get thetas from dataf ---------------------------------------------------
get_see_now <- function(info_list) {
  inv_pen <- solve(info_list$observed_information_penalized)

  sqrt(diag(inv_pen %*% info_list$observed_information_unpenalized %*% inv_pen))
}






# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# BASIC ORIGINAL CODING METHOD --------------------------------------------
# ** Normal Likelihood ----------------------------------------------------
basic_normallike <- function(theta,
                             x1,
                             x2,
                             y) {
  n <- length(y)

  p1 <- ncol(x1) - 1

  beta <- theta[1:(p1 + 1)]
  alpha <- theta[-c(1:(p1 + 1))]

  mu <- x1 %*% beta
  phi <- exp(x2 %*% alpha)

  like <- -(n / 2) * log(2 * pi) - (1 / 2) * (sum(x2 %*% alpha)) - (1 / 2) * (sum((1 / phi) * ((y - mu)^2)))
  like # maximize likelihood
}

# ** Penalized Likelihood -------------------------------------------------
basic_penalizedlike <- function(theta,
                                x1,
                                x2,
                                y,
                                lambda_beta,
                                lambda_alpha,
                                epsilon) {
  # maximise this
  n <- length(y)

  p1 <- ncol(x1) - 1

  beta <- theta[1:(p1 + 1)]
  alpha <- theta[-c(1:(p1 + 1))]

  lambda_beta <- (eval(parse(text = lambda_beta)))
  lambda_alpha <- (eval(parse(text = lambda_alpha)))

  # remove intercepts as won't be penalized
  beta1_p <- beta[-1]
  alpha1_p <- alpha[-1]

  # both sets of coefficients penalized
  p_like <- basic_normallike(theta, x1, x2, y) - ((lambda_beta / 2) * (1 + (sum((beta1_p^2) / (beta1_p^2 + epsilon^2))))) - ((lambda_alpha / 2) * (1 + (sum((alpha1_p^2) / (alpha1_p^2 + epsilon^2)))))
  p_like # divided by two for penalized setup, +1 for estimation of intercepts
}

# ** Negative Penalized Likelihood ----------------------------------------
basic_neg_penalizedlike <- function(...) {
  -basic_penalizedlike(...)
}

# ** NR Penalty -----------------------------------------------------------
basic_nr_normal_penalty <- function(theta,
                                    x1,
                                    x2,
                                    y,
                                    lambda_beta,
                                    lambda_alpha,
                                    epsilon,
                                    initial_step,
                                    max_step_it) {
  n <- length(y)

  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1

  beta <- theta[1:(p1 + 1)]
  alpha <- theta[-c(1:(p1 + 1))]

  beta1_p <- beta[-1] # remove intercepts as won't be penalized
  alpha1_p <- alpha[-1]

  mu <- x1 %*% beta
  phi <- exp(x2 %*% alpha)

  tx1 <- t(x1) # save transpose of x1 as an object
  tx2 <- t(x2) # save transpose of x1 as an object

  lambda_beta <- (eval(parse(text = lambda_beta)))
  lambda_alpha <- (eval(parse(text = lambda_alpha)))

  ## Score
  ## Beta
  z_beta <- ((1 / phi) * (y - mu))
  v_beta <- c(0, (lambda_beta / 2) * ((2 * beta1_p * (epsilon^2)) / ((beta1_p^2 + epsilon^2)^2)))
  score_beta <- (tx1 %*% z_beta) - v_beta

  ## Alpha
  z_alpha <- (-1 / 2) + (1 / (2 * phi)) * (((y - mu)^2))
  v_alpha <- c(0, (lambda_alpha / 2) * ((2 * alpha1_p * (epsilon^2)) / ((alpha1_p^2 + epsilon^2)^2)))
  score_alpha <- (tx2 %*% z_alpha) - v_alpha

  score <- c(score_beta, score_alpha)

  ## Information Matrix
  # Beta
  W_beta <- (1 / phi)
  Sigma_beta <- diag(
    c(
      0,
      (lambda_beta / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
    ),
    nrow = p1 + 1,
    ncol = p1 + 1
  )
  I_beta <- ((t(x1 * c(W_beta))) %*% x1) + Sigma_beta

  # Alpha
  W_alpha <- (1 / (2 * phi)) * ((y - mu)^2)
  Sigma_alpha <- diag(
    c(
      0,
      (lambda_alpha / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (alpha1_p^2))) / ((alpha1_p^2 + epsilon^2)^3))
    ),
    nrow = p2 + 1,
    ncol = p2 + 1
  )
  I_alpha <- ((t(x2 * c(W_alpha))) %*% x2) + Sigma_alpha

  # Cross deriv
  W_beta_alpha <- rep(0, n) # Can change depending on orthogonality of parameters
  # W_beta_alpha <- (1 / phi) * (y - mu) # true derivative
  I_beta_alpha <- t(x1 * c(W_beta_alpha)) %*% x2
  I_alpha_beta <- t(I_beta_alpha)

  information <- rbind(
    cbind(I_beta, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha)
  )

  delta <- solve(information, score)

  ## Step halving
  like_current <- basic_penalizedlike(theta, x1, x2, y, lambda_beta, lambda_alpha, epsilon)
  like_new <- like_current - 1
  j <- 0
  while (like_new < like_current & j < max_step_it) {
    theta_new <- theta + ((initial_step * delta) / (2^j))
    like_new <- basic_penalizedlike(theta_new, x1, x2, y, lambda_beta, lambda_alpha, epsilon)
    j <- j + 1
  }
  max_step_reached <- ifelse(j == max_step_it, TRUE, FALSE)
  return(list(
    "estimate" = theta_new, "maximum" = like_new,
    "steps" = c(max_step_reached, j)
  ))
}

# ** Information Matrices -------------------------------------------------
basic_information_matrices <- function(theta,
                                       x1,
                                       x2, y,
                                       lambda_beta,
                                       lambda_alpha,
                                       epsilon) {
  n <- length(y)

  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1


  beta <- theta[1:(p1 + 1)]
  alpha <- theta[-c(1:(p1 + 1))]

  beta1_p <- beta[-1] # remove intercepts as won't be penalized
  alpha1_p <- alpha[-1]

  mu <- x1 %*% beta
  phi <- exp(x2 %*% alpha)

  lambda_beta <- (eval(parse(text = lambda_beta)))
  lambda_alpha <- (eval(parse(text = lambda_alpha)))

  ## Penalized Observed Information Matrix
  # Beta
  W_beta <- (1 / phi)
  Sigma_beta <- diag(
    c(
      0,
      (lambda_beta / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (beta1_p^2))) / ((beta1_p^2 + epsilon^2)^3))
    ),
    nrow = p1 + 1,
    ncol = p1 + 1
  )
  I_beta <- ((t(x1 * c(W_beta))) %*% x1) + Sigma_beta

  # Alpha
  W_alpha <- (1 / (2 * phi)) * ((y - mu)^2)
  Sigma_alpha <- diag(
    c(
      0,
      (lambda_alpha / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (alpha1_p^2))) / ((alpha1_p^2 + epsilon^2)^3))
    ),
    nrow = p2 + 1,
    ncol = p2 + 1
  )
  I_alpha <- ((t(x2 * c(W_alpha))) %*% x2) + Sigma_alpha

  # Cross deriv
  # W_beta_alpha <- rep(0, n) # Can change depending on orthogonality of parameters
  W_beta_alpha <- (1 / phi) * (y - mu) #  negative derivative wrt to beta and alpha
  I_beta_alpha <- t(x1 * c(W_beta_alpha)) %*% x2
  I_alpha_beta <- t(I_beta_alpha)

  observed_information_penalized <- rbind(
    cbind(I_beta, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha)
  )


  ## Unpenalized Information evaluated at penalized estimates
  I_beta_unpen <- ((t(x1 * c(W_beta))) %*% x1) # no sigma
  I_alpha_unpen <- ((t(x2 * c(W_alpha))) %*% x2) # no sigma

  observed_information_unpenalized <- rbind(
    cbind(I_beta_unpen, I_beta_alpha),
    cbind(I_alpha_beta, I_alpha_unpen)
  )


  list(
    "observed_information_penalized" = observed_information_penalized,
    "observed_information_unpenalized" = observed_information_unpenalized
  )
}

# * basic_telescope -------------------------------------------------------------
basic_telescope <- function(x1,
                            x2,
                            y,
                            theta_init,
                            eps_tele_vec,
                            lambda_beta,
                            lambda_alpha,
                            initial_step,
                            max_step_it,
                            tol,
                            max_it_vec) {
  p1 <- ncol(x1) - 1 # important for getting dimensions of the output matrix
  p2 <- ncol(x2) - 1

  t_initial_guess <- theta_init
  t_res_mat <- matrix(NA,
    nrow = length(eps_tele_vec),
    ncol = (length(t_initial_guess) + 5)
  )
  colnames(t_res_mat) <- c(
    "epsilon",
    paste0("beta_", 0:p1),
    paste0("alpha_", 0:p2),
    "plike_val",
    "step_full",
    "it_full",
    "it"
  )
  ## basic_telescope
  steps_vec <- length(eps_tele_vec)

  for (i in 1:steps_vec) {
    eps_val <- eps_tele_vec[i]
    iter_val <- max_it_vec[i]

    pos <- i

    t_fit <- it_loop(
      theta = t_initial_guess,
      FUN = basic_nr_normal_penalty,
      tol = tol,
      max_it = iter_val,
      x1 = x1,
      x2 = x2,
      y = y,
      lambda_beta = lambda_beta,
      lambda_alpha = lambda_alpha,
      epsilon = eps_val,
      initial_step = initial_step,
      max_step_it = max_step_it
    )
    t_initial_guess <- t_fit$estimate
    steps_full <- ifelse(any(t_fit$max_step_reached == 1), 1, 0) # 1 if true, 0 if false

    t_res_mat[pos, ] <- c(i, t_fit$estimate, t_fit$maximum, steps_full, t_fit$iterations)
  }
  as.data.frame(t_res_mat) # return dataframe
}

basic_telescope_nlm <- function(x1,
                                x2,
                                y,
                                theta_init,
                                eps_tele_vec,
                                lambda_beta,
                                lambda_alpha,
                                max_it_vec,
                                stepmax_nlm) {
  p1 <- ncol(x1) - 1 # important for getting dimensions of the output matrix
  p2 <- ncol(x2) - 1

  t_initial_guess <- theta_init
  t_res_mat <- matrix(NA,
    nrow = length(eps_tele_vec),
    ncol = (length(t_initial_guess) + 5)
  )
  colnames(t_res_mat) <- c(
    "epsilon",
    paste0("beta_", 0:p1),
    paste0("alpha_", 0:p2),
    "plike_val",
    "step_full",
    "it_full",
    "it"
  )

  stepmax_lgl <- ifelse(is.na(stepmax_nlm), "no_stepmax", "yes_stepmax")


  steps_vec <- length(eps_tele_vec)
  for (i in 1:steps_vec) {
    eps_val <- eps_tele_vec[i]
    iter_val <- max_it_vec[i]

    pos <- i

    switch(stepmax_lgl,
      "yes_stepmax" = {
        t_fit <- suppressWarnings({
          nlm(basic_neg_penalizedlike, # nlm minimizes -> likelihood_penalized_neg
            p = t_initial_guess,
            x1 = x1,
            x2 = x2,
            y = y,
            lambda_beta = lambda_beta,
            lambda_alpha = lambda_alpha,
            epsilon = eps_val,
            iterlim = iter_val,
            stepmax = stepmax_nlm
          )
        })
      },
      "no_stepmax" = {
        t_fit <- suppressWarnings({
          nlm(basic_neg_penalizedlike,
            p = t_initial_guess,
            x1 = x1,
            x2 = x2,
            y = y,
            lambda_beta = lambda_beta,
            lambda_alpha = lambda_alpha,
            epsilon = eps_val,
            iterlim = iter_val
          )
        })
      }
    )
    t_initial_guess <- t_fit$estimate
    steps_full <- 0 # set to zero for nlm as not dealing with step halving
    it_full <- ifelse(t_fit$iterations == iter_val, 1, 0)

    t_res_mat[pos, ] <- c(i, t_fit$estimate, -t_fit$minimum, steps_full, it_full, t_fit$iterations)
  }
  as.data.frame(t_res_mat) # return data frame
}

fitting_func_normal <- function(x1,
                                x2,
                                y,
                                optimizer,
                                epsilon_1,
                                epsilon_T,
                                steps_T,
                                lambda_beta,
                                lambda_alpha,
                                max_it_vec,
                                stepmax_nlm,
                                initial_step = 10,
                                max_step_it = 1e3,
                                tol = 1e-8) {
  eps_tele_vec <- rev(exp(seq(log(epsilon_T), log(epsilon_1), length = steps_T)))

  p1 <- ncol(x1) - 1
  p2 <- ncol(x2) - 1

  lm_fit <- lm(y ~ x1[, -1])
  lm_coef <- lm_fit$coefficients
  lm_coef_sig <- c(
    unname(lm_coef),
    log((summary(lm_fit)$sigma)^2)
  )
  theta_init <- c(lm_coef_sig, rep(0, p2))

  switch(optimizer,
    "manual" = {
      fit_out <- basic_telescope(
        x1 = x1,
        x2 = x2,
        y = y,
        theta_init = theta_init,
        eps_tele_vec = eps_tele_vec,
        lambda_beta = lambda_beta,
        lambda_alpha = lambda_alpha,
        initial_step = initial_step,
        max_step_it = max_step_it,
        tol = tol,
        max_it_vec = max_it_vec
      )
    },
    "nlm" = {
      fit_out <- basic_telescope_nlm(
        x1 = x1,
        x2 = x2,
        y = y,
        theta_init = theta_init,
        eps_tele_vec = eps_tele_vec,
        lambda_beta = lambda_beta,
        lambda_alpha = lambda_alpha,
        max_it_vec = max_it_vec,
        stepmax_nlm = stepmax_nlm
      )
    }
  )
  # extract_fit <- extract_theta_plike_val(fit_res = fit_out)
  # extract_fit
  fit_out
}

# My coefmat functions ----------------------------------------------------
printCoefmat_MON <- function(x, digits = max(3L, getOption("digits") - 2L), signif.stars = getOption("show.signif.stars"),
                             signif.legend = signif.stars, dig.tst = max(1L, min(
                               5L,
                               digits - 1L
                             )), cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(),
                             P.values = NULL, has.Pvalue = nc >= 4L && length(cn <- colnames(x)) &&
                               substr(cn[nc], 1L, 3L) %in% c("Pr(", "p-v"), eps.Pvalue = .Machine$double.eps,
                             na.print = "NA", quote = FALSE, right = TRUE, ...) {
  if (is.null(d <- dim(x)) || length(d) != 2L) {
    stop("'x' must be coefficient matrix/data frame")
  }
  nc <- d[2L]
  if (is.null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  } else if (P.values && !has.Pvalue) {
    stop("'P.values' is TRUE, but 'has.Pvalue' is not")
  }
  if (has.Pvalue && !P.values) {
    d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
    nc <- nc - 1
    has.Pvalue <- FALSE
  } else {
    xm <- data.matrix(x)
  }
  k <- nc - has.Pvalue - (if (missing(tst.ind)) {
    1
  } else {
    length(tst.ind)
  })
  if (!missing(cs.ind) && length(cs.ind) > k) {
    stop("wrong k / cs.ind")
  }
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  ok <- !(ina <- is.na(xm))
  for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
  if (length(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs !=
        0])) {
        floor(log10(range(acs[acs != 0], finite = TRUE)))
      } else {
        0
      }
      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits -
        digmin)), digits = digits)
    }
  }
  if (length(tst.ind)) {
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
      digits = digits
    )
  }
  if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc)))) {
    for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
  }
  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue) {
    ok[, -nc]
  } else {
    ok
  }
  x1 <- Cf[okP]
  dec <- getOption("OutDec")
  if (dec != ".") {
    x1 <- chartr(dec, ".", x1)
  }
  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
  if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(
      1L,
      digits - 1L
    ))
  }
  if (any(ina)) {
    Cf[ina] <- na.print
  }
  if (any(inan <- is.nan(xm))) {
    Cf[inan] <- "NaN"
  }
  if (P.values) {
    if (!is.logical(signif.stars) || is.na(signif.stars)) {
      warning("option \"show.signif.stars\" is invalid: assuming TRUE")
      signif.stars <- TRUE
    }
    if (any(okP <- ok[, nc])) {
      pv <- as.vector(xm[, nc])
      Cf[okP, nc] <- format.pval(pv[okP],
        digits = dig.tst,
        eps = eps.Pvalue
      )
      signif.stars <- signif.stars && any(pv[okP] < 0.1)
      if (signif.stars) {
        # MON
        max_pv <- max(pv, na.rm = TRUE)
        if (max_pv < 1) {
          max_cut <- 1
        } else {
          max_cut <- ceiling(max_pv * 2) / 2
        }

        cutpoints_vec <- c(0, 0.001, 0.01, 0.05, 0.1, max_cut)

        Signif <- stats::symnum(pv,
          corr = FALSE, na = FALSE,
          cutpoints = cutpoints_vec,
          symbols = c("***", "**", "*", ".", " ")
        )
        Cf <- cbind(Cf, format(Signif))
      }
    } else {
      signif.stars <- FALSE
    }
  } else {
    signif.stars <- FALSE
  }
  print.default(Cf,
    quote = quote, right = right, na.print = na.print,
    ...
  )
  if (signif.stars && signif.legend) {
    if ((w <- getOption("width")) < nchar(sleg <- attr(
      Signif,
      "legend"
    ))) {
      sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
    }
    cat("---\nSignif. codes:  ", sleg, sep = "", fill = w +
      4 + max(nchar(sleg, "bytes") - nchar(sleg)))
  }
  invisible(x)
}

printCoefmat_MON_stars <- function(x, digits = max(3L, getOption("digits") - 2L), signif.stars = getOption("show.signif.stars"),
                                   signif.legend = signif.stars, dig.tst = max(1L, min(
                                     5L,
                                     digits - 1L
                                   )), cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(),
                                   P.values = NULL, has.Pvalue = nc >= 4L && length(cn <- colnames(x)) &&
                                     substr(cn[nc], 1L, 3L) %in% c("Pr(", "p-v"), eps.Pvalue = .Machine$double.eps,
                                   na.print = "NA", quote = FALSE, right = TRUE, ...) {
  if (is.null(d <- dim(x)) || length(d) != 2L) {
    stop("'x' must be coefficient matrix/data frame")
  }
  nc <- d[2L]
  if (is.null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  } else if (P.values && !has.Pvalue) {
    stop("'P.values' is TRUE, but 'has.Pvalue' is not")
  }
  if (has.Pvalue && !P.values) {
    d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
    nc <- nc - 1
    has.Pvalue <- FALSE
  } else {
    xm <- data.matrix(x)
  }
  k <- nc - has.Pvalue - (if (missing(tst.ind)) {
    1
  } else {
    length(tst.ind)
  })
  if (!missing(cs.ind) && length(cs.ind) > k) {
    stop("wrong k / cs.ind")
  }
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  ok <- !(ina <- is.na(xm))
  for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
  if (length(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs !=
        0])) {
        floor(log10(range(acs[acs != 0], finite = TRUE)))
      } else {
        0
      }
      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits -
        digmin)), digits = digits)
    }
  }
  if (length(tst.ind)) {
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
      digits = digits
    )
  }
  if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc)))) {
    for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
  }
  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue) {
    ok[, -nc]
  } else {
    ok
  }
  x1 <- Cf[okP]
  dec <- getOption("OutDec")
  if (dec != ".") {
    x1 <- chartr(dec, ".", x1)
  }
  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
  if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(
      1L,
      digits - 1L
    ))
  }
  if (any(ina)) {
    Cf[ina] <- na.print
  }
  if (any(inan <- is.nan(xm))) {
    Cf[inan] <- "NaN"
  }
  if (P.values) {
    if (!is.logical(signif.stars) || is.na(signif.stars)) {
      warning("option \"show.signif.stars\" is invalid: assuming TRUE")
      signif.stars <- TRUE
    }
    if (any(okP <- ok[, nc])) {
      pv <- as.vector(xm[, nc])
      Cf[okP, nc] <- format.pval(pv[okP],
        digits = dig.tst,
        eps = eps.Pvalue
      )
      signif.stars <- signif.stars && any(pv[okP] < 0.1)
      if (signif.stars) {
        # MON
        max_pv <- max(pv, na.rm = TRUE)
        if (max_pv < 1) {
          max_cut <- 1
        } else {
          max_cut <- ceiling(max_pv * 2) / 2
        }

        cutpoints_vec <- c(0, 0.001, 0.01, 0.05, 0.1, max_cut)

        Signif <- stats::symnum(pv,
          corr = FALSE, na = FALSE,
          cutpoints = cutpoints_vec,
          symbols = c("***", "**", "*", ".", " ")
        )
        Cf <- cbind(Cf, format(Signif))
      }
    } else {
      signif.stars <- FALSE
    }
  } else {
    signif.stars <- FALSE
  }
  # print.default(Cf,
  #               quote = quote, right = right, na.print = na.print,
  #               ...
  # )
  if (signif.stars && signif.legend) {
    if ((w <- getOption("width")) < nchar(sleg <- attr(
      Signif,
      "legend"
    ))) {
      sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
    }
    cat("---\nSignif. codes:  ", sleg, sep = "", fill = w +
      4 + max(nchar(sleg, "bytes") - nchar(sleg)))
  }
  invisible(x)
}
