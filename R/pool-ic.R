#' Pooling of multiple-imputation PAMM fits
#'
#' Inference for interval-censored PAMMs (\code{\link{pamm_ic}}) pools the
#' \code{m} re-fits by drawing from each fit's empirical-Bayes posterior
#' \eqn{N(\hat\beta^{(m)}, V_\beta^{(m)})} and propagating every draw through the
#' quantity of interest using that fit's \emph{own} design matrix, then taking
#' empirical quantiles of the combined draws. Because mgcv's identifiability
#' constraints make the (centered) spline basis depend on each imputed data set,
#' the design matrix is \emph{not} shared across fits, so each fit must be
#' evaluated with its own \code{lpmatrix}. The combined draws are a sample from
#' the multiple-imputation mixture posterior
#' \eqn{(1/M)\sum_m N(\hat\beta^{(m)}, V_\beta^{(m)})}, whose variance reproduces
#' the within- plus between-imputation (Rubin) variance. Point estimates are the
#' average of the per-fit point estimates (the MI estimate).
#'
#' These methods are dispatched automatically by \code{\link{add_hazard}},
#' \code{\link{add_cumu_hazard}}, \code{\link{add_surv_prob}} and
#' \code{\link{add_cif}} when given a \code{pamm_ic} object.
#'
#' @name pamm_ic_pooling
#' @keywords internal
NULL

# Average a per-fit point estimate (the value column produced by an `add_*`
# default with ci = FALSE) across the imputation fits -> MI point estimate.
pooled_point <- function(object, newdata, adder, value_col, ...) {
  preds <- lapply(object[["fits"]], function(f) {
    adder(newdata, f, ci = FALSE, ...)[[value_col]]
  })
  rowMeans(do.call(cbind, preds))
}

# Per-group cumulative sum of intlen * hazard, applied column-wise to a draw
# matrix (one column per posterior draw). Rows are assumed time-ordered within
# each group, as produced by make_newdata().
ic_group_cumsum <- function(mat, intlen, grp) {
  out <- mat
  for (g in unique(grp)) {
    ix <- which(grp == g)
    out[ix, ] <- apply(mat[ix, , drop = FALSE] * intlen[ix], 2, cumsum)
  }
  out
}

# Pooled posterior draws of hazard / cumulative hazard / survival, evaluated with
# each fit's own design matrix, returning lower/upper quantiles.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef quantile
ic_ci_draws <- function(object, newdata, nsim, kind, alpha, time_var,
  interval_length) {

  fits <- object[["fits"]]
  m    <- length(fits)
  per  <- ceiling(nsim / m)
  tv   <- resolve_time_var(time_var, fits[[1]], newdata)

  nd <- newdata
  intlen <- NULL
  grp    <- NULL
  if (kind != "hazard") {
    nd     <- reconstruct_intlen(nd, time_var = tv, interval_length = interval_length)
    intlen <- nd[[interval_length]]
    grp    <- group_indices(nd)
  }

  cols <- lapply(fits, function(f) {
    X <- predict.gam(f, newdata = nd, type = "lpmatrix")
    B <- rmvnorm(per, mean = coef(f), sigma = f[["Vp"]])
    H <- exp(X %*% t(B))                       # nrow x per hazard draws
    if (kind == "hazard") {
      return(H)
    }
    C <- ic_group_cumsum(H, intlen, grp)
    if (kind == "cumu") C else exp(-C)
  })
  M <- do.call(cbind, cols)

  list(
    lower = apply(M, 1, quantile, probs = alpha / 2,     na.rm = TRUE),
    upper = apply(M, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE))

}

#' @rdname add_hazard
#' @param alpha Significance level for pooled confidence intervals (a
#'   \eqn{(1-\alpha)} interval).
#' @param nsim Total number of pooled posterior draws used for the interval.
#' @export
add_hazard.pamm_ic <- function(
  newdata, object, ci = TRUE, alpha = 0.05, nsim = 500L, time_var = NULL, ...) {

  newdata[["hazard"]] <- pooled_point(object, newdata, add_hazard, "hazard",
    time_var = time_var, ...)
  if (ci) {
    d <- ic_ci_draws(object, newdata, nsim, "hazard", alpha, time_var, "intlen")
    newdata[["ci_lower"]] <- d[["lower"]]
    newdata[["ci_upper"]] <- d[["upper"]]
  }
  newdata

}

#' @rdname add_hazard
#' @export
add_cumu_hazard.pamm_ic <- function(
  newdata, object, ci = TRUE, alpha = 0.05, nsim = 500L, time_var = NULL,
  interval_length = "intlen", ...) {

  newdata[["cumu_hazard"]] <- pooled_point(object, newdata, add_cumu_hazard,
    "cumu_hazard", time_var = time_var, interval_length = interval_length, ...)
  if (ci) {
    d <- ic_ci_draws(object, newdata, nsim, "cumu", alpha, time_var,
      interval_length)
    newdata[["cumu_lower"]] <- d[["lower"]]
    newdata[["cumu_upper"]] <- d[["upper"]]
  }
  newdata

}

#' @rdname add_surv_prob
#' @param alpha Significance level for pooled confidence intervals.
#' @param nsim Total number of pooled posterior draws used for the interval.
#' @export
add_surv_prob.pamm_ic <- function(
  newdata, object, ci = TRUE, alpha = 0.05, nsim = 500L, time_var = NULL,
  interval_length = "intlen", ...) {

  newdata[["surv_prob"]] <- pooled_point(object, newdata, add_surv_prob,
    "surv_prob", time_var = time_var, interval_length = interval_length, ...)
  if (ci) {
    d <- ic_ci_draws(object, newdata, nsim, "surv", alpha, time_var,
      interval_length)
    newdata[["surv_lower"]] <- d[["lower"]]
    newdata[["surv_upper"]] <- d[["upper"]]
  }
  newdata

}

# Pooled CIF draws for a single cause-by-covariate group, evaluated with each
# fit's own design matrices. Returns the [nrow x nsim] matrix of CIF draws.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef predict
ic_cif_draws_group <- function(group_df, object, per, cause_var,
  interval_length) {

  fits         <- object[["fits"]]
  cause_levels <- as.factor(levels(group_df[[cause_var]]))
  cause_data   <- unique(group_df[[cause_var]])
  dt           <- group_df[[interval_length]]

  cols <- lapply(fits, function(f) {
    B <- rmvnorm(per, mean = coef(f), sigma = f[["Vp"]])
    hazards <- lapply(cause_levels, function(cl) {
      dfc <- group_df
      dfc[[cause_var]] <- factor(cl, levels = levels(cause_levels))
      X <- predict(f, dfc, type = "lpmatrix")
      exp(X %*% t(B))                          # nrow x per
    })
    names(hazards) <- as.character(cause_levels)
    total_hazard <- Reduce(`+`, hazards)
    overall_surv <- apply(total_hazard, 2, function(z) exp(-cumsum(z * dt)))
    survival     <- rbind(1, overall_surv[-nrow(overall_surv), , drop = FALSE])
    hazard       <- hazards[[as.character(cause_data)]]
    cif_inc <- (hazard / total_hazard) * survival * (1 - exp(-total_hazard * dt))
    apply(cif_inc, 2, cumsum)
  })

  do.call(cbind, cols)

}

#' @rdname add_cif
#' @param alpha Significance level for pooled confidence intervals.
#' @param nsim Total number of pooled posterior draws used for the interval.
#' @export
add_cif.pamm_ic <- function(
  newdata, object, ci = TRUE, alpha = 0.05, nsim = 500L, cause_var = "cause",
  time_var = NULL, interval_length = "intlen", ...) {

  if (!identical(object[["type"]], "cr")) {
    stop("add_cif() requires a competing-risks `pamm_ic` (see pamm_ic_cr()).",
      call. = FALSE)
  }
  fit1     <- object[["fits"]][[1]]
  m        <- length(object[["fits"]])
  per      <- ceiling(nsim / m)
  time_var <- resolve_time_var(time_var, fit1, newdata)
  joindata <- reconstruct_cutpoints(newdata, fit1, time_var, interval_length)

  joindata <- map_dfr(
    split(joindata, group_indices(joindata)),
    function(.x) {
      .x  <- arrange(.x, .data[[time_var]])
      cif <- ic_cif_draws_group(.x, object, per, cause_var, interval_length)
      .x[["cif"]] <- pmin(pmax(rowMeans(cif), 0), 1)
      if (ci) {
        .x[["cif_lower"]] <- pmin(pmax(
          apply(cif, 1, quantile, alpha / 2, na.rm = TRUE), 0), 1)
        .x[["cif_upper"]] <- pmin(pmax(
          apply(cif, 1, quantile, 1 - alpha / 2, na.rm = TRUE), 0), 1)
      }
      .x
    })

  suppressMessages(newdata %>% left_join(joindata))

}
