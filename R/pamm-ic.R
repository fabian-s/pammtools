#' Fit a PAMM to interval-censored data via multiple imputation
#'
#' Fits a piecewise exponential additive (mixed) model to interval-censored
#' time-to-event data using a multiple-imputation (MI) and re-fit strategy: exact
#' event times are repeatedly drawn from the model-based conditional distribution
#' \eqn{p(T \mid L < T \le R, x, \theta)} (see \code{\link{impute_ic_times}}),
#' each completed data set is transformed to PED format with the standard
#' (right-censored) pipeline and re-fit, and the resulting fits are pooled for
#' inference with the existing \code{add_*} family (see \code{\link{add_surv_prob}}
#' and the \code{pamm_ic} methods).
#'
#' An imputed event time is an exact event time, so once imputation has produced
#' it, the entire downstream pipeline (\code{\link{split_data}} -> \code{\link{pamm}}
#' -> \code{add_*}) is reused unchanged. The interval cut-points are resolved once
#' and shared across all imputations, which keeps the spline bases - and hence the
#' design matrices - identical across fits, a precondition for valid pooling.
#'
#' @param formula A two-sided formula whose left-hand side is an interval-censored
#'   response \code{Surv(L, R, type = "interval2")} and whose right-hand side lists
#'   the covariates to retain (as in \code{\link{as_ped}}).
#' @param data A data frame in standard (one row per subject) format.
#' @param model_formula Optional model formula passed to \code{\link{pamm}}
#'   (e.g.\ \code{ped_status ~ s(tend) + x}). If \code{NULL}, a default
#'   \code{ped_status ~ s(tend) + <covariates>} formula is constructed.
#' @param cut Optional fixed vector of interval cut-points shared across all
#'   imputations. If \code{NULL}, the finite interval endpoints are used.
#' @param max_time Optional cap on the cut-points.
#' @param m Number of imputations (default 10).
#' @param proper Logical; if \code{TRUE} (default, "proper" MI) a coefficient
#'   vector is drawn from the posterior \eqn{N(\hat\beta, V_\beta)} of the
#'   initialiser fit before each imputation, propagating parameter uncertainty.
#'   \code{FALSE} uses the point estimate and is intended for diagnostics only.
#' @param init Initialiser for the first fit: \code{"midpoint"} (default) or
#'   \code{"uniform"} imputation within each interval.
#' @param id Name of the subject identifier column.
#' @param engine Estimation engine passed to \code{\link{pamm}} (\code{"gam"} or
#'   \code{"bam"}).
#' @param ... Further arguments passed to \code{\link{pamm}} / \code{mgcv}.
#' @return An object of class \code{pamm_ic}: a list with the \code{m} pooled
#'   \code{pamm} fits (\code{fits}), the initialiser fit (\code{init_fit}), the
#'   parsed interval bounds (\code{ic}), the shared \code{cut}, and metadata. Each
#'   element of \code{fits} is an ordinary \code{pamm} object, so any single
#'   \code{add_*} call works on it directly.
#' @seealso \code{\link{impute_ic_times}}, \code{\link{add_surv_prob}}
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef as.formula
#' @export
pamm_ic <- function(
  formula,
  data,
  model_formula = NULL,
  cut = NULL,
  max_time = NULL,
  m = 10L,
  proper = TRUE,
  init = c("midpoint", "uniform"),
  id = "id",
  engine = "gam",
  ...
) {
  init <- match.arg(init)
  assert_count(m, positive = TRUE)

  ped0 <- as_ped_ic(data, formula, cut = cut, max_time = max_time, id = id)
  ic <- attr(ped0, "ic")
  cut <- attr(ped0, "breaks")

  if (init == "uniform") {
    # replace midpoint initialiser by a single uniform draw
    t_unif <- draw_uniform_ic(ic, cut)
    ped0 <- build_ic_ped(ic, t_unif, cut, formula, id)
  }

  if (is.null(model_formula)) {
    model_formula <- default_pamm_formula(formula, data = data, id = id)
  }
  fit0 <- pamm(model_formula, data = ped0, engine = engine, ...)
  cache <- ic_pred_cache(fit0, ic, cut)

  fits <- vector("list", m)
  for (mm in seq_len(m)) {
    beta_mm <- if (proper) {
      as.numeric(rmvnorm(1, mean = coef(fit0), sigma = fit0[["Vp"]]))
    } else {
      coef(fit0)
    }
    t_imp <- impute_ic_times(fit0, ic, cut, beta = beta_mm, cache = cache)
    ped_m <- build_ic_ped(ic, t_imp, cut, formula, id)
    fits[[mm]] <- pamm(model_formula, data = ped_m, engine = engine, ...)
  }

  structure(
    list(
      fits = fits,
      init_fit = fit0,
      ic = ic,
      cut = cut,
      formula = formula,
      model_formula = model_formula,
      m = m,
      proper = proper,
      id_var = id,
      type = "single"
    ),
    class = c("pamm_ic", "list")
  )
}

#' Fit a competing-risks PAMM to interval-censored data via multiple imputation
#'
#' Competing-risks extension of \code{\link{pamm_ic}}. The event time is drawn
#' from the all-cause conditional hazard within \eqn{(L, R]} and a cause is
#' assigned: observed causes are retained (with the time drawn so that it follows
#' the cause-specific conditional density, via rejection), unknown causes are
#' sampled with probability proportional to the cause-specific hazards at the
#' imputed time (see \code{\link{impute_ic_cr}}). Each completed data set is
#' transformed with \code{\link{as_ped_cr}} (cause-specific hazards) and re-fit.
#' Cf. Delord & Genin (2016) for MI of interval-censored competing-risks data.
#'
#' @inheritParams pamm_ic
#' @param cause Name of the column in \code{data} giving the observed cause for
#'   events (any factor/character coding). Rows with the censoring code are
#'   treated as right-censored; \code{NA} marks an event with unknown cause.
#' @param censor_code Value of \code{cause} that encodes censoring (default 0).
#' @return An object of class \code{pamm_ic} with \code{type = "cr"}; \code{fits}
#'   are cause-specific (stacked \code{ped_cr}) \code{pamm} objects and
#'   \code{cause_levels} records the competing causes.
#' @seealso \code{\link{pamm_ic}}, \code{\link{add_cif}}
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef as.formula
#' @export
pamm_ic_cr <- function(
  formula,
  data,
  cause,
  model_formula = NULL,
  cut = NULL,
  max_time = NULL,
  m = 10L,
  proper = TRUE,
  censor_code = 0L,
  id = "id",
  engine = "gam",
  ...
) {
  assert_count(m, positive = TRUE)
  assert_string(cause)
  assert_subset(cause, names(data))

  ic <- parse_ic_surv(formula, data, id = id)
  cut <- resolve_ic_cut(ic, cut = cut, max_time = max_time)

  cause_raw <- data[[cause]]
  is_cens <- !is.na(cause_raw) &
    as.character(cause_raw) == as.character(censor_code)
  # for right-censored survival rows, force censoring regardless of cause column
  is_cens <- is_cens | as.character(ic[["ic_kind"]]) == "right"
  cause_known <- ifelse(is_cens, NA, as.character(cause_raw))
  cause_levels <- sort(unique(stats::na.omit(cause_known)))
  if (length(cause_levels) < 2) {
    stop(
      "Fewer than two competing causes found in `",
      cause,
      "`.",
      call. = FALSE
    )
  }

  if (is.null(model_formula)) {
    model_formula <- default_pamm_formula(
      formula,
      data = data,
      id = id,
      by_cause = TRUE,
      exclude = cause
    )
  }

  # initialiser: midpoint times + observed causes (unknown causes drawn from the
  # marginal cause distribution) -> stacked cause-specific PED
  L <- ic[["ic_L"]]
  R <- ic[["ic_R"]]
  t_mid <- pmin(ifelse(ic[["ic_kind"]] == "left", R / 2, (L + R) / 2), max(cut))
  cause0 <- cause_known
  unknown <- which(!is_cens & is.na(cause0))
  if (length(unknown)) {
    cause0[unknown] <- sample(cause_levels, length(unknown), replace = TRUE)
  }
  ped0 <- build_ic_ped_cr(
    ic,
    t_mid,
    cause0,
    is_cens,
    cut,
    formula,
    id,
    cause_levels,
    censor_code,
    cause_var = cause
  )
  fit0 <- pamm(model_formula, data = ped0, engine = engine, ...)
  cache <- ic_pred_cache(fit0, ic, cut, cause_levels = cause_levels)

  fits <- vector("list", m)
  for (mm in seq_len(m)) {
    beta_mm <- if (proper) {
      as.numeric(rmvnorm(1, mean = coef(fit0), sigma = fit0[["Vp"]]))
    } else {
      coef(fit0)
    }
    imp <- impute_ic_cr(
      fit0,
      ic,
      cut,
      beta = beta_mm,
      cache = cache,
      cause_known = cause_known
    )
    ped_m <- build_ic_ped_cr(
      ic,
      imp[["time"]],
      imp[["cause"]],
      is_cens,
      cut,
      formula,
      id,
      cause_levels,
      censor_code,
      cause_var = cause
    )
    fits[[mm]] <- pamm(model_formula, data = ped_m, engine = engine, ...)
  }

  structure(
    list(
      fits = fits,
      init_fit = fit0,
      ic = ic,
      cut = cut,
      formula = formula,
      model_formula = model_formula,
      m = m,
      proper = proper,
      id_var = id,
      cause_levels = cause_levels,
      type = "cr"
    ),
    class = c("pamm_ic", "list")
  )
}

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# Construct a default model formula `ped_status ~ s(tend) [+ by-cause] + covars`.
default_pamm_formula <- function(
  formula,
  data = NULL,
  id = "id",
  by_cause = FALSE,
  exclude = character()
) {
  rhs <- resolve_rhs_vars(formula, data = data, exclude = c(id, exclude))
  base <- if (by_cause) "s(tend, by = cause) + cause" else "s(tend)"
  terms <- c(base, if (length(rhs)) paste(rhs, collapse = " + "))
  stats::as.formula(paste("ped_status ~", paste(terms, collapse = " + ")))
}

# Uniform initial draw within (L, R] (R/2 for left-censored).
draw_uniform_ic <- function(ic, cut) {
  L <- ic[["ic_L"]]
  R <- pmin(ic[["ic_R"]], max(cut))
  u <- stats::runif(nrow(ic))
  ifelse(
    as.character(ic[["ic_kind"]]) %in% c("exact", "right"),
    L,
    L + u * (R - L)
  )
}

# Build a single-event PED from imputed exact times via the standard pipeline.
build_ic_ped <- function(ic, t_imp, cut, formula, id) {
  rhs_vars <- resolve_rhs_vars(formula, data = ic, exclude = id)
  evd <- drop_zero_followup(ic_event_data(ic, t_imp), warn = FALSE)
  ped_form <- stats::as.formula(
    paste0(
      "Surv(.ped_time, .ped_status) ~ ",
      paste0(unique(c(rhs_vars, id)), collapse = " + ")
    )
  )
  split_data(ped_form, data = evd, cut = cut, id = id)
}

# Build a competing-risks PED (stacked ped_cr) from imputed times and causes.
build_ic_ped_cr <- function(
  ic,
  time,
  cause,
  is_cens,
  cut,
  formula,
  id,
  cause_levels,
  censor_code,
  cause_var = NULL
) {
  rhs_vars <- resolve_rhs_vars(formula, data = ic, exclude = c(id, cause_var))
  dat <- ic
  dat[["ic_L"]] <- dat[["ic_R"]] <- dat[["ic_kind"]] <- NULL
  dat[[".ped_time"]] <- ifelse(
    as.character(ic[["ic_kind"]]) %in%
      c("exact", "right"),
    ic[["ic_L"]],
    time
  )
  status_cr <- ifelse(is_cens, censor_code, cause)
  dat[[".status_cr"]] <- factor(
    status_cr,
    levels = c(censor_code, cause_levels)
  )

  keep <- dat[[".ped_time"]] > 0
  dat <- dat[keep, , drop = FALSE]

  cr_form <- stats::as.formula(
    paste0(
      "Surv(.ped_time, .status_cr) ~ ",
      paste0(unique(c(rhs_vars, id)), collapse = " + ")
    )
  )
  as_ped(dat, formula = cr_form, cut = cut, censor_code = censor_code, id = id)
}

#' @rdname pamm_ic
#' @param x,object A \code{pamm_ic} object.
#' @param ... Ignored.
#' @export
print.pamm_ic <- function(x, ...) {
  cat("PAMM fit to interval-censored data via multiple imputation\n")
  cat("  type        :", x[["type"]], "\n")
  cat(
    "  imputations :",
    x[["m"]],
    if (x[["proper"]]) "(proper)" else "(improper)",
    "\n"
  )
  cat(
    "  cut-points  :",
    length(x[["cut"]]),
    "breaks in [",
    format(min(x[["cut"]])),
    ",",
    format(max(x[["cut"]])),
    "]\n"
  )
  n_ic <- sum(as.character(x[["ic"]][["ic_kind"]]) %in% c("interval", "left"))
  cat(
    "  subjects    :",
    nrow(x[["ic"]]),
    "(",
    n_ic,
    "interval/left-censored )\n"
  )
  invisible(x)
}

#' @rdname pamm_ic
#' @export
summary.pamm_ic <- function(object, ...) {
  cat("Initialiser fit:\n")
  print(summary(object[["init_fit"]], ...))
  cat(
    "\nPooled over",
    object[["m"]],
    "imputations. Use add_*() for pooled",
    "quantities of interest.\n"
  )
  invisible(object)
}
