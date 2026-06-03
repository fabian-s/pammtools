#' Slim down a fitted PAMM for storage inside a \code{pamm_ic} object
#'
#' Removes the per-observation slots (model frame, fitted values, residuals,
#' working weights, ...) and the \code{call} (which captures the full PED data),
#' none of which are needed for the downstream multiple-imputation pooling: the
#' pooled \code{add_*} methods only require each fit's \code{coefficients},
#' \code{Vp}/\code{Ve} and the smooth/parametric structure used by
#' \code{predict(type = "lpmatrix")}. Stripping makes the stored size independent
#' of the data set size, so memory does not blow up with many imputations.
#'
#' @param fit A fitted \code{pamm}/\code{gam} object.
#' @return The same object with large per-observation slots removed; class and
#'   everything needed for \code{predict}/\code{coef}/\code{vcov} are retained.
#' @keywords internal
strip_pamm_fit <- function(fit) {
  drop <- c(
    "model", "y", "residuals", "fitted.values", "linear.predictors",
    "weights", "prior.weights", "offset", "wt", "working.weights", "hat",
    "z", "w", "std.rsd", "na.action")
  for (s in intersect(drop, names(fit))) fit[[s]] <- NULL
  fit[["call"]] <- NULL
  fit
}

# Fit one imputation, capture its summary (for median-p pooling) while the full
# object is available, then return a stripped fit + PED row count.
fit_strip_summarise <- function(model_formula, data, engine, ...) {
  f <- pamm(model_formula, data = data, engine = engine, ...)
  list(summary = summary(f), fit = strip_pamm_fit(f), n = nrow(data))
}

# Median over imputations of a per-fit table column (the "median-p rule" for
# pooling significance across MI, shown to work well for GAMs by Bolt et al.
# 2022, BMC Med Res Methodol, doi:10.1186/s12874-022-01613-w). Other columns
# (edf, statistics) are averaged.
pool_param_table <- function(ptabs, est, V, nsdf) {
  if (nsdf < 1L || is.null(ptabs[[1]]) || nrow(ptabs[[1]]) == 0L) return(NULL)
  rn   <- rownames(ptabs[[1]])
  pcol <- ncol(ptabs[[1]])
  pmat <- vapply(ptabs, function(t) t[, pcol], numeric(nrow(ptabs[[1]])))
  pmed <- if (is.matrix(pmat)) apply(pmat, 1, stats::median) else stats::median(pmat)
  ind  <- seq_len(nsdf)
  Est  <- est[ind]
  SE   <- sqrt(diag(V))[ind]
  out  <- cbind("Estimate" = Est, "Std. Error" = SE,
    "z value" = Est / SE, "Pr(>|z|)" = pmed)
  rownames(out) <- rn
  out
}

pool_smooth_table <- function(stabs) {
  if (is.null(stabs[[1]]) || nrow(stabs[[1]]) == 0L) return(NULL)
  rn  <- rownames(stabs[[1]])
  arr <- simplify2array(stabs)                 # rows x cols x m
  if (length(dim(arr)) == 2L) arr <- array(arr, c(dim(arr), 1L))
  pc  <- dim(arr)[2]
  out <- cbind(
    "edf"      = apply(arr[, 1, , drop = FALSE], 1, mean),
    "Ref.df"   = apply(arr[, 2, , drop = FALSE], 1, mean),
    "statistic" = apply(arr[, 3, , drop = FALSE], 1, mean),
    "p-value"  = apply(arr[, pc, , drop = FALSE], 1, stats::median))
  rownames(out) <- rn
  out
}

#' Pool a list of (stripped) imputation fits into a single pooled fit
#'
#' Combines the \code{m} imputation fits with Rubin's rules: the pooled
#' coefficients are the average \eqn{\bar Q}, and the pooled covariances inflate
#' the within-imputation covariance by the between-imputation variance,
#' \eqn{V = \bar W + (1 + 1/m) B} (applied to both \code{Vp} and \code{Ve}), so
#' the reported standard errors include the additional MI variability. Term
#' p-values are pooled with the median-p rule (see references in
#' \code{\link{strip_pamm_fit}}). Returns a pooled \code{gam}-like object plus the
#' pooled coefficient/smooth tables and per-coefficient MI diagnostics
#' (relative increase in variance, fraction of missing information).
#'
#' @param fits List of stripped imputation fits.
#' @param smry List of \code{summary.gam} objects, one per fit (computed before
#'   stripping).
#' @keywords internal
#' @importFrom stats coef cov median
pool_pamm_fits <- function(fits, smry) {
  m  <- length(fits)
  p  <- length(coef(fits[[1]]))
  cf <- vapply(fits, coef, numeric(p))                       # p x m
  Qbar <- rowMeans(cf)

  Wp <- Reduce(`+`, lapply(fits, function(f) f[["Vp"]])) / m
  We <- Reduce(`+`, lapply(fits, function(f)
    if (!is.null(f[["Ve"]])) f[["Ve"]] else f[["Vp"]])) / m
  B  <- if (m > 1) stats::cov(t(cf)) else matrix(0, p, p)
  Vp <- Wp + (1 + 1 / m) * B
  Ve <- We + (1 + 1 / m) * B

  edf_mat <- vapply(fits, function(f) as.numeric(f[["edf"]]), numeric(p))
  edf <- rowMeans(edf_mat)

  # pooled gam-like skeleton (structure from fit 1; coefficients/cov replaced)
  g <- fits[[1]]
  g[["coefficients"]] <- Qbar
  g[["Vp"]] <- Vp
  g[["Ve"]] <- Ve
  g[["edf"]] <- edf

  # per-coefficient MI diagnostics
  wdiag <- diag(Wp); bdiag <- diag(B)
  riv <- ifelse(wdiag > 0, (1 + 1 / m) * bdiag / wdiag, NA_real_)
  dfbr <- (m - 1) * (1 + 1 / riv)^2
  fmi <- (riv + 2 / (dfbr + 3)) / (riv + 1)
  names(riv) <- names(fmi) <- names(Qbar)

  list(
    gam      = g,
    coefficients = Qbar,
    Vp       = Vp,
    Ve       = Ve,
    edf      = edf,
    p.table  = pool_param_table(lapply(smry, `[[`, "p.table"), Qbar, Vp,
      g[["nsdf"]]),
    s.table  = pool_smooth_table(lapply(smry, `[[`, "s.table")),
    riv      = riv,
    fmi      = fmi,
    family   = g[["family"]])
}

#' @rdname pamm_ic
#' @param x,object A \code{pamm_ic} object.
#' @param ... Passed on (ignored for \code{print}).
#' @export
print.pamm_ic <- function(x, ...) {
  cat("Piecewise exponential additive model for interval-censored data\n")
  cat("  via multiple imputation (", x[["m"]],
    if (x[["proper"]]) " proper" else " improper", " imputations)\n", sep = "")
  cat("  task        : ",
    if (identical(x[["type"]], "cr")) "competing risks" else "single event", "\n",
    sep = "")
  cat("  model       : ", deparse(x[["model_formula"]]), "\n", sep = "")
  n_ic <- sum(as.character(x[["ic"]][["ic_kind"]]) %in% c("interval", "left"))
  cat("  subjects    : ", nrow(x[["ic"]]), " (", n_ic,
    " interval/left-censored)\n", sep = "")
  cat("  cut-points  : ", length(x[["cut"]]), " breaks in [",
    format(min(x[["cut"]])), ", ", format(max(x[["cut"]])), "]\n\n", sep = "")
  pt <- x[["pooled"]][["p.table"]]
  if (!is.null(pt)) {
    cat("Pooled parametric coefficients:\n")
    stats::printCoefmat(pt, signif.stars = FALSE, digits = 3L)
  }
  cat("\nUse summary() for the full pooled fit and add_*() for pooled",
    "quantities of interest.\n")
  invisible(x)
}

#' @rdname pamm_ic
#' @export
summary.pamm_ic <- function(object, ...) {
  p <- object[["pooled"]]
  structure(
    list(
      type = object[["type"]],
      m = object[["m"]],
      proper = object[["proper"]],
      model_formula = object[["model_formula"]],
      family = p[["family"]][["family"]],
      n_obs = object[["n_obs"]],
      n_subj = nrow(object[["ic"]]),
      n_ic = sum(as.character(object[["ic"]][["ic_kind"]]) %in%
        c("interval", "left")),
      cut = object[["cut"]],
      p.table = p[["p.table"]],
      s.table = p[["s.table"]],
      fmi = p[["fmi"]]),
    class = "summary.pamm_ic")
}

#' @rdname pamm_ic
#' @export
print.summary.pamm_ic <- function(x, ...) {
  cat("Pooled PAMM summary (multiple imputation for interval-censored data)\n\n")
  cat("Task           :",
    if (identical(x[["type"]], "cr")) "competing risks" else "single event", "\n")
  cat("Family         :", x[["family"]], "\n")
  cat("Model          :", deparse(x[["model_formula"]]), "\n")
  cat("Imputations    :", x[["m"]],
    if (x[["proper"]]) "(proper)" else "(improper)", "\n")
  cat("Subjects       :", x[["n_subj"]], "(", x[["n_ic"]],
    "interval/left-censored ); PED rows per fit:", x[["n_obs"]], "\n\n")

  if (!is.null(x[["p.table"]])) {
    cat("Parametric coefficients (Rubin-pooled estimates & SEs, median-p):\n")
    stats::printCoefmat(x[["p.table"]], has.Pvalue = TRUE, signif.stars = TRUE,
      digits = 3L)
  }
  if (!is.null(x[["s.table"]])) {
    cat("\nApproximate significance of smooth terms",
      "(mean edf, median-p over imputations):\n")
    stats::printCoefmat(x[["s.table"]], has.Pvalue = TRUE, signif.stars = TRUE,
      digits = 3L)
  }
  fmi <- x[["fmi"]][is.finite(x[["fmi"]])]
  if (length(fmi)) {
    cat("\nFraction of missing information (from interval censoring):",
      sprintf("median %.2f, max %.2f", stats::median(fmi), max(fmi)), "\n")
  }
  cat("\nStandard errors include within- + between-imputation variance",
    "(Rubin's rules);\np-values are medians over imputations.\n")
  invisible(x)
}
