# Shared functions for the mixed model for repeated measures (MMRM) and its multiple-imputation
# sensitivity analyses. Sourced by mmrm_fit.R and mmrm_mi.R. Base R + mmrm + jsonlite only.
suppressMessages({library(mmrm); library(jsonlite)})

# Covariance structures, tried in order: unstructured, heterogeneous Toeplitz, heterogeneous first-order
# ante-dependence, heterogeneous compound symmetry, compound symmetry. The first structure that converges
# with a positive definite estimate is used.
COV_STRUCTURES <- c("us", "toeph", "adh", "csh", "cs")
# A converged estimate whose smallest eigenvalue is a negligible fraction of its largest is numerically
# singular (its conditional distributions are degenerate, so it cannot be refitted on imputed datasets)
# and is treated as not positive definite.
PD_TOL <- 1e-6

# Long data: one row per (patient, visit) with columns pid, visit ("week_4" ... "week_52"),
# uefm, uefm_base. Visit levels are ordered by week; the baseline covariate is centred at the
# analysis-set mean so that each visit coefficient is directly the LS mean at the mean baseline.
read_long <- function(path, base_mean = NULL) {
  d <- read.csv(path, colClasses = c(pid = "character", visit = "character"))
  weeks <- as.numeric(sub("week_", "", d$visit))
  d$visit <- factor(d$visit, levels = unique(d$visit)[order(unique(weeks))])
  d$pid <- factor(d$pid)
  if (is.null(base_mean)) base_mean <- mean(tapply(d$uefm_base, d$pid, function(x) x[1]))
  d$uefm_base_c <- d$uefm_base - base_mean
  attr(d, "base_mean") <- base_mean
  d
}

week_of <- function(visit) as.numeric(sub("week_", "", as.character(visit)))

# Fit Y_it = mu_t + beta * (Y_i0 - mean Y_0) + e_it, REML, requested covariance structure.
fit_mmrm <- function(d, cov, df_method = "Kenward-Roger") {
  form <- as.formula(sprintf("uefm ~ 0 + visit + uefm_base_c + %s(visit | pid)", cov))
  fit <- tryCatch(mmrm(form, data = d, reml = TRUE, control = mmrm_control(method = df_method)),
                  error = function(e) e)
  if (inherits(fit, "error")) return(list(ok = FALSE, reason = gsub("\\s+", " ", conditionMessage(fit))))
  ev <- eigen(VarCorr(fit), symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= PD_TOL * max(ev))
    return(list(ok = FALSE, reason = sprintf("covariance estimate not positive definite (smallest/largest eigenvalue %.1e)", min(ev) / max(ev))))
  list(ok = TRUE, fit = fit)
}

# Try the structures in order; return the first acceptable fit and a one-line record of each attempt.
fit_first_acceptable <- function(d, structures = COV_STRUCTURES, df_method = "Kenward-Roger") {
  attempts <- character(0)
  for (cov in structures) {
    r <- fit_mmrm(d, cov, df_method)
    attempts <- c(attempts, sprintf("%s: %s", cov, if (r$ok) "accepted" else r$reason))
    if (r$ok) return(list(fit = r$fit, cov = cov, attempts = attempts))
  }
  stop("No covariance structure converged with a positive definite estimate: ", paste(attempts, collapse = "; "))
}

# Estimate k' beta with its df (Kenward-Roger, per the fit) and 95% CI.
# `null` shifts the tested value: p tests k' beta = null.
contrast <- function(fit, k, null = 0) {
  x <- df_1d(fit, k)
  t <- (x$est - null) / x$se
  tc <- qt(0.975, x$df)
  data.frame(estimate = x$est - null, se = x$se, df = x$df,
             lower = x$est - null - tc * x$se, upper = x$est - null + tc * x$se,
             p = 2 * pt(-abs(t), x$df))
}

# Zero coefficient vector with the fit's names; set entries by name.
k_vector <- function(fit, ...) {
  k <- setNames(rep(0, length(coef(fit))), names(coef(fit)))
  set <- list(...)
  for (nm in names(set)) k[nm] <- set[[nm]]
  k
}

# Change from baseline at every visit: the visit's LS mean at the mean baseline minus the mean baseline.
visit_changes <- function(fit, base_mean, visits) {
  do.call(rbind, lapply(visits, function(v) {
    k <- k_vector(fit); k[paste0("visit", v)] <- 1
    cbind(contrast = v, week = week_of(v), contrast(fit, k, null = base_mean))
  }))
}

# Difference between the week-52 and week-12 LS means.
week52_minus_week12 <- function(fit) {
  cbind(contrast = "week_52 - week_12", week = NA,
        contrast(fit, k_vector(fit, visitweek_12 = -1, visitweek_52 = 1)))
}
