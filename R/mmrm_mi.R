#!/usr/bin/env Rscript
# Multiple-imputation sensitivity analyses for the MMRM: impute each assessment that was due but not
# recorded, refit the primary model on every completed dataset, and pool by Rubin's rules
# (Barnard-Rubin degrees of freedom).
#   imputation=locf  carry-forward: the patient's last observed value + N(0, sigma_visit), sigma_visit = the
#                    primary fit's standard deviation at that visit.
#   imputation=mar   draws from the primary model's conditional distribution given the patient's observed
#                    values (missing at random; fixed effects drawn from their sampling distribution,
#                    covariance at its REML estimate); a penalty delta is subtracted from every imputed value
#                    after week 12 (tipping-point analysis).
#   estimate=trajectory  pooled change from baseline at every visit (run at delta = 0)
#   estimate=contrast    pooled week-52 minus week-12 difference at each delta; the sweep stops one point past
#                        the first delta whose 95% CI includes 0
# Usage: Rscript R/mmrm_mi.R observed=<csv> cells=<csv> primary=<summary.json> imputation=<locf|mar>
#                            estimate=<trajectory|contrast> M=50 seed=<int> deltas=0,0.5,... out=<csv>
#   cells columns: pid, visit, week, post12, last_obs_uefm, uefm_base
# Rscript passes the script path with spaces encoded as "~+~"
script_dir <- dirname(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE))
source(file.path(script_dir, "mmrm_functions.R"))
raw <- commandArgs(trailingOnly = TRUE)
arg <- setNames(as.list(sub("^[^=]*=", "", raw)), sub("=.*$", "", raw))   # key=value; the value may contain "="
for (key in c("observed", "cells", "primary", "imputation", "estimate", "M", "seed", "deltas", "out"))
  if (is.null(arg[[key]])) stop("missing argument ", key, "=")
if (!arg$imputation %in% c("locf", "mar")) stop("imputation must be locf or mar")
if (!arg$estimate %in% c("trajectory", "contrast")) stop("estimate must be trajectory or contrast")
primary <- fromJSON(arg$primary)
obs <- read_long(arg$observed, base_mean = primary$base_mean)
cells <- read.csv(arg$cells, colClasses = c(pid = "character", visit = "character"))
M <- as.integer(arg$M); deltas <- as.numeric(strsplit(arg$deltas, ",")[[1]])
df_method <- "Kenward-Roger"
set.seed(as.integer(arg$seed), kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
visits <- levels(obs$visit)

# --- imputation models ---------------------------------------------------------------------------------
sigma <- unlist(primary$visit_sd)
impute_locf <- function(delta) cells$last_obs_uefm + rnorm(nrow(cells), 0, sigma[cells$visit]) - delta * cells$post12

fit0 <- fit_mmrm(obs, primary$covariance, df_method)$fit          # the primary model on the observed data
beta_hat <- coef(fit0); L_beta <- t(chol(vcov(fit0)))
Sigma <- VarCorr(fit0); dimnames(Sigma) <- list(visits, visits)
cell_rows <- split(seq_len(nrow(cells)), cells$pid)
obs_rows <- split(obs, as.character(obs$pid))
# Conditional-distribution matrices depend only on which visits are observed and which are imputed, so they
# are computed once per observed/imputed visit pattern: A maps observed deviations to the imputed means,
# L is the Cholesky factor of the conditional covariance.
pattern <- sapply(names(cell_rows), function(id)
  paste(c(as.character(obs_rows[[id]]$visit), "|", cells$visit[cell_rows[[id]]]), collapse = ","))
pre <- lapply(unique(pattern), function(p) {
  id <- names(pattern)[match(p, pattern)]
  ov <- as.character(obs_rows[[id]]$visit); mv <- cells$visit[cell_rows[[id]]]
  A <- Sigma[mv, ov, drop = FALSE] %*% solve(Sigma[ov, ov, drop = FALSE])
  list(ov = ov, mv = mv, A = A, L = t(chol(Sigma[mv, mv, drop = FALSE] - A %*% Sigma[ov, mv, drop = FALSE])))
})
names(pre) <- unique(pattern)
impute_mar <- function(delta) {
  beta <- setNames(as.vector(beta_hat + L_beta %*% rnorm(length(beta_hat))), names(beta_hat))
  y <- numeric(nrow(cells))
  for (id in names(cell_rows)) {
    rows <- cell_rows[[id]]; o <- obs_rows[[id]]; pp <- pre[[pattern[[id]]]]
    mu <- setNames(beta[paste0("visit", visits)] + beta["uefm_base_c"] * (cells$uefm_base[rows[1]] - primary$base_mean), visits)
    y[rows] <- mu[pp$mv] + pp$A %*% (o$uefm - mu[pp$ov]) + pp$L %*% rnorm(length(pp$mv))
  }
  y - delta * cells$post12
}
impute <- if (arg$imputation == "mar") impute_mar else impute_locf
estimate <- if (arg$estimate == "trajectory") function(fit) visit_changes(fit, primary$base_mean, visits) else week52_minus_week12

# --- Rubin's rules with Barnard-Rubin df -----------------------------------------------------------------
pool <- function(q, se, df_complete) {
  m <- length(q); Q <- mean(q); U <- mean(se^2); B <- var(q); Tv <- U + (1 + 1 / m) * B
  r <- (1 + 1 / m) * B / U; nu_old <- (m - 1) * (1 + 1 / r)^2
  nu_obs <- (df_complete + 1) / (df_complete + 3) * df_complete * U / Tv
  nu <- 1 / (1 / nu_old + 1 / nu_obs); tc <- qt(0.975, nu)
  data.frame(estimate = Q, se = sqrt(Tv), df = nu, lower = Q - tc * sqrt(Tv), upper = Q + tc * sqrt(Tv),
             p = 2 * pt(-abs(Q / sqrt(Tv)), nu), fraction_missing_info = (1 + 1 / m) * B / Tv, m_used = m)
}

# --- delta sweep ------------------------------------------------------------------------------------------
results <- data.frame(); tipping <- NA
for (delta in deltas) {
  t0 <- Sys.time(); draws <- list()
  for (m in seq_len(M)) {
    imputed <- data.frame(pid = cells$pid, visit = cells$visit, uefm = impute(delta), uefm_base = cells$uefm_base)
    d <- rbind(obs[, c("pid", "visit", "uefm", "uefm_base")], imputed)
    d$pid <- factor(d$pid); d$visit <- factor(d$visit, levels = visits)
    d$uefm_base_c <- d$uefm_base - primary$base_mean
    r <- fit_mmrm(d, primary$covariance, df_method)
    if (r$ok) draws[[length(draws) + 1]] <- estimate(r$fit)
  }
  if (!length(draws)) stop(sprintf("no completed-dataset fit converged under %s at delta %.1f", primary$covariance, delta))
  if (length(draws) < M) cat(sprintf("note: %d of %d completed-dataset fits converged at delta %.1f\n", length(draws), M, delta))
  for (name in unique(draws[[1]]$contrast)) {
    rows <- do.call(rbind, lapply(draws, function(x) x[x$contrast == name, ]))
    results <- rbind(results, cbind(delta = delta, contrast = name, week = rows$week[1],
                                    pool(rows$estimate, rows$se, mean(rows$df))))
  }
  this <- results[results$delta == delta, ]
  cat(sprintf("%s delta %4.1f (%d fits, %.0f s)\n", arg$imputation, delta, length(draws), as.numeric(Sys.time() - t0, units = "secs")))
  if (arg$estimate == "contrast") {
    if (is.na(tipping) && delta > 0 && this$lower <= 0) tipping <- delta
    if (!is.na(tipping) && delta >= tipping + 1) break               # one point past the crossing, then stop
  }
}
write.csv(results, arg$out, row.names = FALSE)
cat("wrote", arg$out, "\n")
