#!/usr/bin/env Rscript
# Mixed model for repeated measures: Y_it = mu_t + beta * Y_i0 + e_it, follow-up week as an unordered
# factor, baseline UEFM as a covariate, REML, Kenward-Roger degrees of freedom.
# Usage: Rscript R/mmrm_fit.R <input.csv> <output_prefix>
#   input columns: pid, visit, uefm, uefm_base   (one row per observed post-baseline visit)
#   writes <prefix>_lsmeans.csv (LS mean and change from baseline per week), <prefix>_baseline_coef.csv,
#   <prefix>_summary.json (covariance structure, mean baseline, per-week SD used by the imputation draws)
# Rscript passes the script path with spaces encoded as "~+~"
script_dir <- dirname(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE))
source(file.path(script_dir, "mmrm_functions.R"))
args <- commandArgs(trailingOnly = TRUE)
d <- read_long(args[1]); out <- args[2]
base_mean <- attr(d, "base_mean")
visits <- levels(d$visit)

res <- fit_first_acceptable(d)
fit <- res$fit
cat("covariance structures tried:", paste(res$attempts, collapse = "; "), "\n")

# LS mean at each visit (at the mean baseline) and the implied change from baseline.
lsm <- do.call(rbind, lapply(visits, function(v) {
  k <- k_vector(fit); k[paste0("visit", v)] <- 1
  ls <- contrast(fit, k); ch <- contrast(fit, k, null = base_mean)
  data.frame(visit = v, week = week_of(v), n_obs = sum(d$visit == v),
             lsmean = ls$estimate, se = ls$se, df = ls$df, lower = ls$lower, upper = ls$upper,
             change = ch$estimate, change_lower = ch$lower, change_upper = ch$upper, change_p = ch$p)
}))
write.csv(lsm, paste0(out, "_lsmeans.csv"), row.names = FALSE)
write.csv(contrast(fit, k_vector(fit, uefm_base_c = 1)), paste0(out, "_baseline_coef.csv"), row.names = FALSE)

summary <- list(covariance = res$cov, n_obs = nrow(d), n_patients = nlevels(d$pid), base_mean = base_mean,
                visit_sd = setNames(as.list(sqrt(diag(VarCorr(fit)))), visits))
writeLines(toJSON(summary, auto_unbox = TRUE, digits = 8, pretty = TRUE), paste0(out, "_summary.json"))
cat(sprintf("MMRM: covariance %s | %d observations, %d patients | wrote %s_*\n", res$cov, nrow(d), nlevels(d$pid), out))
