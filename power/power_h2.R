# =============================================================================
# H2  Stress exposure and chronic stress markers
# DFG proposal "Sharing stress in the family" — J. Blasberg
#
# Stress exposure is associated with an indicator of chronic stress (hair cortisol,
# perceived stress) between families (stable differences) and within families
# (change across sessions). Per role (parents or adolescents): N = 140 persons,
# 3 sessions, ~15% attrition until T3 (monotone, MCAR), ICC = .50 for both
# variables. Separate tests of the within- and between-family association
# (alpha = .05, two-sided).
#
# Run: Rscript power_h2.R [reps]   -> results/power_h2.csv
# =============================================================================
suppressPackageStartupMessages({ library(lme4); library(data.table); library(future.apply) })
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args)) as.integer(args[1]) else 1000
plan(multisession, workers = 2); set.seed(20261116)

N_REC  <- 140                  # persons per role recruited at T1
RETAIN <- c(1, 0.925, 0.85)    # cumulative retention T1, T2, T3 (MCAR, monotone)
ALPHA  <- 0.05

# monotone dropout: last session attended per person
last_session <- function(n) {
  u <- runif(n)
  ifelse(u > RETAIN[2], 1L, ifelse(u > RETAIN[3], 2L, 3L))
}
pval_z <- function(est, se) 2 * pnorm(-abs(est / se))

# rho = correlation of exposure and indicator, at both the between and the within level
sim_H2 <- function(rho) {
  ls <- last_session(N_REC)
  S <- matrix(c(1, rho, rho, 1), 2)
  B <- MASS::mvrnorm(N_REC, c(0, 0), S) * sqrt(0.5)                 # between-person part
  d <- CJ(id = 1:N_REC, session = 1:3)[session <= ls[id]]
  W <- MASS::mvrnorm(nrow(d), c(0, 0), S) * sqrt(0.5)               # session-specific part
  d[, `:=`(x = B[id, 1] + W[, 1], y = B[id, 2] + W[, 2])]
  d[, `:=`(x_b = mean(x), x_w = x - mean(x)), by = id]              # person-mean centering
  fit <- suppressMessages(lmer(y ~ x_w + x_b + (1 | id), data = d, REML = FALSE))
  cf <- summary(fit)$coefficients
  c(within  = pval_z(cf["x_w", 1], cf["x_w", 2]) < ALPHA,
    between = pval_z(cf["x_b", 1], cf["x_b", 2]) < ALPHA)
}

dir.create("results", showWarnings = FALSE)
res <- rbindlist(lapply(c(.10, .15, .20, .25, .30), function(rho) {
  r <- future_replicate(REPS, sim_H2(rho), future.seed = TRUE)
  data.table(test = "H2 exposure-indicator", value = rho, outcome = rownames(r),
             power = rowMeans(r))
}))
print(res, digits = 3)
fwrite(res, "results/power_h2.csv")
