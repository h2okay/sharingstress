# =============================================================================
# Power and reliability for the family stress conversation (H1, transmission indices)
# DFG proposal "Sharing stress in the family" — J. Blasberg, 17.11.2026
#
# One direction of transmission (e.g., parent -> adolescent); the other is identical.
# Per session: 3 stressor trials (narrator = parent, listener = adolescent) and
# 2 comparison trials (neutral, positive); 3 sessions; 140 dyads, ~15% attrition to T3.
# Outcome: listener's mean heart rate in the narration phase (bpm).
#
# Anchor: in Blasberg, Kanske, & Engert (2026), adolescents' HR increased by 8.1 bpm (SD 10.7)
# when interacting with their stressed parent (TSST, baseline-to-peak). We assume about
# one fifth of this mean effect for recounted stress and mean (not peak) heart rate.
#
# Run: Rscript power_paradigm.R [reps]   -> results/power_paradigm.csv
# =============================================================================
suppressPackageStartupMessages({ library(lme4); library(data.table); library(future.apply) })
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args)) as.integer(args[1]) else 500
plan(multisession, workers = 2); set.seed(20261117)

# ---- assumptions (bpm) -------------------------------------------------------
BASE <- list(
  delta    = 1.5,   # listener HR: stressor vs comparison trials (H1 i)
  sd_react = 3,     # between-dyad SD of the listener's reactivity
  sd_trial = 4,     # trial-to-trial SD within a session (residual)
  sd_sess  = 3,     # session-to-session SD of the listener's HR level
  mu_narr  = 4,     # narrator's mean HR response in stressor trials
  sd_narr  = 3,     # trial-to-trial SD of the narrator's response within a dyad
  me_narr  = 3,     # measurement noise in the narrator's observed response
  beta     = 0.15,  # within-dyad slope: listener response ~ narrator response (H1 ii)
  sd_beta  = 0.15,  # between-dyad SD of this slope (transmission strength)
  r_h3     = 0.30,  # true correlation of transmission strength with chronic family stress (H3)
  n_dyads  = 140, retain = c(1, .925, .85)
)

simulate <- function(P) {
  n <- P$n_dyads
  u <- runif(n); last <- ifelse(u > P$retain[2], 1L, ifelse(u > P$retain[3], 2L, 3L))
  d <- CJ(dyad = 1:n, session = 1:3, trial = 1:5)[session <= last[dyad]]
  d[, stress := as.integer(trial <= 3)]
  react <- rnorm(n, P$delta, P$sd_react)                 # dyad-specific listener reactivity
  b_d   <- rnorm(n, P$beta, P$sd_beta)                   # dyad-specific transmission slope
  lvl   <- rnorm(n, 75, 8); sess <- rnorm(n * 3, 0, P$sd_sess)
  d[, narr_true := ifelse(stress == 1, rnorm(.N, P$mu_narr, P$sd_narr), 0)]
  d[, narr_obs  := narr_true + rnorm(.N, 0, P$me_narr)]
  d[, narr_c    := narr_true - mean(narr_true[stress == 1]), by = dyad]   # within-dyad deviation
  d[, listener  := lvl[dyad] + sess[(dyad - 1) * 3 + session] + stress * react[dyad] +
                   stress * b_d[dyad] * narr_c + rnorm(.N, 0, P$sd_trial)]
  list(d = d, react = react, b_d = b_d)
}

one_rep <- function(P) {
  s <- simulate(P); d <- s$d
  # H1 (i): stressor vs comparison trials
  f1 <- suppressMessages(lmer(listener ~ stress + (1 + stress || dyad) + (1 | dyad:session),
                              data = d, REML = FALSE))
  c1 <- summary(f1)$coefficients["stress", ]
  # H1 (ii): listener response increases with the narrator's (stressor trials, within-dyad)
  ds <- d[stress == 1]
  ds[, `:=`(nw = narr_obs - mean(narr_obs), nb = mean(narr_obs)), by = dyad]
  f2 <- suppressMessages(lmer(listener ~ nw + nb + (1 + nw | dyad) + (1 | dyad:session),
                              data = ds, REML = FALSE))
  c2 <- summary(f2)$coefficients["nw", ]
  # reliability of dyad-level indices (observed, two-step): correlation^2 with true values
  idx <- d[, .(react_hat = mean(listener[stress == 1]) - mean(listener[stress == 0])), by = dyad]
  sl  <- ds[, .(slope_hat = if (.N >= 3 && var(nw) > 0) coef(lm(listener ~ nw))[2] else NA_real_), by = dyad]
  # H3: family stress correlates r_h3 with the TRUE index; tested with the observed index
  z <- function(x) (x - mean(x)) / sd(x)
  fam_r <- P$r_h3 * z(s$react) + sqrt(1 - P$r_h3^2) * rnorm(P$n_dyads)
  fam_b <- P$r_h3 * z(s$b_d)   + sqrt(1 - P$r_h3^2) * rnorm(P$n_dyads)
  p_r <- cor.test(idx$react_hat, fam_r[idx$dyad])$p.value
  p_b <- cor.test(sl$slope_hat, fam_b[sl$dyad])$p.value
  c(h1_i  = unname(abs(c1[1] / c1[2]) > 1.96),
    h1_ii = unname(abs(c2[1] / c2[2]) > 1.96),
    rel_react = cor(idx$react_hat, s$react[idx$dyad])^2,
    rel_slope = cor(sl$slope_hat, s$b_d[sl$dyad], use = "complete.obs")^2,
    h3_react = p_r < .05, h3_slope = p_b < .05)
}

run <- function(label, param, grid, reps = REPS) {
  rbindlist(lapply(grid, function(v) {
    P <- BASE; P[[param]] <- v
    r <- future_replicate(reps, one_rep(P), simplify = "array", future.seed = TRUE)
    data.table(condition = label, param = param, value = v,
               power_H1_i = mean(r["h1_i", ]), power_H1_ii = mean(r["h1_ii", ]),
               rel_reactivity = mean(r["rel_react", ]), rel_slope = mean(r["rel_slope", ]),
               power_H3_reactivity = mean(r["h3_react", ]), power_H3_slope = mean(r["h3_slope", ]))
  }))
}

dir.create("results", showWarnings = FALSE)
res <- rbindlist(list(
  run("H1 (i): effect size",            "delta",    c(0.5, 1, 1.5, 2)),
  run("H1 (ii): slope",                 "beta",     c(0.05, 0.10, 0.15, 0.20)),
  run("Sensitivity: trial noise",       "sd_trial", c(3, 5)),
  run("Sensitivity: reactivity SD",     "sd_react", c(2, 5)),
  run("Sensitivity: slope SD",          "sd_beta",  c(0.10, 0.30)),
  run("H3: true correlation",           "r_h3",     c(0.20, 0.30, 0.40))
))
print(res, digits = 2)
fwrite(res, "results/power_paradigm.csv")
