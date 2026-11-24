# =============================================================================
# H1 (ii) pooled across both directions, H4a and H4b
# DFG proposal "Sharing stress in the family" — J. Blasberg
#
# (1) H1 (ii) pooled across both directions of transmission (direction as moderator),
#     using the assumptions and simulate() from power_paradigm.R
# (2) H4a / H4b across three sessions with transmission strength = dyad-level listener
#     reactivity (reliability ~ .63 from power_paradigm.R), including a variant in which
#     more burdened families drop out more often (MAR dropout).
#
# Run: Rscript power_h1ii_h4.R [reps]   -> results/power_h1ii_h4.csv
#      (requires power_paradigm.R in the same folder)
# =============================================================================
suppressPackageStartupMessages({ library(lme4); library(data.table); library(future.apply) })
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args)) as.integer(args[1]) else 500

src <- readLines("power_paradigm.R")
eval(parse(text = src[1:(grep("^one_rep <- function", src) - 1)]))   # BASE, simulate()
plan(multisession, workers = 2); set.seed(20261118)   # after sourcing, which sets its own seed

# ---- (1) H1 (ii) pooled ------------------------------------------------------
h1ii_pooled <- function(beta) {
  P <- BASE; P$beta <- beta
  a <- simulate(P)$d[stress == 1][, direction := -0.5]               # parent -> adolescent
  b <- simulate(P)$d[stress == 1][, direction :=  0.5]               # adolescent -> parent
  ds <- rbind(a, b)
  ds[, `:=`(nw = narr_obs - mean(narr_obs), nb = mean(narr_obs)), by = .(dyad, direction)]
  # listener-specific intercepts (parent and adolescent differ in HR level); slope varies by dyad
  f <- suppressMessages(lmer(listener ~ nw * direction + nb + (0 + nw | dyad) + (1 | dyad:direction) +
                               (1 | dyad:direction:session), data = ds, REML = FALSE))
  cf <- summary(f)$coefficients["nw", ]
  unname(abs(cf[1] / cf[2]) > 1.96)
}

# ---- (2) H4a / H4b, three sessions -------------------------------------------
N <- 140; RETAIN <- c(1, .925, .85)
dropout <- function(score, mar) {
  # 15% drop out by T3 (half after T1, half after T2); under MAR, probability ~ exp(0.7 * score)
  w <- if (mar) exp(0.7 * score) else rep(1, length(score))
  drop <- sample(seq_along(score), round(.15 * length(score)), prob = w)
  last <- rep(3L, length(score)); last[drop] <- sample(1:2, length(drop), replace = TRUE)
  last
}
gen_h4 <- function(p, i = 0, mar = FALSE, actor = .30, r_x = .30, rel_x = .80, rel_m = .63, m = .10) {
  S <- matrix(c(1, r_x, r_x, 1), 2)
  Xb <- MASS::mvrnorm(N, c(0, 0), S) * sqrt(.5)
  last <- dropout(rowMeans(Xb), mar)
  d <- CJ(dyad = 1:N, session = 1:3)[session <= last[dyad]]
  Xw <- MASS::mvrnorm(nrow(d), c(0, 0), S) * sqrt(.5)
  Xa <- Xb[d$dyad, 1] + Xw[, 1]; Xp <- Xb[d$dyad, 2] + Xw[, 2]
  M_true <- rnorm(N)                                        # dyad-level transmission strength
  M_obs  <- sqrt(rel_m) * M_true + sqrt(1 - rel_m) * rnorm(N)
  M <- M_true[d$dyad]
  lp <- actor * Xa + p * Xp + m * M + i * Xp * M
  ev <- max(1 - var(lp), .2)
  y  <- lp + rnorm(N, 0, sqrt(ev / 2))[d$dyad] + rnorm(nrow(d), 0, sqrt(ev / 2))
  d[, `:=`(y = y, xa = sqrt(rel_x) * Xa + sqrt(1 - rel_x) * rnorm(.N),
           xp = sqrt(rel_x) * Xp + sqrt(1 - rel_x) * rnorm(.N), mo = M_obs[dyad])]
  d
}
pz <- function(cf, nm) 2 * pnorm(-abs(cf[nm, 1] / cf[nm, 2]))
h4a <- function(p, mar) {
  cf <- summary(suppressMessages(lmer(y ~ xa + xp + (1 | dyad), data = gen_h4(p, mar = mar), REML = FALSE)))$coefficients
  pv <- pz(cf, "xp"); c(a05 = pv < .05, a005 = pv < .005)
}
h4b <- function(i, mar) {
  cf <- summary(suppressMessages(lmer(y ~ xa + xp * mo + (1 | dyad), data = gen_h4(.20, i = i, mar = mar), REML = FALSE)))$coefficients
  pz(cf, "xp:mo") < .05
}

out <- list()
for (b in c(.10, .125, .15, .20)) out[[length(out) + 1]] <- data.table(test = "H1 (ii) pooled", value = b, variant = "MCAR", outcome = "power",
  power = mean(future_replicate(REPS, h1ii_pooled(b), future.seed = TRUE)))
for (mar in c(FALSE, TRUE)) {
  for (p in c(.10, .15, .20)) {
    r <- future_replicate(REPS, h4a(p, mar), future.seed = TRUE)
    out[[length(out) + 1]] <- data.table(test = "H4a partner effect", value = p, variant = ifelse(mar, "MAR", "MCAR"),
                                         outcome = c("alpha .05", "alpha .005"), power = rowMeans(r))
  }
  for (i in c(.10, .15, .20, .25)) out[[length(out) + 1]] <- data.table(test = "H4b moderation", value = i,
    variant = ifelse(mar, "MAR", "MCAR"), outcome = "power", power = mean(future_replicate(REPS, h4b(i, mar), future.seed = TRUE)))
}
res <- rbindlist(out); print(res, digits = 2)
fwrite(res, "results/power_h1ii_h4.csv")
