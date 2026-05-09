# R/01_simulate.R
# Synthetic dose-response dataset:
# 3 stressors x 2 species x 7 doses x 6 reps = 252 plants.
# Response = relative root length (% of solvent control), so the curves
# decline from ~100% at low dose to a floor near 5-15% at high dose.

set.seed(20250601)
dir.create("data", showWarnings = FALSE)

species   <- c("Tomato", "Cress")
stressors <- c("NaCl", "PEG", "Zn")                 # mM, % (w/v), uM
dose_grid <- list(
  NaCl = c(0,   25,   50,  100, 200, 400, 800),     # mM
  PEG  = c(0,  2.5,    5,   10,  15,  20,  25),     # % (w/v) PEG-6000
  Zn   = c(0,   25,   50,  100, 250, 500, 1000)     # uM ZnSO4
)
trials <- 1:3                                        # repeated experiments
reps   <- 2                                          # reps per trial (3 x 2 = 6 total)
# trial-to-trial multiplicative drift on ED50 (small batch effect)
trial_drift <- c(`1` = 0.92, `2` = 1.00, `3` = 1.10)

# 4-parameter log-logistic mean function (LL.4 in drc):
#   y = c + (d - c) / (1 + exp(b * (log(x) - log(e))))
ll4 <- function(x, b, c, d, e) {
  out <- ifelse(x <= 0, d, c + (d - c) / (1 + exp(b * (log(x) - log(e)))))
  out
}

# species x stressor "true" parameters (b slope, c lower, d upper, e ED50)
# Story:
#  - Tomato = osmotic-sensitive (lower ED50 for NaCl & PEG)
#  - Cress  = heavy-metal-sensitive (lower ED50 for Zn)
truth <- list(
  Tomato = list(NaCl = c(b = 1.4, c = 8,  d = 100, e = 95),    # mM
                PEG  = c(b = 2.0, c = 6,  d = 100, e = 11),    # %
                Zn   = c(b = 1.6, c = 7,  d = 100, e = 280)),  # uM
  Cress  = list(NaCl = c(b = 1.2, c = 12, d = 100, e = 160),   # mM
                PEG  = c(b = 1.8, c = 8,  d = 100, e = 17),    # %
                Zn   = c(b = 2.0, c = 5,  d = 100, e = 90))    # uM
)

rows <- list()
for (sp in species) {
  for (st in stressors) {
    p <- truth[[sp]][[st]]
    for (tr in trials) {
      e_tr <- p["e"] * trial_drift[as.character(tr)]
      for (d_ in dose_grid[[st]]) {
        mu <- ll4(d_, p["b"], p["c"], p["d"], e_tr)
        # heteroscedastic noise: a bit larger near the inflection
        noise_sd <- 5 + 4 * exp(-((log(max(d_, 1)) - log(e_tr))^2) / 2)
        y <- mu + rnorm(reps, 0, noise_sd)
        y <- pmax(y, 0)
        rows[[length(rows) + 1L]] <- data.frame(
          species  = sp,
          stressor = st,
          trial    = tr,
          dose     = d_,
          rep      = seq_len(reps),
          response = y
        )
      }
    }
  }
}
dat <- do.call(rbind, rows)
dat$species  <- factor(dat$species,  levels = species)
dat$stressor <- factor(dat$stressor, levels = stressors)
dat$trial    <- factor(dat$trial,    levels = trials)

write.csv(dat, "data/dose_response.csv", row.names = FALSE)
saveRDS(dat, "data/dose_response.rds")

# stash the truth for the report (so we can compare fitted ED50 vs true)
saveRDS(truth, "data/truth.rds")

cat(sprintf("Wrote %d rows: %d species x %d stressors x %d trials x 7 doses x %d reps\n",
            nrow(dat), length(species), length(stressors), length(trials), reps))
