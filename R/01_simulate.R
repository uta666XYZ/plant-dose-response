# R/01_simulate.R
# Synthetic dose-response dataset:
# 3 stressors x 2 species x 7 doses x 6 reps = 252 plants.
# Response = relative root length (% of solvent control), so the curves
# decline from ~100% at low dose to a floor near 5-15% at high dose.

set.seed(20240601)
dir.create("data", showWarnings = FALSE)

species   <- c("Lettuce", "Cress")
stressors <- c("NaCl", "Cu", "Glyphosate")          # mM, uM, uM
dose_grid <- list(
  NaCl       = c(0, 25,  50, 100, 200, 400, 800),    # mM
  Cu         = c(0,  5,  10,  25,  50, 100, 200),    # uM
  Glyphosate = c(0,  1,   3,  10,  30, 100, 300)     # uM
)
reps <- 6

# 4-parameter log-logistic mean function (LL.4 in drc):
#   y = c + (d - c) / (1 + exp(b * (log(x) - log(e))))
ll4 <- function(x, b, c, d, e) {
  out <- ifelse(x <= 0, d, c + (d - c) / (1 + exp(b * (log(x) - log(e)))))
  out
}

# species x stressor "true" parameters (b slope, c lower, d upper, e ED50)
truth <- list(
  Lettuce = list(NaCl       = c(b = 1.4, c = 8,  d = 100, e = 110),
                 Cu         = c(b = 2.0, c = 5,  d = 100, e = 28),
                 Glyphosate = c(b = 1.8, c = 6,  d = 100, e = 18)),
  Cress   = list(NaCl       = c(b = 1.2, c = 12, d = 100, e = 165),
                 Cu         = c(b = 1.6, c = 8,  d = 100, e = 45),
                 Glyphosate = c(b = 2.2, c = 4,  d = 100, e = 9))
)

rows <- list()
for (sp in species) {
  for (st in stressors) {
    p <- truth[[sp]][[st]]
    for (d_ in dose_grid[[st]]) {
      mu <- ll4(d_, p["b"], p["c"], p["d"], p["e"])
      # heteroscedastic noise: a bit larger near the inflection
      noise_sd <- 5 + 4 * exp(-((log(max(d_, 1)) - log(p["e"]))^2) / 2)
      y <- mu + rnorm(reps, 0, noise_sd)
      y <- pmax(y, 0)
      rows[[length(rows) + 1L]] <- data.frame(
        species  = sp,
        stressor = st,
        dose     = d_,
        rep      = seq_len(reps),
        response = y
      )
    }
  }
}
dat <- do.call(rbind, rows)
dat$species  <- factor(dat$species,  levels = species)
dat$stressor <- factor(dat$stressor, levels = stressors)

write.csv(dat, "data/dose_response.csv", row.names = FALSE)
saveRDS(dat, "data/dose_response.rds")

# stash the truth for the report (so we can compare fitted ED50 vs true)
saveRDS(truth, "data/truth.rds")

cat(sprintf("Wrote %d rows across %d species x %d stressors x 7 doses x %d reps\n",
            nrow(dat), length(species), length(stressors), reps))
