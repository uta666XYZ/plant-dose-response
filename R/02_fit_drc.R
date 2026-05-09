# R/02_fit_drc.R
# Fit a 4-parameter log-logistic curve per (species, stressor) using drc,
# extract ED50 with 95% CI, and write the comparison plots.

suppressPackageStartupMessages({
  library(drc); library(dplyr); library(tidyr); library(ggplot2)
})

dir.create("figures", showWarnings = FALSE)
dat <- readRDS("data/dose_response.rds")

# fit one model per stressor with species as the curve identifier
fit_one <- function(stressor_name) {
  d <- dat[dat$stressor == stressor_name, ]
  drm(response ~ dose, curveid = species, data = d,
      fct = LL.4(names = c("slope", "lower", "upper", "ED50")))
}

fits <- lapply(setNames(levels(dat$stressor), levels(dat$stressor)), fit_one)

# ---- ED50 table ----------------------------------------------------------
ed50 <- do.call(rbind, lapply(names(fits), function(st) {
  e <- ED(fits[[st]], 50, interval = "delta", display = FALSE)
  out <- as.data.frame(e)
  out$stressor <- st
  out$species  <- sub(".*:(.*):.*", "\\1", rownames(e))
  out
}))
names(ed50)[1:4] <- c("ED50", "SE", "Lower95", "Upper95")
ed50 <- ed50[, c("stressor", "species", "ED50", "SE", "Lower95", "Upper95")]
rownames(ed50) <- NULL
write.csv(ed50, "figures/ED50_table.csv", row.names = FALSE)

# ---- predicted curves on a fine grid -------------------------------------
grid_one <- function(st) {
  m <- fits[[st]]
  xr <- range(dat$dose[dat$stressor == st & dat$dose > 0])
  xs <- exp(seq(log(xr[1]), log(xr[2]), length.out = 200))
  preds <- predict(m, newdata = expand.grid(dose = xs,
                                            species = levels(dat$species)),
                   interval = "confidence")
  out <- expand.grid(dose = xs, species = levels(dat$species))
  out$fit   <- preds[, 1]
  out$lwr   <- preds[, 2]
  out$upr   <- preds[, 3]
  out$stressor <- st
  out
}
curves <- do.call(rbind, lapply(levels(dat$stressor), grid_one))

# ---- plot 1: dose-response curves ----------------------------------------
units <- c(NaCl = "(mM)", PEG = "(% w/v)", Zn = "(uM)")
curves$panel <- sprintf("%s %s", curves$stressor, units[as.character(curves$stressor)])
dat$panel    <- sprintf("%s %s", dat$stressor,    units[as.character(dat$stressor)])

p_curves <- ggplot() +
  geom_ribbon(data = curves,
              aes(dose, ymin = lwr, ymax = upr, fill = species),
              alpha = 0.18) +
  geom_line(data = curves,
            aes(dose, fit, colour = species), linewidth = 0.8) +
  geom_jitter(data = dat[dat$dose > 0, ],
              aes(dose, response, colour = species),
              width = 0.05, height = 0, alpha = 0.7, size = 1.4) +
  facet_wrap(~ panel, scales = "free_x") +
  scale_x_log10() +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Dose (log scale)",
       y = "Response (% of control)",
       title = "Dose-response curves (4-parameter log-logistic, drc::drm)",
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave("figures/curves.png", p_curves, width = 9, height = 4, dpi = 150)

# ---- plot 2: ED50 comparison with 95% CI ---------------------------------
p_ed50 <- ggplot(ed50, aes(species, ED50, fill = species)) +
  geom_col(width = 0.6, colour = "grey25") +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95),
                width = 0.2, linewidth = 0.6) +
  facet_wrap(~ stressor, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = "ED50 (dose units of each stressor)",
       title = "Sensitivity comparison — ED50 with 95% CI") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

ggsave("figures/ed50_bars.png", p_ed50, width = 9, height = 3.6, dpi = 150)

# expose objects
assign("fits",   fits,   envir = globalenv())
assign("ed50",   ed50,   envir = globalenv())
assign("curves", curves, envir = globalenv())

cat("Wrote figures/curves.png, figures/ed50_bars.png, figures/ED50_table.csv\n")
