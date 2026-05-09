# Plant dose-response — drc fits & ED50 comparison

Portfolio piece — a synthetic recreation of a dose-response analysis I
helped run during a colleague's MSc project. Plate-level **% of
control** response across a range of doses and three stressors
(NaCl, Cu, glyphosate); four-parameter log-logistic fits with
**`drc::drm()`**; ED50 (= IC50) extraction with delta-method 95 % CIs;
species-sensitivity comparison.

[![Live report](https://img.shields.io/badge/Live-report-blue?logo=github)](https://uta666xyz.github.io/plant-dose-response/)

The original measurements are not shareable, so the dataset is
generated from log-logistic curves with hand-chosen parameters plus
Gaussian noise. The qualitative pattern (Cress more sensitive to
glyphosate, Lettuce more sensitive to Cu) follows what was observed
on the real plates.

## What's here

* `R/01_simulate.R` — generates `data/dose_response.{csv,rds}`
  (2 species × 3 stressors × 7 doses × 6 reps = 252 plants).
* `R/02_fit_drc.R` — fits `LL.4()` per (species × stressor),
  writes `figures/curves.png`, `figures/ed50_bars.png`, and
  `figures/ED50_table.csv`.
* `analysis.Rmd` — knit to a self-contained HTML report containing
  the curve panel, the ED50 table & bar chart, a sanity check against
  the simulator's ground truth, and pairwise ED50-ratio tests
  (`drc::EDcomp`).

## Reproducing

```r
install.packages(c("drc", "ggplot2", "dplyr", "tidyr", "rmarkdown"))
source("R/01_simulate.R")            # writes data/*.{csv,rds}
source("R/02_fit_drc.R")             # writes figures/*.png + ED50 table
rmarkdown::render("analysis.Rmd")    # writes analysis.html
```

## What the figures show

* **Dose-response curves** — points are individual replicates, lines
  are `LL.4()` fits, ribbons are 95 % confidence bands. X-axis on
  log scale.
* **ED50 bar chart** — ED50 with delta-method 95 % CIs per
  (species × stressor); lower bars mean greater sensitivity.
* **Truth comparison table** — fitted vs simulator-true ED50 with
  relative error (sanity check).

## Disclaimer

Numbers are simulated. The repo demonstrates the **drc-based
dose-response workflow** (fit → ED50 → CI → species comparison),
not a finding about real chemicals or species.
