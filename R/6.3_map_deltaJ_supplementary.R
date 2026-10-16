#!/usr/bin/env Rscript
# delta |J| maps for all metrics (not used in the manuscript)

library(terra)
library(mapdata)
library(tikzDevice)

compile_tikz <- function(tex_path) {
  old_wd <- setwd(dirname(tex_path))
  on.exit(setwd(old_wd))
  # underscore package breaks \pgfimage lookup
  tx <- readLines(basename(tex_path), warn = FALSE)
  tx <- append(tx, "\\let\\pgfimage\\includegraphics",
               after = max(grep("\\usepackage", tx, fixed = TRUE)))
  writeLines(tx, basename(tex_path))
  system(paste("pdflatex -interaction=nonstopmode", basename(tex_path)),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  setwd(old_wd)
}

BASE_DIR   <- "../figures/1_base/"
IN_SHAPE   <- "../raw_data/alpestris_shapefile/data_0.shp"

RES  <- 1 / 6
CRS  <- "EPSG:4326"
XLIM <- c(-10, 40)
YLIM <- c(35, 70)

METRICS <- c("MaxWgt", "MaxLen", "MinRes",
             "Bsk", "Forage", "Food",
             "Tmax", "Tmin", "MeanDailyDef",
             "LifeSpan", "ReproEnergy", "R0")

source("utils.R")

add_coastline <- function(col = "grey40", lwd = 0.3) {
  maps::map("worldHires", xlim = XLIM, ylim = YLIM,
            col = col, lwd = lwd, interior = FALSE, add = TRUE)
}

add_shapefile <- function(col = "grey20", lwd = 0.3) {
  plot(shp, col = NA, border = col, lwd = lwd, add = TRUE)
}

cat("Loading data\n")
shp <- vect(IN_SHAPE)
template <- rast(xmin = XLIM[1], xmax = XLIM[2],
                 ymin = YLIM[1], ymax = YLIM[2],
                 resolution = RES, crs = CRS)

cmp <- readRDS("../results/jensen_comparison.rds")
per_cell <- cmp$per_cell
# males: supplementary
per_cell <- per_cell[per_cell$sex == "M", ]
cat(sprintf("  Per-cell (male): %d rows\n", nrow(per_cell)))

# delta J rasters
cat("Building ΔJ rasters\n")
dz <- lapply(METRICS, function(m) {
  sub <- per_cell[per_cell$metric == m, ]
  rasterize(cbind(sub$lon, sub$lat), template,
            values = sub$delta_abs_J_w4, fun = mean, na.rm = TRUE)
})
names(dz) <- METRICS

# delta |J| maps
cat("\nDelta |J| maps (all metrics)\n")

pal_div <- colorRampPalette(c(rgb(.10,.10,.44), "white", rgb(.80,.22,.18)))(21)

METRIC_LABS <- c(MaxWgt = "Maximum body mass",
                 MaxLen = "Maximum length",
                 MinRes = "Minimum reserve",
                 Bsk    = "Basking time",
                 Forage = "Foraging time",
                 Food   = "Food intake",
                 Tmax   = expression(italic(T[b])~"max"),
                 Tmin   = expression(italic(T[b])~"min"),
                 MeanDailyDef = "Daily water need (g~day$^{-1}$)",
                 LifeSpan = "Life span",
                 ReproEnergy = "Reproductive energy (J)",
                 R0     = expression(italic(R)[0]))

tikz(paste0(BASE_DIR, "figureS1_deltaJ_all.tex"), width = 9, height = 6.75, standAlone = TRUE)
layout(matrix(1:12, nrow = 3, byrow = TRUE), respect = FALSE)

for (i in seq_along(METRICS)) {
  m <- METRICS[i]
  r <- dz[[m]]
  vals <- terra::values(r)
  is_zero <- all(vals == 0, na.rm = TRUE)

  if (is_zero) {
    plot(r, col = "grey95", legend = FALSE, axes = FALSE,
         mar = c(1, 1, 1, 2))
  } else {
    rng <- compute_limits(vals, robust = FALSE, type = "symmax")
    if (is.null(rng)) {
      plot(r, col = "grey95", legend = FALSE, axes = FALSE,
           mar = c(1, 1, 1, 2))
    } else {
      plot(r, col = pal_div, range = rng,
           legend = TRUE, axes = FALSE, mar = c(1, 1, 1, 2),
           plg = list(cex = 1.8, title = ""))
    }
  }
  add_coastline()
  add_shapefile()
}

dev.off()
compile_tikz(paste0(BASE_DIR, "figureS1_deltaJ_all.tex"))

cat("\nDone.\n")
cat(sprintf("  %sfigureS1_deltaJ_all.tex\n", BASE_DIR))
