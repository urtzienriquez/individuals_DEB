#!/usr/bin/env Rscript
# Jensen bias amplification under +4°C for selected metrics (superseded by 7.1; not used in the manuscript)

library(terra)
library(mapdata)
library(tikzDevice)

BASE_DIR   <- "../figures/1_base/"
IN_SHAPE   <- "../raw_data/alpestris_shapefile/data_0.shp"

RES  <- 1 / 6
CRS  <- "EPSG:4326"
XLIM <- c(-10, 40)
YLIM <- c(35, 70)

SELECTED   <- c("Bsk", "Forage", "R0", "Fecundity")
SELECTED_M <- c("Bsk", "Forage", "R0", "ReproEnergy")

source("utils.R")

METRIC_LABS <- c(MaxWgt = "Max body mass",
                 MaxLen = "Max length",
                 MinRes = "Min reserve",
                 Bsk    = "Basking time",
                 Forage = "Foraging time",
                 Food   = "Food intake",
                 Tmax   = expression(italic(T[b])~"max"),
                 Tmin   = expression(italic(T[b])~"min"),
                 MeanDailyDef = "Daily water need (g~day$^{-1}$)",
                 LifeSpan = "Life span",
                 Fecundity = "Fecundity",
                 Clutches = "Clutches",
                 ReproEnergy = "Reproductive energy (J)",
                 R0     = expression(italic(R)[0]))

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

add_coastline <- function(col = "grey40", lwd = 0.3) {
  maps::map("worldHires", xlim = XLIM, ylim = YLIM,
            col = col, lwd = lwd, interior = FALSE, add = TRUE)
}

add_shapefile <- function(col = "grey20", lwd = 0.3) {
  plot(shp, col = NA, border = col, lwd = lwd, add = TRUE)
}

rasterize_deltaJ <- function(df, metric_col, template) {
  sub <- df[df$metric == metric_col, ]
  rasterize(cbind(sub$lon, sub$lat), template,
            values = sub$delta_abs_J_w4, fun = mean, na.rm = TRUE)
}

cat("Loading data\n")
shp <- vect(IN_SHAPE)
template <- rast(xmin = XLIM[1], xmax = XLIM[2],
                 ymin = YLIM[1], ymax = YLIM[2],
                 resolution = RES, crs = CRS)

cmp <- readRDS("../results/jensen_comparison.rds")
stats <- cmp$per_metric
per_cell <- cmp$per_cell
cat(sprintf("  Per-metric stats: %d rows\n", nrow(stats)))
cat(sprintf("  Per-cell: %d rows\n", nrow(per_cell)))

# females: main text
per_cell_F <- per_cell[per_cell$sex == "F", ]
stats_F    <- stats[stats$sex == "F", ]
# males: supplementary
per_cell_M <- per_cell[per_cell$sex == "M", ]
stats_M    <- stats[stats$sex == "M", ]

# delta J rasters
cat("Building ΔJ rasters\n")
dz <- lapply(SELECTED, function(m) rasterize_deltaJ(per_cell_F, m, template))
names(dz) <- SELECTED

# relative delta J (%)
stats_F$rel_deltaJ <- stats_F$mean_delta_abs_J_w4 / (stats_F$mean_absJ_cur + 1e-10) * 100

cat("Relative ΔJ values (%, sorted):\n")
for (i in order(stats_F$rel_deltaJ)) {
  cat(sprintf("  %-12s %7.1f%% (pct_inc = %.0f%%)\n",
      stats_F$metric[i], stats_F$rel_deltaJ[i], stats_F$pct_increase_w4[i]))
}

# females
cat("\nJensen amplification maps (female)\n")

pal_div <- colorRampPalette(c("#440154", "white", "#FDE725"))(21)

bar_cols <- ifelse(stats_F$rel_deltaJ > 0, "#FDE725", "#440154")

tikz(paste0(BASE_DIR, "figure5_jensen_comparison.tex"), width = 8, height = 6, standAlone = TRUE)

layout(matrix(c(1, 1, 1, 1,   # a1
                2, 2, 2, 2,   # a2
                3, 4, 5, 6), nrow = 3, byrow = TRUE),
       heights = c(1.1, 0.9, 1.5))

# panel a1
par(mar = c(1, 5, 3, 1))
ord <- order(stats_F$rel_deltaJ)
bar_vals <- stats_F$rel_deltaJ[ord]
x <- barplot(bar_vals, col = bar_cols[ord],
             border = NA, axes = FALSE,
             xlab = NULL, ylab = "",
             names.arg = NULL,
             ylim = c(min(bar_vals) * 1.5, max(bar_vals) * 1.25))
ticks <- pretty(par("usr")[3:4])
axis(2, at = ticks, las = 1, cex.axis = 1.3)
mtext("$\\Delta|J|$ (\\%)", side = 2, line = 2.5, cex = 1.1)
abline(h = 0, lwd = 0.8)

# panel a2
par(mar = c(5.5, 5, 1, 1))
pct_vals <- stats_F$pct_increase_w4[ord]
x <- barplot(pct_vals, col = bar_cols[ord],
             border = NA, axes = FALSE,
             xlab = NULL, ylab = "",
             names.arg = NULL,
             ylim = c(0, 100))
x_labs <- sub("^Tmax$", "italic(T[b])~max",
              sub("^Tmin$", "italic(T[b])~min",
                  sub("^R0$", "italic(R)[0]",
                      stats_F$metric[ord])))
text(x, par("usr")[3] - diff(par("usr")[3:4]) * 0.065,
     labels = parse(text = x_labs), srt = 35, adj = 1, xpd = TRUE, cex = 1.1)
axis(2, at = seq(0, 100, 25), las = 1, cex.axis = 1.3)
mtext("Cells with increased $|J|$ (\\%)", side = 2, line = 2.5, cex = 1.1)
abline(h = 50, lwd = 0.5, lty = 3)

# panel b
for (i in seq_along(SELECTED)) {
  m <- SELECTED[i]
  r <- dz[[m]]
  vals <- terra::values(r)
  is_zero <- all(vals == 0, na.rm = TRUE)

  if (is_zero) {
    plot(r, col = "grey95", legend = FALSE, axes = FALSE, mar = c(3, 1, 1.5, 2))
  } else {
    rng <- compute_limits(vals, robust = TRUE, type = "symmax")
    if (is.null(rng)) {
      plot(r, col = "grey95", legend = FALSE, axes = FALSE, mar = c(3, 1, 1.5, 2))
    } else {
      plot(r, col = pal_div, range = rng,
           legend = TRUE, axes = FALSE, mar = c(3, 1, 1.5, 2),
           plg = list(cex = 1.4, title = "", x = "bottom"))
    }
  }
  add_coastline()
  add_shapefile()
}

dev.off()

compile_tikz(paste0(BASE_DIR, "figure5_jensen_comparison.tex"))

# males
cat("\nSupplementary: male Jensen amplification\n")

dz_M <- lapply(SELECTED_M, function(m) rasterize_deltaJ(per_cell_M, m, template))
names(dz_M) <- SELECTED_M

stats_M$rel_deltaJ <- stats_M$mean_delta_abs_J_w4 / (stats_M$mean_absJ_cur + 1e-10) * 100
bar_cols_M <- ifelse(stats_M$rel_deltaJ > 0, "#FDE725", "#440154")

tikz(paste0(BASE_DIR, "figureS5_jensen_comparison.tex"), width = 8, height = 6, standAlone = TRUE)

layout(matrix(c(1, 1, 1, 1,
                2, 2, 2, 2,
                3, 4, 5, 6), nrow = 3, byrow = TRUE),
       heights = c(1.1, 0.9, 1.5))

# panel a1
par(mar = c(1, 5, 3, 1))
ord_M <- order(stats_M$rel_deltaJ)
bar_vals_M <- stats_M$rel_deltaJ[ord_M]
x <- barplot(bar_vals_M, col = bar_cols_M[ord_M],
             border = NA, axes = FALSE,
             xlab = NULL, ylab = "",
             names.arg = NULL,
             ylim = c(min(bar_vals_M) * 1.5, max(bar_vals_M) * 1.25))
ticks <- pretty(par("usr")[3:4])
axis(2, at = ticks, las = 1, cex.axis = 1.3)
mtext("$\\Delta|J|$ (\\%)", side = 2, line = 2.5, cex = 1.1)
abline(h = 0, lwd = 0.8)

# panel a2
par(mar = c(5.5, 5, 1, 1))
pct_vals_M <- stats_M$pct_increase_w4[ord_M]
x <- barplot(pct_vals_M, col = bar_cols_M[ord_M],
             border = NA, axes = FALSE,
             xlab = NULL, ylab = "",
             names.arg = NULL,
             ylim = c(0, 100))
x_labs_M <- sub("^Tmax$", "italic(T[b])~max",
               sub("^Tmin$", "italic(T[b])~min",
                   sub("^R0$", "italic(R)[0]",
                       stats_M$metric[ord_M])))
text(x, par("usr")[3] - diff(par("usr")[3:4]) * 0.065,
     labels = parse(text = x_labs_M), srt = 35, adj = 1, xpd = TRUE, cex = 1.1)
axis(2, at = seq(0, 100, 25), las = 1, cex.axis = 1.3)
mtext("Cells with increased $|J|$ (\\%)", side = 2, line = 2.5, cex = 1.1)
abline(h = 50, lwd = 0.5, lty = 3)

# panel b
for (i in seq_along(SELECTED_M)) {
  m <- SELECTED_M[i]
  r <- dz_M[[m]]
  vals <- terra::values(r)
  is_zero <- all(vals == 0, na.rm = TRUE)

  if (is_zero) {
    plot(r, col = "grey95", legend = FALSE, axes = FALSE, mar = c(3, 1, 1.5, 2))
  } else {
    rng <- compute_limits(vals, robust = TRUE, type = "symmax")
    if (is.null(rng)) {
      plot(r, col = "grey95", legend = FALSE, axes = FALSE, mar = c(3, 1, 1.5, 2))
    } else {
      plot(r, col = pal_div, range = rng,
           legend = TRUE, axes = FALSE, mar = c(3, 1, 1.5, 2),
           plg = list(cex = 1.4, title = "", x = "bottom"))
    }
  }
  add_coastline()
  add_shapefile()
}

dev.off()
compile_tikz(paste0(BASE_DIR, "figureS5_jensen_comparison.tex"))

cat("\nDone.\n")
cat(sprintf("  %sfigure5_jensen_comparison.pdf (female)\n", BASE_DIR))
cat(sprintf("  %sfigureS5_jensen_comparison.pdf (male)\n", BASE_DIR))
