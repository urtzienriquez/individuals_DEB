#!/usr/bin/env Rscript
# figs. 4 and S3.5-S3.7: population-average vs individual-level predictions

library(terra)
library(mapdata)
library(tikzDevice)

compile_tikz <- function(tex_path) {
  old_wd <- setwd(dirname(tex_path))
  on.exit(setwd(old_wd))
  tx <- readLines(basename(tex_path), warn = FALSE)
  tx <- append(tx, "\\let\\pgfimage\\includegraphics",
               after = max(grep("\\usepackage", tx, fixed = TRUE)))
  engine <- if (nzchar(Sys.which("lualatex"))) "lualatex" else "pdflatex"
  writeLines(tx, basename(tex_path))
  system(paste(engine, "-interaction=nonstopmode", basename(tex_path)),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
}

BASE_DIR   <- "../figures/1_base/"
IN_SHAPE   <- "../raw_data/alpestris_shapefile/data_0.shp"

RES  <- 1 / 6
CRS  <- "EPSG:4326"
XLIM <- c(-10, 40)
YLIM <- c(35, 70)

METRICS <- c("LifeSpan", "MaxWgt", "Fecundity")

METRICS_M <- c("LifeSpan", "MaxWgt", "ReproEnergy")

source("utils.R")

PAL_CUR      <- colorRampPalette(c("white", "darkgreen"))(100)
PAL_DIVERGING <- colorRampPalette(c("#b2182b", "white", "#2166ac"))(100)

COL_NAMES <- c("P_species_cur", "P_ind_cur",
               "P_species_w2",  "P_ind_w2",
               "P_species_w4",  "P_ind_w4")
DELTA_REF <- c(NA, NA, 1, 2, 1, 2)

rasterize_metric <- function(df, template) {
  rasterize(cbind(df$lon, df$lat), template,
            values = df$value, fun = mean, na.rm = TRUE)
}

add_coastline <- function(col = "grey40", lwd = 0.3) {
  maps::map("worldHires", xlim = XLIM, ylim = YLIM,
            col = col, lwd = lwd, interior = FALSE, add = TRUE)
}

add_shapefile <- function(col = "grey20", lwd = 0.3) {
  if (!is.null(shp)) {
    plot(shp, col = NA, border = col, lwd = lwd, add = TRUE)
  }
}

cat("Loading data\n")
shp <- NULL
if (file.exists(IN_SHAPE)) {
  shp <- vect(IN_SHAPE)
  shp <- crop(shp, ext(XLIM[1], XLIM[2], YLIM[1], YLIM[2]))
  cat("  shapefile loaded\n")
}

template <- rast(xmin = XLIM[1], xmax = XLIM[2],
                 ymin = YLIM[1], ymax = YLIM[2],
                 resolution = RES, crs = CRS)

cmp <- readRDS("../results/jensen_comparison.rds")
per_cell <- cmp$per_cell
cat(sprintf("  Comparison data: %d rows\n", nrow(per_cell)))

# females: main text
per_cell_F <- per_cell[per_cell$sex == "F", ]
# males: supplementary
per_cell_M <- per_cell[per_cell$sex == "M", ]

cat("Building rasters\n")

rasters <- list()
ranges_raw <- list()

for (m in METRICS) {
  rasters[[m]] <- list()
  sub <- per_cell_F[per_cell_F$metric == m, ]
  cat(sprintf("  %s: %d cells\n", m, nrow(sub)))

  for (j in seq_along(COL_NAMES)) {
    df <- sub[, c("lon", "lat")]
    df$value <- sub[[COL_NAMES[j]]]
    rasters[[m]][[j]] <- rasterize_metric(df, template)
  }

  all_vals <- unlist(lapply(rasters[[m]], function(r) terra::values(r)))
  rng <- compute_limits(all_vals, robust = FALSE, type = "range")
  if (is.null(rng)) rng <- c(0, 1)
  ranges_raw[[m]] <- rng

  cat(sprintf("    shared range: [%.4f, %.4f]\n", rng[1], rng[2]))
}

max_abs <- list()
for (m in METRICS) {
  delta_vals <- c()
  for (j in 3:6) {
    delta_vals <- c(delta_vals,
      as.vector(rasters[[m]][[j]] - rasters[[m]][[DELTA_REF[j]]]))
  }
  rng <- compute_limits(delta_vals, robust = FALSE, type = "symmax")
  max_abs[[m]] <- if (is.null(rng)) 1 else max(abs(rng))
}

# fig. 4 (female)
cat("\nFig. 4: population-average vs individual predictions (female)\n")

tikz(paste0(BASE_DIR, "fig4-species-vs-ind-female.tex"),
     width = 12, height = 12, standAlone = TRUE)

lmat <- matrix(1:42, nrow = 6, ncol = 7, byrow = TRUE)
layout(lmat, widths = c(rep(1, 6), 0.4), heights = rep(1, 6))

draw_legend <- function(pal, zrange, is_diverging, m) {
  par(mar = c(0.5, 0.4, 0.5, 0.6))
  plot.new()
  plot.window(xlim = c(0, 0.40), ylim = range(zrange))
  if (is_diverging) {
    z <- seq(-1, 1, length.out = 100)
    lab_at  <- c(-1, 0, 1)
    lab_lbl <- c("-1", "0", "+1")
  } else {
    z <- seq(zrange[1], zrange[2], length.out = 100)
    lab_at <- pretty(zrange, 4)
    lab_at <- lab_at[lab_at >= zrange[1] & lab_at <= zrange[2]]
    lab_lbl <- format(lab_at, digits = 2, trim = TRUE)
  }
  ncol <- length(pal)
  idx <- as.integer(cut(z, breaks = ncol, include.lowest = TRUE))
  cols <- pal[idx]
  y <- seq(zrange[1], zrange[2], length.out = 100)
  rect(rep(0, 99), y[1:99], rep(0.1, 99), y[2:100],
       col = cols[1:99], border = NA)
  rect(0, zrange[1], 0.1, zrange[2], border = "grey60", lwd = 0.3, col = NA)
  axis(4, at = lab_at, labels = lab_lbl, las = 1, cex.axis = 1.0,
       tcl = -0.2, mgp = c(0, 0.4, 0), pos = 0.1)
}

par(oma = c(0, 0, 0, 0))
par(cex = 1)

for (k in 1:42) {
  ri <- (k - 1) %/% 7 + 1
  ci <- (k - 1) %% 7 + 1
  metric_idx <- (ri + 1) %/% 2
  is_bottom  <- ri %% 2 == 0
  m <- METRICS[metric_idx]

  if (ci == 7) {
    if (is_bottom) {
      draw_legend(PAL_DIVERGING, c(-1, 1), is_diverging = TRUE, m)
    } else {
      draw_legend(PAL_CUR, ranges_raw[[m]], is_diverging = FALSE, m)
    }
  } else if (is_bottom && ci <= 2) {
    par(mar = c(0, 0, 0, 0))
    plot.new()
  } else if (!is_bottom) {
    r <- rasters[[m]][[ci]]
    plot(r, col = PAL_CUR, range = ranges_raw[[m]], legend = FALSE,
         axes = FALSE, mar = c(0.2, 0, 0.2, 0))
    add_coastline()
    add_shapefile()
  } else {
    r_delta <- rasters[[m]][[ci]] - rasters[[m]][[DELTA_REF[ci]]]
    r_norm <- r_delta / max_abs[[m]]
    plot(r_norm, col = PAL_DIVERGING, range = c(-1, 1), legend = FALSE,
         axes = FALSE, mar = c(0.2, 0, 0.2, 0))
    add_coastline()
    add_shapefile()
  }
}

dev.off()
compile_tikz(paste0(BASE_DIR, "fig4-species-vs-ind-female.tex"))

# fig. S3.5 (male)
cat("\nSupplementary: male suitability maps\n")

rasters_M <- list()
ranges_raw_M <- list()

for (m in METRICS_M) {
  rasters_M[[m]] <- list()
  sub <- per_cell_M[per_cell_M$metric == m, ]
  cat(sprintf("  %s: %d cells\n", m, nrow(sub)))

  for (j in seq_along(COL_NAMES)) {
    df <- sub[, c("lon", "lat")]
    df$value <- sub[[COL_NAMES[j]]]
    rasters_M[[m]][[j]] <- rasterize_metric(df, template)
  }

  all_vals <- unlist(lapply(rasters_M[[m]], function(r) terra::values(r)))
  rng <- compute_limits(all_vals, robust = FALSE, type = "range")
  if (is.null(rng)) rng <- c(0, 1)
  ranges_raw_M[[m]] <- rng
}

max_abs_M <- list()
for (m in METRICS_M) {
  delta_vals <- c()
  for (j in 3:6) {
    delta_vals <- c(delta_vals,
      as.vector(rasters_M[[m]][[j]] - rasters_M[[m]][[DELTA_REF[j]]]))
  }
  rng <- compute_limits(delta_vals, robust = FALSE, type = "symmax")
  max_abs_M[[m]] <- if (is.null(rng)) 1 else max(abs(rng))
}

tikz(paste0(BASE_DIR, "figS2-5-species-vs-ind-male.tex"),
     width = 12, height = 12, standAlone = TRUE)

layout(lmat, widths = c(rep(1, 6), 0.4), heights = rep(1, 6))

par(oma = c(0, 0, 0, 0))
par(cex = 1)

for (k in 1:42) {
  ri <- (k - 1) %/% 7 + 1
  ci <- (k - 1) %% 7 + 1
  metric_idx <- (ri + 1) %/% 2
  is_bottom  <- ri %% 2 == 0
  m <- METRICS_M[metric_idx]

  if (ci == 7) {
    if (is_bottom) {
      draw_legend(PAL_DIVERGING, c(-1, 1), is_diverging = TRUE, m)
    } else {
      draw_legend(PAL_CUR, ranges_raw_M[[m]], is_diverging = FALSE, m)
    }
  } else if (is_bottom && ci <= 2) {
    par(mar = c(0, 0, 0, 0))
    plot.new()
  } else if (!is_bottom) {
    r <- rasters_M[[m]][[ci]]
    plot(r, col = PAL_CUR, range = ranges_raw_M[[m]], legend = FALSE,
         axes = FALSE, mar = c(0.2, 0, 0.2, 0))
    add_coastline()
    add_shapefile()
  } else {
    r_delta <- rasters_M[[m]][[ci]] - rasters_M[[m]][[DELTA_REF[ci]]]
    r_norm <- r_delta / max_abs_M[[m]]
    plot(r_norm, col = PAL_DIVERGING, range = c(-1, 1), legend = FALSE,
         axes = FALSE, mar = c(0.2, 0, 0.2, 0))
    add_coastline()
    add_shapefile()
  }
}

dev.off()
compile_tikz(paste0(BASE_DIR, "figS2-5-species-vs-ind-male.tex"))

# percentage difference maps

cat("\nPercentage difference maps\n")

METRICS_PCT_F <- c("LifeSpan", "MaxWgt", "Fecundity")
METRICS_PCT_M <- c("LifeSpan", "MaxWgt", "ReproEnergy")


PAL_PCT <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(21)

SCENARIOS <- c("cur", "w2", "w4")
SCOL_LABELS <- c("Current climate", "+2\\,°C", "+4\\,°C")

pct_diff_raster <- function(df_sub, sc) {
  p_sp <- df_sub[[paste0("P_species_", sc)]]
  p_ind <- df_sub[[paste0("P_ind_", sc)]]
  pct <- (p_sp - p_ind) / abs(p_ind) * 100
  pct[!is.finite(pct)] <- NA
  vals <- data.frame(lon = df_sub$lon, lat = df_sub$lat, value = pct)
  rasterize_metric(vals, template)
}

build_pct_rasters <- function(per_cell_sex, metrics) {
  rasters <- list()
  all_vals <- c()
  for (m in metrics) {
    sub <- per_cell_sex[per_cell_sex$metric == m, ]
    rasters[[m]] <- list()
    for (sc in SCENARIOS) {
      rasters[[m]][[sc]] <- pct_diff_raster(sub, sc)
      all_vals <- c(all_vals, terra::values(rasters[[m]][[sc]]))
    }
  }
  list(rasters = rasters, range = compute_limits(all_vals, robust = TRUE, type = "symmax"))
}

pct_F <- build_pct_rasters(per_cell_F, METRICS_PCT_F)
pct_M <- build_pct_rasters(per_cell_M, METRICS_PCT_M)

draw_pct_figure <- function(pct_obj, metrics, fig_name) {
  tikz(paste0(BASE_DIR, fig_name), width = 9, height = 6.75, standAlone = TRUE)

  n_rows <- length(metrics)
  n_cols <- 4
  layout(matrix(seq_len(n_rows * n_cols), nrow = n_rows, byrow = TRUE),
         widths = c(rep(1, 3), 0.35))

  par(oma = c(0, 0, 0, 0))
  par(cex = 1)

  rng <- pct_obj$range

  for (ri in seq_along(metrics)) {
    m <- metrics[ri]
    for (ci in 1:3) {
      sc <- SCENARIOS[ci]
      r <- pct_obj$rasters[[m]][[sc]]
      is_na <- all(is.na(terra::values(r)))

      if (is_na) {
        par(mar = c(0.1, 0.1, 0.1, 0.1))
        plot.new()
        plot.window(xlim = 0:1, ylim = 0:1)
        text(0.5, 0.5, "no data", cex = 0.8, col = "grey60")
      } else {
        plot(r, col = PAL_PCT, range = rng, legend = FALSE,
             axes = FALSE, mar = c(0.1, 0.1, 0.1, 0.1))
        add_coastline()
        add_shapefile()
      }
    }

    par(mar = c(0.5, 0.2, 0.5, 2.0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = rng)
    z <- seq(rng[1], rng[2], length.out = 100)
    ncol <- length(PAL_PCT)
    idx <- as.integer(cut(z, breaks = ncol, include.lowest = TRUE))
    cols <- PAL_PCT[idx]
    y <- seq(rng[1], rng[2], length.out = 100)
    rect(rep(0.15, 99), y[1:99], rep(0.40, 99), y[2:100],
         col = cols[1:99], border = NA)
    rect(0.15, rng[1], 0.40, rng[2], border = "grey60", lwd = 0.3, col = NA)
    lab_at <- pretty(rng, 5)
    lab_at <- lab_at[lab_at >= rng[1] & lab_at <= rng[2]]
    axis(4, at = lab_at, labels = paste0(lab_at, "\\%"), las = 1,
         cex.axis = 1.0, tcl = -0.2, mgp = c(0, 0.5, 0), pos = 0.40)
  }

  dev.off()
  compile_tikz(paste0(BASE_DIR, fig_name))
}

draw_pct_figure(pct_F, METRICS_PCT_F, "figS2-6-relative-bias-female.tex")
cat("  female pct diff done\n")

draw_pct_figure(pct_M, METRICS_PCT_M, "figS2-7-relative-bias-male.tex")
cat("  male pct diff done\n")

cat("\nDone.\n")
cat(sprintf("  %sfig4-species-vs-ind-female.tex (female)\n", BASE_DIR))
cat(sprintf("  %sfigS2-5-species-vs-ind-male.tex (male)\n", BASE_DIR))
cat(sprintf("  %sfigS2-6-relative-bias-female.tex (female pct diff)\n", BASE_DIR))
cat(sprintf("  %sfigS2-7-relative-bias-male.tex (male pct diff)\n", BASE_DIR))
