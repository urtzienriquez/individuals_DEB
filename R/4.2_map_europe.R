#!/usr/bin/env Rscript
# base pdfs for fig. 2 and figs. S3.1-S3.3 (current climate)

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

IN_SUMMARY <- "../results/europe_summary.rds"
IN_JENSEN  <- "../results/europe_jensen.rds"
IN_SHAPE   <- "../raw_data/alpestris_shapefile/data_0.shp"
BASE_DIR   <- "../figures/1_base/"

RES  <- 1 / 6
CRS  <- "EPSG:4326"
XLIM <- c(-10, 40)
YLIM <- c(35, 70)

METRICS <- c("MaxWgt", "MaxLen", "MinRes",
             "Bsk", "Forage", "Food",
             "Tmax", "Tmin", "MeanDailyDef",
             "LifeSpan", "Fecundity", "R0")

METRICS_M <- c("MaxWgt", "MaxLen", "MinRes",
               "Bsk", "Forage", "Food",
               "Tmax", "Tmin", "MeanDailyDef",
               "LifeSpan", "ReproEnergy", "R0")

source("utils.R")

rasterize_metric <- function(summary, col, template) {
  rasterize(cbind(summary$lon, summary$lat), template,
            values = summary[[col]], fun = mean, na.rm = TRUE)
}

add_coastline <- function(col = "grey40", lwd = 0.3) {
  maps::map("worldHires", xlim = XLIM, ylim = YLIM,
            col = col, lwd = lwd, interior = FALSE, add = TRUE)
}

add_shapefile <- function(col = "grey20", lwd = 0.3) {
  plot(shp, col = NA, border = col, lwd = lwd, add = TRUE)
}

cat("Loading data\n")
summary <- readRDS(IN_SUMMARY)
jensen  <- readRDS(IN_JENSEN)
cat(sprintf("  Summary: %d cells\n", nrow(summary)))
cat(sprintf("  Jensen:  %d rows\n", nrow(jensen)))

# females: main text
summary <- summary[summary$sex == "F", ]
jensen  <- jensen[jensen$sex == "F", ]
cat(sprintf("  Female summary: %d cells\n", nrow(summary)))
cat(sprintf("  Female Jensen:  %d rows\n", nrow(jensen)))

template <- rast(xmin = XLIM[1], xmax = XLIM[2],
                 ymin = YLIM[1], ymax = YLIM[2],
                 resolution = RES, crs = CRS)

shp <- vect(IN_SHAPE)

# per-metric prediction maps (female)
cat("\nPrediction maps (female)\n")

for (i in seq_along(METRICS)) {
  m <- METRICS[i]
  sp_col <- paste0("sp_", m)
  mn_col <- paste0(m, "_mean")
  cv_col <- paste0(m, "_cv")

  r_sp <- rasterize_metric(summary, sp_col, template)
  r_mn <- rasterize_metric(summary, mn_col, template)

  is_zero <- all(terra::values(r_sp) == 0, na.rm = TRUE) &&
             all(terra::values(r_mn) == 0, na.rm = TRUE)

  pal <- hcl.colors(20, "Blues 2")

  if (!is_zero) {
    all_vals <- c(terra::values(r_sp), terra::values(r_mn))
    global_range <- compute_limits(all_vals, robust = TRUE, type = "range")
    if (is.null(global_range)) { is_zero <- TRUE }
  }

  tex_file <- paste0(BASE_DIR, "fig2-pred-", m, ".tex")
  tikz(tex_file, width = 10, height = 5, standAlone = TRUE)

  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.12))

  if (is_zero) {
    plot(r_sp, col = "grey95", legend = FALSE, axes = FALSE,
         mar = c(2, 2.5, 2, 4))
  } else {
    plot(r_sp, col = pal, range = global_range, legend = FALSE,
         axes = FALSE, mar = c(2, 2.5, 2, 4))
  }
  add_coastline()
  add_shapefile()

  if (is_zero) {
    plot(r_sp, col = "grey95", legend = FALSE, axes = FALSE,
         mar = c(2, 2.5, 2, 4))
  } else {
    plot(r_mn, col = pal, range = global_range, legend = FALSE,
         axes = FALSE, mar = c(2, 2.5, 2, 4))
  }
  add_coastline()
  add_shapefile()

  if (!is_zero) {
    par(mar = c(2, 0.3, 2, 2))
    global_breaks <- seq(global_range[1], global_range[2], length.out = 21)
    plot(NA, xlim = 0:1, ylim = global_range, axes = FALSE, xlab = "", ylab = "")
    rect(0, global_breaks[-21], 1, global_breaks[-1], col = pal, border = NA)
    ticks <- pretty(global_range, n = 5)
    ticks <- ticks[ticks >= global_range[1] & ticks <= global_range[2]]
    axis(4, at = ticks, mgp = c(3, 0.2, 0), cex.axis = 2.0, tcl = -0.2, las = 1)
  } else {
    par(mar = rep(0, 4))
    plot.new()
  }

  dev.off()
  compile_tikz(tex_file)

  cat(sprintf("  [%d/%d] %s\n", i, length(METRICS), m))
}

# fig. 2: Jensen bias maps (female)
cat("\nFig. 2: Jensen bias maps (female)\n")

pal_div  <- colorRampPalette(c(rgb(.10,.10,.44), "white", rgb(.80,.22,.18)))(21)

J_rasters <- lapply(METRICS, function(m) {
  sub <- jensen[jensen$metric == m, ]
  rasterize(cbind(sub$lon, sub$lat), template,
            values = sub$J, fun = mean, na.rm = TRUE)
})
names(J_rasters) <- METRICS

tex_file_fig3 <- paste0(BASE_DIR, "fig2-jensen-female.tex")
tikz(tex_file_fig3, width = 9, height = 6, standAlone = TRUE)
layout(matrix(1:12, nrow = 3, byrow = TRUE))

for (i in seq_along(METRICS)) {
  m <- METRICS[i]
  r <- J_rasters[[m]]
  is_zero <- all(terra::values(r) == 0, na.rm = TRUE)

  if (is_zero) {
    plot(r, col = "grey95", legend = FALSE, axes = FALSE,
         mar = c(1, 0, 1, 5.5))
  } else {
    rng <- compute_limits(terra::values(r), robust = TRUE, type = "symmax")
    if (is.null(rng)) {
      plot(r, col = "grey95", legend = FALSE, axes = FALSE,
           mar = c(1, 0, 1, 5.5))
    } else {
      max_abs <- max(abs(rng))
      plot(r, col = pal_div, range = c(-max_abs, max_abs), legend = TRUE,
           axes = FALSE, mar = c(1, 0, 1, 5.5),
           plg = list(cex = 1.8, title = ""))
    }
  }
  add_coastline()
  add_shapefile()
}

dev.off()
compile_tikz(tex_file_fig3)

# fig. S3.2: CV maps (female)
cat("\nFig. S3.2: CV maps (female)\n")

cv_cols <- paste0(METRICS, "_cv")
CV_rasters <- lapply(cv_cols, function(col) {
  rasterize_metric(summary, col, template)
})
names(CV_rasters) <- METRICS

pal_cv <- colorRampPalette(c("white", "red"))(20)

tex_file_fig4 <- paste0(BASE_DIR, "figS2-2-cv-female.tex")
tikz(tex_file_fig4, width = 9, height = 6, standAlone = TRUE)
layout(matrix(1:12, nrow = 3, byrow = TRUE))

for (i in seq_along(METRICS)) {
  m <- METRICS[i]
  r <- CV_rasters[[m]]
  is_zero <- all(terra::values(r) == 0, na.rm = TRUE)

  if (is_zero) {
    plot(r, col = "grey95", legend = FALSE, axes = FALSE,
         mar = c(1, 1, 1, 3))
  } else {
    pmax <- compute_limits(terra::values(r), robust = FALSE, type = "posmax")
    if (is.null(pmax)) {
      plot(r, col = "grey95", legend = FALSE, axes = FALSE,
           mar = c(1, 1, 1, 3))
    } else {
      plot(r, col = pal_cv, range = c(0, pmax), legend = TRUE,
           axes = FALSE, mar = c(1, 1, 1, 3),
           plg = list(cex = 1.8, title = ""))
    }
  }
  add_coastline()
  add_shapefile()
}

dev.off()
compile_tikz(tex_file_fig4)

cat("\nSupplementary: male figures\n")

summary_M <- readRDS(IN_SUMMARY)
jensen_M  <- readRDS(IN_JENSEN)
summary_M <- summary_M[summary_M$sex == "M", ]
jensen_M  <- jensen_M[jensen_M$sex == "M", ]
cat(sprintf("  Male summary: %d cells\n", nrow(summary_M)))
cat(sprintf("  Male Jensen:  %d rows\n", nrow(jensen_M)))

# per-metric prediction maps (male)
cat("\nPrediction maps (male)\n")

for (i in seq_along(METRICS_M)) {
  m <- METRICS_M[i]
  sp_col <- paste0("sp_", m)
  mn_col <- paste0(m, "_mean")

  r_sp <- rasterize_metric(summary_M, sp_col, template)
  r_mn <- rasterize_metric(summary_M, mn_col, template)

  is_zero <- all(terra::values(r_sp) == 0, na.rm = TRUE) &&
             all(terra::values(r_mn) == 0, na.rm = TRUE)

  pal <- hcl.colors(20, "Blues 2")

  if (!is_zero) {
    all_vals_M <- c(terra::values(r_sp), terra::values(r_mn))
    global_range_M <- compute_limits(all_vals_M, robust = TRUE, type = "range")
    if (is.null(global_range_M)) { is_zero <- TRUE }
  }

  tex_file <- paste0(BASE_DIR, "figS2-1-pred-", m, ".tex")
  tikz(tex_file, width = 10, height = 5, standAlone = TRUE)

  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.12))

  if (is_zero) {
    plot(r_sp, col = "grey95", legend = FALSE, axes = FALSE,
         mar = c(2, 2.5, 2, 4))
  } else {
    plot(r_sp, col = pal, range = global_range_M, legend = FALSE,
         axes = FALSE, mar = c(2, 2.5, 2, 4))
  }
  add_coastline()
  add_shapefile()

  if (is_zero) {
    plot(r_sp, col = "grey95", legend = FALSE, axes = FALSE,
         mar = c(2, 2.5, 2, 4))
  } else {
    plot(r_mn, col = pal, range = global_range_M, legend = FALSE,
         axes = FALSE, mar = c(2, 2.5, 2, 4))
  }
  add_coastline()
  add_shapefile()

  if (!is_zero) {
    par(mar = c(2, 0.3, 2, 2))
    global_breaks_M <- seq(global_range_M[1], global_range_M[2], length.out = 21)
    plot(NA, xlim = 0:1, ylim = global_range_M, axes = FALSE, xlab = "", ylab = "")
    rect(0, global_breaks_M[-21], 1, global_breaks_M[-1], col = pal, border = NA)
    ticks <- pretty(global_range_M, n = 5)
    ticks <- ticks[ticks >= global_range_M[1] & ticks <= global_range_M[2]]
    axis(4, at = ticks, mgp = c(3, 0.2, 0), cex.axis = 2.0, tcl = -0.2, las = 1)
  } else {
    par(mar = rep(0, 4))
    plot.new()
  }

  dev.off()
  compile_tikz(tex_file)
  cat(sprintf("  [%d/%d] %s\n", i, length(METRICS_M), m))
}

# fig. S3.1: Jensen bias maps (male)
cat("\nFig. S3.1: Jensen bias maps (male)\n")

J_rasters_M <- lapply(METRICS_M, function(m) {
  sub <- jensen_M[jensen_M$metric == m, ]
  rasterize(cbind(sub$lon, sub$lat), template,
            values = sub$J, fun = mean, na.rm = TRUE)
})
names(J_rasters_M) <- METRICS_M

tex_file_figS3 <- paste0(BASE_DIR, "figS2-1-jensen-male.tex")
tikz(tex_file_figS3, width = 9, height = 6, standAlone = TRUE)
layout(matrix(1:12, nrow = 3, byrow = TRUE))

for (i in seq_along(METRICS_M)) {
  m <- METRICS_M[i]
  r <- J_rasters_M[[m]]
  is_zero <- all(terra::values(r) == 0, na.rm = TRUE)

  if (is_zero) {
    plot(r, col = "grey95", legend = FALSE, axes = FALSE,
         mar = c(1, 0, 1, 5.5))
  } else {
    rng <- compute_limits(terra::values(r), robust = TRUE, type = "symmax")
    if (is.null(rng)) {
      plot(r, col = "grey95", legend = FALSE, axes = FALSE,
           mar = c(1, 0, 1, 5.5))
    } else {
      max_abs <- max(abs(rng))
      plot(r, col = pal_div, range = c(-max_abs, max_abs), legend = TRUE,
           axes = FALSE, mar = c(1, 0, 1, 5.5),
           plg = list(cex = 1.8, title = ""))
    }
  }
  add_coastline()
  add_shapefile()
}

dev.off()
compile_tikz(tex_file_figS3)

# fig. S3.3: CV maps (male)
cat("\nFig. S3.3: CV maps (male)\n")

cv_cols <- paste0(METRICS_M, "_cv")
CV_rasters_M <- lapply(cv_cols, function(col) {
  rasterize_metric(summary_M, col, template)
})
names(CV_rasters_M) <- METRICS_M

tex_file_figS4 <- paste0(BASE_DIR, "figS2-3-cv-male.tex")
tikz(tex_file_figS4, width = 9, height = 6, standAlone = TRUE)
layout(matrix(1:12, nrow = 3, byrow = TRUE))

for (i in seq_along(METRICS_M)) {
  m <- METRICS_M[i]
  r <- CV_rasters_M[[m]]
  is_zero <- all(terra::values(r) == 0, na.rm = TRUE)

  if (is_zero) {
    plot(r, col = "grey95", legend = FALSE, axes = FALSE,
         mar = c(1, 1, 1, 3))
  } else {
    pmax <- compute_limits(terra::values(r), robust = FALSE, type = "posmax")
    if (is.null(pmax)) {
      plot(r, col = "grey95", legend = FALSE, axes = FALSE,
           mar = c(1, 1, 1, 3))
    } else {
      plot(r, col = pal_cv, range = c(0, pmax), legend = TRUE,
           axes = FALSE, mar = c(1, 1, 1, 3),
           plg = list(cex = 1.8, title = ""))
    }
  }
  add_coastline()
  add_shapefile()
}

dev.off()
compile_tikz(tex_file_figS4)

cat("\nDone.\n")
