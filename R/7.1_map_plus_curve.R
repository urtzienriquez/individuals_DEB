#!/usr/bin/env Rscript
# figs. 3 and S3.4: delta |J| maps and temperature curves

library(terra)
library(dplyr)
library(tidyr)
library(mgcv)
library(maps)
library(mapdata)
library(tikzDevice)
library(purrr)

compile_tikz <- function(tex_path) {
  old_wd <- setwd(dirname(tex_path))
  on.exit(setwd(old_wd))
  tx <- readLines(basename(tex_path), warn = FALSE)
  tx <- append(tx, "\\let\\pgfimage\\includegraphics",
               after = max(grep("\\usepackage", tx, fixed = TRUE)))
  # lualatex for more memory
  engine <- if (nzchar(Sys.which("lualatex"))) "lualatex" else "pdflatex"
  writeLines(tx, basename(tex_path))
  system(paste(engine, "-interaction=nonstopmode", basename(tex_path)),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
}

RES  <- 1 / 6
CRS  <- "EPSG:4326"
XLIM <- c(-10, 40)
YLIM <- c(35, 70)

viridis_max <- "#FDE725"
viridis_min <- "#440154"
pal_div <- colorRampPalette(c(viridis_min, "white", viridis_max))(100)

# mean annual temperature
mat_file <- "../results/mat_europe.tif"
if (!file.exists(mat_file)) {
  folder <- "~/Documents/gcfolder"
  global_climate <- rast(paste0(folder, "/global_climate.nc"))
  eu_climate <- crop(global_climate, ext(XLIM, YLIM))
  tmean_monthly <- (eu_climate[[38:49]] + eu_climate[[50:61]]) / 20
  mat_raster <- app(tmean_monthly, mean, na.rm = TRUE)
  names(mat_raster) <- "mat"
  writeRaster(mat_raster, mat_file, overwrite = TRUE)
} else {
  mat_raster <- rast(mat_file)
}

comp <- readRDS("../results/jensen_comparison.rds")
d <- comp$per_cell
rm(comp)

# females: main text
d_F <- d[d$sex == "F", ]
# males: supplementary
d_M <- d[d$sex == "M", ]

cells <- distinct(d, lon, lat)
cells$mat <- as.vector(terra::extract(mat_raster, cells[, c("lon", "lat")])[, "mat"])
d_F <- left_join(d_F, cells, by = c("lon", "lat")) %>% filter(!is.na(mat))
d_M <- left_join(d_M, cells, by = c("lon", "lat")) %>% filter(!is.na(mat))

# rasterize
template <- rast(xmin = XLIM[1], xmax = XLIM[2],
                 ymin = YLIM[1], ymax = YLIM[2],
                 resolution = RES, crs = CRS)

metrics <- c("MaxWgt", "MaxLen", "MinRes",
             "Bsk", "Forage", "Food",
             "Tmax", "Tmin", "MeanDailyDef",
             "LifeSpan", "Fecundity", "R0")

metrics_M <- c("MaxWgt", "MaxLen", "MinRes",
               "Bsk", "Forage", "Food",
               "Tmax", "Tmin", "MeanDailyDef",
               "LifeSpan", "ReproEnergy", "R0")

source("utils.R")

cat("Rasterizing maps...\n")
map_rasters <- lapply(metrics, function(m) {
  sub <- filter(d_F, metric == m)
  rasterize(cbind(sub$lon, sub$lat), template,
            values = sub$delta_abs_J_w4, fun = mean, na.rm = TRUE)
})
names(map_rasters) <- metrics

# binned means
bin_width <- 1
d_F <- d_F %>% mutate(mat_bin = round(mat / bin_width) * bin_width)

binned <- d_F %>%
  group_by(metric, mat_bin) %>%
  summarise(n = n(),
            mean_delta = mean(delta_abs_J_w4, na.rm = TRUE),
            se_delta = sd(delta_abs_J_w4, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  filter(n >= 10)

# GAM fits
cat("Fitting GAMs...\n")
mat_grid <- seq(round(min(d_F$mat)), round(max(d_F$mat)), by = 0.1)

fit_gam_safe <- function(df, mg) {
  fit <- tryCatch(
    gam(delta_abs_J_w4 ~ s(mat, bs = "tp", k = 12),
        data = df, method = "REML", select = TRUE),
    error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble(mat = mg, pred = NA_real_, se = NA_real_))
  }
  p <- predict(fit, newdata = tibble(mat = mg), se.fit = TRUE)
  tibble(mat = mg, pred = as.vector(p$fit), se = as.vector(p$se.fit))
}

predictions <- d_F %>%
  nest(.by = metric) %>%
  mutate(
    pred_data = map(data, ~ fit_gam_safe(.x, mat_grid))
  ) %>%
  unnest(pred_data) %>%
  arrange(metric, mat) %>%
  select(metric, mat, pred, se)

shp <- NULL
shp_file <- "../raw_data/alpestris_shapefile/data_0.shp"
if (file.exists(shp_file)) {
  shp <- vect(shp_file)
  shp <- crop(shp, ext(XLIM[1], XLIM[2], YLIM[1], YLIM[2]))
  cat("  shapefile loaded\n")
}

add_coastline <- function(col = "grey40", lwd = 0.3) {
  maps::map("worldHires", xlim = XLIM, ylim = YLIM,
            col = col, lwd = lwd, interior = FALSE, add = TRUE)
}

add_shapefile <- function() {
  if (!is.null(shp)) {
    plot(shp, col = NA, border = "grey20", lwd = 0.3, add = TRUE)
  }
}

# figure
cat("Plotting...\n")
outfile_tex <- "../figures/1_base/fig3-warming-female.tex"
tikz(outfile_tex, width = 12, height = 7, standAlone = TRUE)

n_metrics <- length(metrics)
n_cells <- 12
gap <- 0.008
cw <- (1 - 3 * gap) / 4
ch <- (1 - 2 * gap) / 3

for (i in seq_len(n_cells)) {
  col_idx <- (i - 1) %% 4 + 1
  row_idx <- ceiling(i / 4)

  x_left  <- (col_idx - 1) * (cw + gap)
  x_right <- x_left + cw
  y_bottom <- (3 - row_idx) * (ch + gap)
  y_top    <- y_bottom + ch

  cw <- x_right - x_left
  ch <- y_top - y_bottom

  map_right  <- x_left + cw * 0.6
  curve_left <- x_left + cw * 0.6

  leg_height <- ch * 0.10
  new_flag   <- (i > 1)

  if (i <= n_metrics) {
    m <- metrics[i]
    r <- map_rasters[[m]]

    rng <- compute_limits(terra::values(r), robust = TRUE, type = "symmax")
    if (is.null(rng)) rng <- c(-1, 1)
    per_lim <- max(abs(rng))
    if (per_lim == 0) per_lim <- 1

    cd <- filter(predictions, metric == m)
    bd <- filter(binned, metric == m)
    max_abs <- max(abs(cd$pred), abs(bd$mean_delta), na.rm = TRUE) * 1.05
    if (is.na(max_abs) || max_abs == 0) max_abs <- 1

    # map
    vals <- terra::values(r)[, 1]
    breaks <- seq(-per_lim, per_lim, length.out = 101)
    cuts <- findInterval(vals, breaks, all.inside = TRUE)
    cols <- pal_div[cuts]
    cols[is.na(cols)] <- "white"
    r_mat <- matrix(cols, nrow = nrow(r), ncol = ncol(r), byrow = TRUE)

    par(fig = c(x_left, map_right, y_bottom, y_top),
        new = new_flag, mar = c(2.2, 0.3, 0.5, 0))
    plot.new()
    plot.window(xlim = XLIM, ylim = YLIM, xaxs = "i", yaxs = "i")
    rasterImage(r_mat, XLIM[1], YLIM[1], XLIM[2], YLIM[2],
                interpolate = FALSE)
    add_coastline()
    add_shapefile()

    z_vals <- seq(-1, 1, length.out = 100)
    idx <- findInterval(z_vals, seq(-1, 1, length.out = 101),
                        all.inside = TRUE)
    grad_img <- as.raster(matrix(pal_div[idx], nrow = 1))
    par(fig = c(x_left, map_right, y_bottom + leg_height * 0.4, y_bottom + leg_height * 1.4),
        new = TRUE, mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    rasterImage(grad_img, 0.12, 0.4, 0.88, 0.8, interpolate = TRUE)
    text(0.12, 0.08, round(-per_lim, 2), cex = 0.8, adj = c(0.5, 0),
         col = "grey40")
    text(0.88, 0.08, round(per_lim, 2),  cex = 0.8, adj = c(0.5, 0),
         col = "grey40")

    # curve
    par(fig = c(curve_left, x_right, y_bottom, y_top),
        new = TRUE, mar = c(2.2, 0, 0.5, 1.8))

    plot(cd$pred, cd$mat, type = "n",
         xlim = c(-max_abs, max_abs),
         ylim = rev(range(cd$mat)),
         axes = FALSE, xlab = "", ylab = "")

    abline(v = 0, lty = 2, col = "grey40", lwd = 0.3)

    pos <- which(cd$pred > 0)
    if (length(pos) > 0) {
      runs <- split(pos, cumsum(c(1, diff(pos) != 1)))
      for (run in runs) {
        pp <- cd[run, ]
        polygon(
          c(0, pp$pred, 0),
          c(pp$mat[1], pp$mat, pp$mat[length(pp$mat)]),
          col = adjustcolor(viridis_max, alpha.f = 0.35), border = NA
        )
      }
    }

    neg <- which(cd$pred < 0)
    if (length(neg) > 0) {
      runs <- split(neg, cumsum(c(1, diff(neg) != 1)))
      for (run in runs) {
        nn <- cd[run, ]
        polygon(
          c(nn$pred, 0, 0),
          c(nn$mat, nn$mat[length(nn$mat)], nn$mat[1]),
          col = adjustcolor(viridis_min, alpha.f = 0.35), border = NA
        )
      }
    }

    points(bd$mean_delta, bd$mat_bin, pch = 16, cex = 0.25,
           col = "grey20")

    lines(cd$pred, cd$mat, lwd = 0.4, col = "black")

    axis(1, cex.axis = 0.8, tcl = -0.2, mgp = c(0, 0.12, 0))
    axis(3, labels = FALSE, tcl = 0.2)
    axis(4, cex.axis = 0.8, tcl = -0.2, mgp = c(0, 0.3, 0), las = 1)

  } else {
    par(fig = c(x_left, x_right, y_bottom, y_top), new = new_flag,
        mar = rep(0, 4))
    plot.new()
  }
}

dev.off()
compile_tikz(outfile_tex)
cat("Saved to", outfile_tex, "\n")

# males: supplementary
cat("\nSupplementary: male map+curve\n")

map_rasters_M <- lapply(metrics_M, function(m) {
  sub <- filter(d_M, metric == m)
  rasterize(cbind(sub$lon, sub$lat), template,
            values = sub$delta_abs_J_w4, fun = mean, na.rm = TRUE)
})
names(map_rasters_M) <- metrics_M

d_M <- d_M %>% mutate(mat_bin = round(mat / bin_width) * bin_width)

binned_M <- d_M %>%
  group_by(metric, mat_bin) %>%
  summarise(n = n(),
            mean_delta = mean(delta_abs_J_w4, na.rm = TRUE),
            se_delta = sd(delta_abs_J_w4, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  filter(n >= 10)

mat_grid_M <- seq(round(min(d_M$mat)), round(max(d_M$mat)), by = 0.1)

predictions_M <- d_M %>%
  nest(.by = metric) %>%
  mutate(
    pred_data = map(data, ~ fit_gam_safe(.x, mat_grid_M))
  ) %>%
  unnest(pred_data) %>%
  arrange(metric, mat) %>%
  select(metric, mat, pred, se)

outfile_tex_M <- "../figures/1_base/figS2-4-warming-male.tex"
tikz(outfile_tex_M, width = 12, height = 7, standAlone = TRUE)

n_metrics_M <- length(metrics_M)

for (i in seq_len(n_cells)) {
  col_idx <- (i - 1) %% 4 + 1
  row_idx <- ceiling(i / 4)

  x_left  <- (col_idx - 1) * (cw + gap)
  x_right <- x_left + cw
  y_bottom <- (3 - row_idx) * (ch + gap)
  y_top    <- y_bottom + ch

  cw_local <- x_right - x_left
  ch_local <- y_top - y_bottom

  map_right  <- x_left + cw_local * 0.6
  curve_left <- x_left + cw_local * 0.6

  leg_height <- ch_local * 0.10
  new_flag   <- (i > 1)

  if (i <= n_metrics_M) {
    m <- metrics_M[i]
    r <- map_rasters_M[[m]]

    rng <- compute_limits(terra::values(r), robust = TRUE, type = "symmax")
    if (is.null(rng)) rng <- c(-1, 1)
    per_lim <- max(abs(rng))
    if (per_lim == 0) per_lim <- 1

    cd <- filter(predictions_M, metric == m)
    bd <- filter(binned_M, metric == m)
    max_abs <- max(abs(cd$pred), abs(bd$mean_delta), na.rm = TRUE) * 1.05
    if (is.na(max_abs) || max_abs == 0) max_abs <- 1

    vals <- terra::values(r)[, 1]
    breaks <- seq(-per_lim, per_lim, length.out = 101)
    cuts <- findInterval(vals, breaks, all.inside = TRUE)
    cols <- pal_div[cuts]
    cols[is.na(cols)] <- "white"
    r_mat <- matrix(cols, nrow = nrow(r), ncol = ncol(r), byrow = TRUE)

    par(fig = c(x_left, map_right, y_bottom, y_top),
        new = new_flag, mar = c(2.2, 0.3, 0.5, 0))
    plot.new()
    plot.window(xlim = XLIM, ylim = YLIM, xaxs = "i", yaxs = "i")
    rasterImage(r_mat, XLIM[1], YLIM[1], XLIM[2], YLIM[2],
                interpolate = FALSE)
    add_coastline()
    add_shapefile()

    z_vals <- seq(-1, 1, length.out = 100)
    idx <- findInterval(z_vals, seq(-1, 1, length.out = 101),
                        all.inside = TRUE)
    grad_img <- as.raster(matrix(pal_div[idx], nrow = 1))
    par(fig = c(x_left, map_right, y_bottom + leg_height * 0.4, y_bottom + leg_height * 1.4),
        new = TRUE, mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    rasterImage(grad_img, 0.12, 0.4, 0.88, 0.8, interpolate = TRUE)
    text(0.12, 0.08, round(-per_lim, 2), cex = 0.8, adj = c(0.5, 0),
         col = "grey40")
    text(0.88, 0.08, round(per_lim, 2),  cex = 0.8, adj = c(0.5, 0),
         col = "grey40")

    par(fig = c(curve_left, x_right, y_bottom, y_top),
        new = TRUE, mar = c(2.2, 0, 0.5, 1.8))

    plot(cd$pred, cd$mat, type = "n",
         xlim = c(-max_abs, max_abs),
         ylim = rev(range(cd$mat)),
         axes = FALSE, xlab = "", ylab = "")

    abline(v = 0, lty = 2, col = "grey40", lwd = 0.3)

    pos <- which(cd$pred > 0)
    if (length(pos) > 0) {
      runs <- split(pos, cumsum(c(1, diff(pos) != 1)))
      for (run in runs) {
        pp <- cd[run, ]
        polygon(
          c(0, pp$pred, 0),
          c(pp$mat[1], pp$mat, pp$mat[length(pp$mat)]),
          col = adjustcolor(viridis_max, alpha.f = 0.35), border = NA
        )
      }
    }

    neg <- which(cd$pred < 0)
    if (length(neg) > 0) {
      runs <- split(neg, cumsum(c(1, diff(neg) != 1)))
      for (run in runs) {
        nn <- cd[run, ]
        polygon(
          c(nn$pred, 0, 0),
          c(nn$mat, nn$mat[length(nn$mat)], nn$mat[1]),
          col = adjustcolor(viridis_min, alpha.f = 0.35), border = NA
        )
      }
    }

    points(bd$mean_delta, bd$mat_bin, pch = 16, cex = 0.25,
           col = "grey20")

    lines(cd$pred, cd$mat, lwd = 0.4, col = "black")

    axis(1, cex.axis = 0.8, tcl = -0.2, mgp = c(0, 0.12, 0))
    axis(3, labels = FALSE, tcl = 0.2)
    axis(4, cex.axis = 0.8, tcl = -0.2, mgp = c(0, 0.3, 0), las = 1)

  } else {
    par(fig = c(x_left, x_right, y_bottom, y_top), new = new_flag,
        mar = rep(0, 4))
    plot.new()
  }
}

dev.off()
compile_tikz(outfile_tex_M)
cat("Saved to", outfile_tex_M, "\n")
