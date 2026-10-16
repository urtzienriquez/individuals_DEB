# load results for dynamic reporting in main text
cmp <- readRDS("../results/jensen_comparison.rds")
stats <- cmp$per_metric
cv_sorted <- readRDS("../results/param_cv.rds")

# sign distribution, current climate
per_cell_F <- cmp$per_cell[cmp$per_cell$sex == "F", ]
sign_tab_F <- table(per_cell_F$metric, per_cell_F$sign_cur)
stats$pct_pos <- as.vector(round(sign_tab_F[stats$metric, "overestimate"] / rowSums(sign_tab_F) * 100))
stats$pct_neg <- as.vector(round(sign_tab_F[stats$metric, "underestimate"] / rowSums(sign_tab_F) * 100))

# metric display names and units
METRIC_NAMES <- c(
  MaxWgt = "maximum body mass", MaxLen = "maximum body length",
  Tmax = "maximum body temperature", Tmin = "minimum body temperature",
  MinRes = "minimum reserve", Bsk = "basking time",
  Forage = "foraging time", Food = "food intake",
  Fecundity = "lifetime fecundity", Clutches = "number of clutches",
  MeanDailyDef = "daily water need", LifeSpan = "lifespan",
  R0 = "$R_0$", ReproEnergy = "reproductive energy"
)
METRIC_UNITS <- c(
  MaxWgt = "\\si{\\g}", MaxLen = "\\si{\\cm}",
  Tmax = "${}^{\\circ}$C", Tmin = "${}^{\\circ}$C",
  MinRes = "J", Bsk = "h", Forage = "h", Food = "\\si{\\g}",
  Fecundity = "eggs", Clutches = "",
  MeanDailyDef = "g\\,day$^{-1}$", LifeSpan = "years",
  R0 = "", ReproEnergy = "J"
)

# metrics in the figures
METRICS_SHOWN_F <- c("MaxWgt", "MaxLen", "Tmax", "Tmin", "MinRes",
                      "Bsk", "Forage", "Food", "Fecundity",
                      "MeanDailyDef", "LifeSpan", "R0")
METRICS_SHOWN_M <- c("MaxWgt", "MaxLen", "Tmax", "Tmin", "MinRes",
                      "Bsk", "Forage", "Food", "ReproEnergy",
                      "MeanDailyDef", "LifeSpan", "R0")

# mean J, current climate
signed_J <- aggregate(J_cur ~ metric, data = per_cell_F, FUN = mean)
names(signed_J)[2] <- "mean_J_cur"
stats <- merge(stats, signed_J, by = "metric", all.x = TRUE)
stats_F <- stats[stats$sex == "F", ]
stats_F <- stats_F[stats_F$metric %in% METRICS_SHOWN_F, ]

stats_F_absJ <- stats_F[order(-abs(stats_F$mean_J_cur)), ]

# amplified/attenuated at +4°C
amp_w4 <- stats_F[stats_F$direction_w4 == "amplified", ]
amp_w4 <- amp_w4[order(-amp_w4$rel_deltaJ_w4), ]
att_w4 <- stats_F[stats_F$direction_w4 == "attenuated", ]
att_w4 <- att_w4[order(att_w4$rel_deltaJ_w4), ]

# +2°C
amp_w2 <- stats_F[stats_F$direction_w2 == "amplified", ]
amp_w2 <- amp_w2[order(-amp_w2$rel_deltaJ_w2), ]
att_w2 <- stats_F[stats_F$direction_w2 == "attenuated", ]

# amplified at +2°C, attenuated at +4°C
switchers <- stats_F[stats_F$direction_w2 == "amplified" &
                      stats_F$direction_w4 == "attenuated", ]

top_over <- stats_F[order(-stats_F$mean_J_cur), ][1:3, ]
top_under <- stats_F[order(stats_F$mean_J_cur), ][1:3, ]

format_metric_list <- function(nms, vals) {
  n <- length(nms)
  if (n == 0) return("")
  items <- sapply(seq_len(n), function(i) {
    paste0(nms[i], " ($", round(vals[i], 0), "\\%$)")
  })
  if (n == 1) return(items)
  if (n == 2) return(paste(items, collapse = " and "))
  paste(paste(items[-n], collapse = ", "), "and", items[n])
}

att_text_w4 <- format_metric_list(
  METRIC_NAMES[att_w4$metric], att_w4$rel_deltaJ_w4)

amp_w2_min_name <- METRIC_NAMES[amp_w2$metric[which.min(amp_w2$rel_deltaJ_w2)]]
amp_w2_min_val  <- round(min(amp_w2$rel_deltaJ_w2), 0)
amp_w2_max_name <- METRIC_NAMES[amp_w2$metric[which.max(amp_w2$rel_deltaJ_w2)]]
amp_w2_max_val  <- round(max(amp_w2$rel_deltaJ_w2), 0)
att_text_w2 <- format_metric_list(
  METRIC_NAMES[att_w2$metric], att_w2$rel_deltaJ_w2)
switchers_text <- paste(METRIC_NAMES[switchers$metric], collapse = ", ")

top_sign <- stats_F[order(-stats_F$pct_sign_change_w4), ][1:3, ]

# males
per_cell_M <- cmp$per_cell[cmp$per_cell$sex == "M", ]
sign_tab_M <- table(per_cell_M$metric, per_cell_M$sign_cur)
stats$pct_pos_M <- as.vector(round(sign_tab_M[stats$metric, "overestimate"] / rowSums(sign_tab_M) * 100))
stats$pct_neg_M <- as.vector(round(sign_tab_M[stats$metric, "underestimate"] / rowSums(sign_tab_M) * 100))

signed_J_M <- aggregate(J_cur ~ metric, data = per_cell_M, FUN = mean)
names(signed_J_M)[2] <- "mean_J_cur_M"
stats <- merge(stats, signed_J_M, by = "metric", all.x = TRUE)
stats_M <- stats[stats$sex == "M", ]
stats_M <- stats_M[stats_M$metric %in% METRICS_SHOWN_M, ]
stats_M_absJ <- stats_M[order(-abs(stats_M$mean_J_cur_M)), ]

amp_w4_M <- stats_M[stats_M$direction_w4 == "amplified", ]
amp_w4_M <- amp_w4_M[order(-amp_w4_M$rel_deltaJ_w4), ]
att_w4_M <- stats_M[stats_M$direction_w4 == "attenuated", ]
att_w4_M <- att_w4_M[order(att_w4_M$rel_deltaJ_w4), ]

amp_w2_M <- stats_M[stats_M$direction_w2 == "amplified", ]
amp_w2_M <- amp_w2_M[order(-amp_w2_M$rel_deltaJ_w2), ]
att_w2_M <- stats_M[stats_M$direction_w2 == "attenuated", ]

switchers_M <- stats_M[stats_M$direction_w2 == "amplified" &
                        stats_M$direction_w4 == "attenuated", ]

top_over_M <- stats_M[order(-stats_M$mean_J_cur_M), ][1:3, ]
top_under_M <- stats_M[order(stats_M$mean_J_cur_M), ][1:3, ]

top_sign_M <- stats_M[order(-stats_M$pct_sign_change_w4), ][1:3, ]

att_text_w4_M <- format_metric_list(METRIC_NAMES[att_w4_M$metric], att_w4_M$rel_deltaJ_w4)
amp_w2_min_name_M <- METRIC_NAMES[amp_w2_M$metric[which.min(amp_w2_M$rel_deltaJ_w2)]]
amp_w2_min_val_M <- round(min(amp_w2_M$rel_deltaJ_w2), 0)
amp_w2_max_name_M <- METRIC_NAMES[amp_w2_M$metric[which.max(amp_w2_M$rel_deltaJ_w2)]]
amp_w2_max_val_M <- round(max(amp_w2_M$rel_deltaJ_w2), 0)
att_text_w2_M <- format_metric_list(METRIC_NAMES[att_w2_M$metric], att_w2_M$rel_deltaJ_w2)
switchers_text_M <- paste(METRIC_NAMES[switchers_M$metric], collapse = ", ")

# relative bias (%)

per_cell_F$rel_bias <- (per_cell_F$P_species_cur - per_cell_F$P_ind_cur) /
                        abs(per_cell_F$P_ind_cur) * 100
per_cell_F$rel_bias[!is.finite(per_cell_F$rel_bias)] <- NA

per_cell_M$rel_bias <- (per_cell_M$P_species_cur - per_cell_M$P_ind_cur) /
                        abs(per_cell_M$P_ind_cur) * 100
per_cell_M$rel_bias[!is.finite(per_cell_M$rel_bias)] <- NA

pct95_rel_bias_F <- aggregate(rel_bias ~ metric, data = per_cell_F,
                              FUN = function(x) round(quantile(abs(x), 0.95, na.rm = TRUE), 0))
names(pct95_rel_bias_F)[2] <- "pct95"

pct95_rel_bias_M <- aggregate(rel_bias ~ metric, data = per_cell_M,
                              FUN = function(x) round(quantile(abs(x), 0.95, na.rm = TRUE), 0))
names(pct95_rel_bias_M)[2] <- "pct95"

# fold changes
fold_tmin <- round(subset(amp_w4, metric == "Tmin")$rel_deltaJ_w4 / 100, 1)
fold_tmax <- round(subset(amp_w4, metric == "Tmax")$rel_deltaJ_w4 / 100, 1)
fold_r0   <- round(subset(amp_w4, metric == "R0")$rel_deltaJ_w4 / 100, 1)

att_pct_lifespan <- round(subset(att_w4, metric == "LifeSpan")$rel_deltaJ_w4, 0)
att_pct_minres   <- round(subset(att_w4, metric == "MinRes")$rel_deltaJ_w4, 0)

att_pct_lifespan_M <- round(subset(att_w4_M, metric == "LifeSpan")$rel_deltaJ_w4, 0)

fold_tmin_M <- round(subset(amp_w4_M, metric == "Tmin")$rel_deltaJ_w4, 0)
fold_tmax_M <- round(subset(amp_w4_M, metric == "Tmax")$rel_deltaJ_w4, 0)
fold_r0_M   <- round(subset(amp_w4_M, metric == "R0")$rel_deltaJ_w4, 0)

n_amp_w2_F <- sum(stats_F$direction_w2 == "amplified")
n_amp_w2_M <- sum(stats_M$direction_w2 == "amplified")
n_metrics_F <- nrow(stats_F)
n_metrics_M <- nrow(stats_M)

max_sign_change_w4 <- round(max(stats_F$pct_sign_change_w4), 0)

# normalised |J|
norm_absJ <- function(pc) {
  sapply(split(pc, pc$metric), function(x)
    round(mean(x$absJ_cur, na.rm = TRUE) / mean(abs(x$P_ind_cur), na.rm = TRUE) * 100, 1))
}
norm_absJ_F <- norm_absJ(per_cell_F)
norm_absJ_M <- norm_absJ(per_cell_M)
norm_absJ_F_top <- sort(norm_absJ_F[METRICS_SHOWN_F], decreasing = TRUE)
norm_absJ_M_top <- sort(norm_absJ_M[METRICS_SHOWN_M], decreasing = TRUE)

# CV by latitude band
eu_summary <- readRDS("../results/europe_summary.rds")
eu_summary$band <- cut(eu_summary$lat, c(-Inf, 40, 60, Inf),
                       labels = c("south", "central", "north"))
cv_bands <- do.call(rbind, lapply(union(METRICS_SHOWN_F, METRICS_SHOWN_M), function(m) {
  do.call(rbind, lapply(c("F", "M"), function(sx) {
    d <- eu_summary[eu_summary$sex == sx, ]
    v <- tapply(d[[paste0(m, "_cv")]], d$band, median, na.rm = TRUE)
    data.frame(metric = m, sex = sx, south = v[["south"]],
               central = v[["central"]], north = v[["north"]])
  }))
}))
cv_bands$ratio <- pmax(cv_bands$south, cv_bands$north) / cv_bands$central
cv_bands$margin <- ifelse(cv_bands$south >= cv_bands$north, "south", "north")
rownames(cv_bands) <- NULL
cv_ratio <- function(m, sx) round(cv_bands$ratio[cv_bands$metric == m & cv_bands$sex == sx], 1)
rm(eu_summary)

# consequences for forecasts
CONSEQ_METRICS <- list(F = c("LifeSpan", "MaxWgt", "Fecundity"),
                       M = c("LifeSpan", "MaxWgt", "ReproEnergy"))
conseq <- do.call(rbind, lapply(c("F", "M"), function(sx) {
  pc <- cmp$per_cell[cmp$per_cell$sex == sx, ]
  do.call(rbind, lapply(CONSEQ_METRICS[[sx]], function(m) {
    d <- pc[pc$metric == m, ]
    do.call(rbind, lapply(c("cur", "w2", "w4"), function(sc) {
      ps <- d[[paste0("P_species_", sc)]]
      pi <- d[[paste0("P_ind_", sc)]]
      rb <- (ps - pi) / abs(pi) * 100
      rb[!is.finite(rb)] <- NA
      agree <- if (sc == "cur") NA else
        mean(sign(ps - d$P_species_cur) == sign(pi - d$P_ind_cur), na.rm = TRUE) * 100
      data.frame(sex = sx, metric = m, scenario = sc,
                 pct_under = mean(ps < pi - 1e-10, na.rm = TRUE) * 100,
                 pct_over  = mean(ps > pi + 1e-10, na.rm = TRUE) * 100,
                 med_rel_bias = median(rb, na.rm = TRUE),
                 pct_agree = agree)
    }))
  }))
}))
cq <- function(sx, m, sc, col, digits = 0) {
  round(conseq[[col]][conseq$sex == sx & conseq$metric == m & conseq$scenario == sc], digits)
}
min_agree_w4 <- round(min(conseq$pct_agree[conseq$scenario == "w4"]), 0)
