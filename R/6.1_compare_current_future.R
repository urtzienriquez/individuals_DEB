#!/usr/bin/env Rscript
# jensen comparison across climate scenarios

METRICS <- c("MaxWgt", "MaxLen", "Tmax", "Tmin", "MinRes",
             "Bsk", "Forage", "Food", "Fecundity", "Clutches",
             "MeanDailyDef", "LifeSpan", "R0", "ReproEnergy")

cat("Loading Jensen results\n")

j_cur <- readRDS("../results/europe_jensen.rds")
j_w2  <- readRDS("../results/europe_jensen_w2.rds")
j_w4  <- readRDS("../results/europe_jensen_w4.rds")

cat(sprintf("  Current: %d rows\n", nrow(j_cur)))
cat(sprintf("  w2:      %d rows\n", nrow(j_w2)))
cat(sprintf("  w4:      %d rows\n", nrow(j_w4)))

cat("\nAligning scenarios by cell\n")

comparison <- merge(
  j_cur[, c("lon", "lat", "sex", "metric", "J", "J_rel", "sign", "P_species", "P_ind")],
  j_w2[,  c("lon", "lat", "sex", "metric", "J", "J_rel", "sign", "P_species", "P_ind")],
  by = c("lon", "lat", "sex", "metric"),
  suffixes = c("_cur", "_w2")
)
comparison <- merge(
  comparison,
  j_w4[, c("lon", "lat", "sex", "metric", "J", "J_rel", "sign", "P_species", "P_ind")],
  by = c("lon", "lat", "sex", "metric")
)
names(comparison)[grep("^J$", names(comparison))] <- "J_w4"
names(comparison)[grep("^J_rel$", names(comparison))] <- "J_rel_w4"
names(comparison)[grep("^sign$", names(comparison))] <- "sign_w4"
names(comparison)[grep("^P_species$", names(comparison))] <- "P_species_w4"
names(comparison)[grep("^P_ind$", names(comparison))] <- "P_ind_w4"

cat(sprintf("  Aligned: %d rows (%d cells × %d metrics × 2 sexes)\n",
    nrow(comparison), nrow(comparison) / length(METRICS) / 2, length(METRICS)))

# delta |J| = |J_future| - |J_current|
cat("\nComputing Δ|J| values\n")

comparison$absJ_cur   <- abs(comparison$J_cur)
comparison$absJ_w2    <- abs(comparison$J_w2)
comparison$absJ_w4    <- abs(comparison$J_w4)

comparison$delta_abs_J_w2  <- comparison$absJ_w2 - comparison$absJ_cur
comparison$delta_abs_J_w4  <- comparison$absJ_w4 - comparison$absJ_cur

comparison$raw_delta_J_w2  <- comparison$J_w2 - comparison$J_cur
comparison$raw_delta_J_w4  <- comparison$J_w4 - comparison$J_cur

comparison$deltaJ_w2 <- comparison$delta_abs_J_w2
comparison$deltaJ_w4 <- comparison$delta_abs_J_w4

# sign changes
comparison$sign_change_w2 <- comparison$sign_cur != comparison$sign_w2
comparison$sign_change_w4 <- comparison$sign_cur != comparison$sign_w4

comparison$direction_w2 <- ifelse(comparison$delta_abs_J_w2 > 1e-10, "increase",
                           ifelse(comparison$delta_abs_J_w2 < -1e-10, "decrease", "unchanged"))
comparison$direction_w4 <- ifelse(comparison$delta_abs_J_w4 > 1e-10, "increase",
                           ifelse(comparison$delta_abs_J_w4 < -1e-10, "decrease", "unchanged"))

# per-metric stats by sex
cat("\nPer-metric statistics (by sex)\n")

stat_rows <- lapply(METRICS, function(m) {
  lapply(c("F", "M"), function(sx) {
    sub <- comparison[comparison$metric == m & comparison$sex == sx, ]

    data.frame(
      sex              = sx,
      metric           = m,
      n_cells          = nrow(sub),
      mean_absJ_cur       = mean(sub$absJ_cur, na.rm = TRUE),
      mean_absJ_w2        = mean(sub$absJ_w2, na.rm = TRUE),
      mean_absJ_w4        = mean(sub$absJ_w4, na.rm = TRUE),
      mean_delta_abs_J_w2 = mean(sub$delta_abs_J_w2, na.rm = TRUE),
      mean_delta_abs_J_w4 = mean(sub$delta_abs_J_w4, na.rm = TRUE),
      median_delta_abs_J_w2 = median(sub$delta_abs_J_w2, na.rm = TRUE),
      median_delta_abs_J_w4 = median(sub$delta_abs_J_w4, na.rm = TRUE),
      pct_increase_w2  = mean(sub$delta_abs_J_w2 > 1e-10, na.rm = TRUE) * 100,
      pct_increase_w4  = mean(sub$delta_abs_J_w4 > 1e-10, na.rm = TRUE) * 100,
      pct_decrease_w2  = mean(sub$delta_abs_J_w2 < -1e-10, na.rm = TRUE) * 100,
      pct_decrease_w4  = mean(sub$delta_abs_J_w4 < -1e-10, na.rm = TRUE) * 100,
      pct_sign_change_w2 = mean(sub$sign_change_w2, na.rm = TRUE) * 100,
      pct_sign_change_w4 = mean(sub$sign_change_w4, na.rm = TRUE) * 100,
      direction_w2     = ifelse(mean(sub$delta_abs_J_w2, na.rm = TRUE) > 1e-10, "amplified",
                          ifelse(mean(sub$delta_abs_J_w2, na.rm = TRUE) < -1e-10, "attenuated", "unchanged")),
      direction_w4     = ifelse(mean(sub$delta_abs_J_w4, na.rm = TRUE) > 1e-10, "amplified",
                          ifelse(mean(sub$delta_abs_J_w4, na.rm = TRUE) < -1e-10, "attenuated", "unchanged")),
      stringsAsFactors = FALSE
    )
  })
})
stats <- do.call(rbind, unlist(stat_rows, recursive = FALSE))

cat("\n  Summary table (Δ|J| = |J_warmed| − |J_current|):\n")
cat(sprintf("  %-4s %-12s %8s %8s %8s %8s\n",
    "Sex", "Metric", "Δ|J|_w2", "Δ|J|_w4", "%Inc_w2", "%Inc_w4"))
for (i in seq_len(nrow(stats))) {
  cat(sprintf("  %-4s %-12s %8.4f %8.4f %7.1f%% %7.1f%%\n",
      stats$sex[i], stats$metric[i],
      stats$mean_delta_abs_J_w2[i], stats$mean_delta_abs_J_w4[i],
      stats$pct_increase_w2[i], stats$pct_increase_w4[i]))
}

# relative delta |J| (%)
stats$rel_deltaJ_w2 <- stats$mean_delta_abs_J_w2 / (stats$mean_absJ_cur + 1e-10) * 100
stats$rel_deltaJ_w4 <- stats$mean_delta_abs_J_w4 / (stats$mean_absJ_cur + 1e-10) * 100

OUT <- "../results/jensen_comparison.rds"
saveRDS(list(per_cell = comparison, per_metric = stats), OUT)
cat(sprintf("\nSaved to %s\n", OUT))
