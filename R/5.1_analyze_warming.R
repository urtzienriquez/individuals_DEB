#!/usr/bin/env Rscript
# post-process europe-grid warming simulations (w2 + w4)

SCENARIOS <- c("w2", "w4")

metrics <- c("MaxWgt", "MaxLen", "Tmax", "Tmin", "MinRes",
             "Bsk", "Forage", "Food", "Fecundity", "Clutches",
             "MeanDailyDef", "LifeSpan", "R0", "ReproEnergy")

for (sc in SCENARIOS) {
  cat(sprintf("\n========== SCENARIO: %s ==========\n", sc))

  INPUT       <- sprintf("../results/europe_results_%s.rds", sc)
  OUT_CLEANED <- sprintf("../results/europe_cleaned_%s.rds", sc)
  OUT_SUMMARY <- sprintf("../results/europe_summary_%s.rds", sc)
  OUT_JENSEN  <- sprintf("../results/europe_jensen_%s.rds", sc)

  # load and clean
  cat("\nSection 1: Loading and cleaning results\n")
  cat(sprintf("  Input: %s\n", INPUT))

  results <- readRDS(INPUT)
  cat(sprintf("  Loaded: %d rows, %d columns\n", nrow(results), ncol(results)))

  results <- results[results$type != "error", ]
  cat(sprintf("  After dropping errors: %d rows\n", nrow(results)))
  cat(sprintf("    species(F): %d  species(M): %d  |  individual(F): %d  individual(M): %d\n",
      sum(results$type == "species" & results$sex == "F"),
      sum(results$type == "species" & results$sex == "M"),
      sum(results$type == "individual" & results$sex == "F"),
      sum(results$type == "individual" & results$sex == "M")))

  # survival flags from CauseDeath
  results$survived    <- results$CauseDeath == 0
  results$dehydrated  <- results$CauseDeath == 3
  results$temp_stress <- results$CauseDeath == 4

  # reproductive energy (J)
  cat("  Computing E_0 for all individuals ...\n")
  if (!exists("build_E0_lookup", inherits = FALSE)) source("utils.R")

  ind_ids_sc <- sort(unique(results$id[results$type == "individual"]))
  cat(sprintf("    %d unique individual IDs\n", length(ind_ids_sc)))

  if (!exists("sp_par_f", inherits = FALSE)) {
    sp_par_f <- build_par_list(
      mat_path         = "../sp/results_Ichthyosaura_alpestris.mat",
      res_html_path    = "../sp/Ichthyosaura_alpestris_res.html",
      report_html_path = "../sp/report_Ichthyosaura_alpestris.html",
      sex = "F", CT_min = -2.5, CT_max = 36
    )
  }

  e0_lookup <- build_E0_lookup(ind_ids_sc, sp_par_m = sp_par_f, base_dir = "../ind")
  cat(sprintf("    Computed E_0 for %d individuals\n", nrow(e0_lookup)))

  results$E_0 <- NA_real_
  is_ind <- results$type == "individual"
  match_idx <- match(results$id[is_ind], e0_lookup$id)
  results$E_0[is_ind] <- e0_lookup$E_0[match_idx]
  results$E_0[results$type == "species"] <- sp_par_f$E_0

  results$ReproEnergy <- results$Fecundity * results$E_0
  cat(sprintf("    ReproEnergy range: %.1f – %.1f J\n",
      min(results$ReproEnergy, na.rm = TRUE), max(results$ReproEnergy, na.rm = TRUE)))

  saveRDS(results, OUT_CLEANED)
  cat(sprintf("  Cleaned data -> %s\n", OUT_CLEANED))

  # per-cell summaries
  cat("\nSection 2: Computing per-cell summaries\n")

  ind <- results[results$type == "individual", ]
  sp  <- results[results$type == "species", ]

  # individual stats per cell and sex
  cat("  Aggregating individual statistics (by sex) ...\n")

  agg_mean <- aggregate(as.formula(paste0("cbind(", paste(metrics, collapse = ","),
      ") ~ lon + lat + sex")), data = ind, FUN = mean, na.rm = TRUE)
  names(agg_mean)[-(1:3)] <- paste0(metrics, "_mean")
  cat("    mean done\n")

  agg_median <- aggregate(as.formula(paste0("cbind(", paste(metrics, collapse = ","),
      ") ~ lon + lat + sex")), data = ind, FUN = median, na.rm = TRUE)
  names(agg_median)[-(1:3)] <- paste0(metrics, "_median")
  cat("    median done\n")

  agg_sd <- aggregate(as.formula(paste0("cbind(", paste(metrics, collapse = ","),
      ") ~ lon + lat + sex")), data = ind, FUN = sd, na.rm = TRUE)
  names(agg_sd)[-(1:3)] <- paste0(metrics, "_sd")
  cat("    sd done\n")

  agg_p5 <- aggregate(as.formula(paste0("cbind(", paste(metrics, collapse = ","),
      ") ~ lon + lat + sex")), data = ind, FUN = function(x)
          quantile(x, 0.05, na.rm = TRUE))
  names(agg_p5)[-(1:3)] <- paste0(metrics, "_p5")
  cat("    p5 done\n")

  agg_p95 <- aggregate(as.formula(paste0("cbind(", paste(metrics, collapse = ","),
      ") ~ lon + lat + sex")), data = ind, FUN = function(x)
          quantile(x, 0.95, na.rm = TRUE))
  names(agg_p95)[-(1:3)] <- paste0(metrics, "_p95")
  cat("    p95 done\n")

  agg_min <- aggregate(as.formula(paste0("cbind(", paste(metrics, collapse = ","),
      ") ~ lon + lat + sex")), data = ind, FUN = min, na.rm = TRUE)
  names(agg_min)[-(1:3)] <- paste0(metrics, "_min")
  cat("    min done\n")

  agg_max <- aggregate(as.formula(paste0("cbind(", paste(metrics, collapse = ","),
      ") ~ lon + lat + sex")), data = ind, FUN = max, na.rm = TRUE)
  names(agg_max)[-(1:3)] <- paste0(metrics, "_max")
  cat("    max done\n")

  cat("  Merging summary tables ...\n")
  summary <- agg_mean
  summary <- merge(summary, agg_median, by = c("lon", "lat", "sex"))
  summary <- merge(summary, agg_sd,     by = c("lon", "lat", "sex"))
  summary <- merge(summary, agg_p5,     by = c("lon", "lat", "sex"))
  summary <- merge(summary, agg_p95,    by = c("lon", "lat", "sex"))
  summary <- merge(summary, agg_min,    by = c("lon", "lat", "sex"))
  summary <- merge(summary, agg_max,    by = c("lon", "lat", "sex"))

  for (m in metrics) {
    m_mean <- summary[[paste0(m, "_mean")]]
    m_sd   <- summary[[paste0(m, "_sd")]]
    eps <- 1e-10
    summary[[paste0(m, "_cv")]] <- m_sd / (abs(m_mean) + eps)
  }
  cat("    cv computed\n")

  cat("  Aggregating survival fractions ...\n")
  binary_metrics <- c("survived", "dehydrated", "temp_stress")
  for (bm in binary_metrics) {
    fmla <- as.formula(paste0(bm, " ~ lon + lat + sex"))
    agg_bin <- aggregate(fmla, data = ind, FUN = mean, na.rm = TRUE)
    names(agg_bin)[4] <- paste0("frac_", bm)
    summary <- merge(summary, agg_bin, by = c("lon", "lat", "sex"))
  }
  cat("    fraction summaries done\n")

  # add population-average values
  cat("  Merging population-average values ...\n")
  sp_cols <- c("lon", "lat", "sex", metrics)
  sp_clean <- sp[, sp_cols]
  names(sp_clean)[-(1:3)] <- paste0("sp_", metrics)
  summary <- merge(summary, sp_clean, by = c("lon", "lat", "sex"), all.x = TRUE)
  cat(sprintf("    %d rows in summary (%d cells x 2 sexes)\n",
      nrow(summary), nrow(summary) / 2))

  saveRDS(summary, OUT_SUMMARY)
  cat(sprintf("  Summary -> %s\n", OUT_SUMMARY))

  # jensen bias
  cat("\nSection 3: Computing Jensen's inequality\n")

  J_list <- lapply(metrics, function(m) {
    sp_col <- paste0("sp_", m)
    ind_col <- paste0(m, "_mean")

    sp_val  <- summary[[sp_col]]
    ind_val <- summary[[ind_col]]

    J_raw <- sp_val - ind_val
    eps <- 1e-10
    J_rel <- J_raw / (abs(ind_val) + eps)

    sign_label <- ifelse(J_raw > eps, "overestimate",
                  ifelse(J_raw < -eps, "underestimate", "none"))

    data.frame(
      lon       = summary$lon,
      lat       = summary$lat,
      sex       = summary$sex,
      metric    = m,
      P_species = sp_val,
      P_ind     = ind_val,
      J         = J_raw,
      J_rel     = J_rel,
      sign      = sign_label,
      stringsAsFactors = FALSE
    )
  })
  jensen <- do.call(rbind, J_list)
  jensen$sign <- factor(jensen$sign,
    levels = c("underestimate", "none", "overestimate"))

  cat(sprintf("  %d rows (%d metrics x %d cells x 2 sexes)\n",
      nrow(jensen), length(metrics), nrow(summary) / 2))

  sign_tab <- table(jensen$metric, jensen$sign)
  cat("  Sign distribution (rows = metrics, cols = sign):\n")
  print(sign_tab)

  saveRDS(jensen, OUT_JENSEN)
  cat(sprintf("  Jensen -> %s\n", OUT_JENSEN))

  cat(sprintf("\nDone with %s.\n", sc))
}

cat("\nAll scenarios complete.\n")
