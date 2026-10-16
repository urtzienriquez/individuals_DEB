#!/usr/bin/env Rscript
# per-individual DEB growth trajectory plots

library(NicheMapR)

source("utils.R")
source("params_definition.R")

run_mode <- "all" # "sp" or "all"
fig_dir <- "../figures/1_base/"

all_pars <- c(list(sp_par_f, sp_par_m), ind_pars)
all_pars[[1]]$id    <- "sp"
all_pars[[1]]$label <- "Population-average (female)"
all_pars[[2]]$id    <- "sp_male"
all_pars[[2]]$label <- "Population-average (male)"

for (i in seq_along(ind_pars)) {
  idx <- i + 2
  all_pars[[idx]]$label <- sprintf("PID %s (%s)", all_pars[[idx]]$id,
                                    all_pars[[idx]]$sex)
}

if (run_mode == "sp") {
  all_pars <- all_pars[1:2]
}

# microclimate for two sites
cat("Running micro_global for cold (Czechia) and warm (Spain) sites ... ")
flush.console()

micro_cold_in <- micro_global(
  loc = c(15.0, 49.5), timeinterval = 365, nyears = 3,
  runmoist = 1, snowmodel = 1, runshade = 1, soiltype = 4, warm = 0
)
micro_cold <- extend_micro(micro_cold_in)

micro_warm_in <- micro_global(
  loc = c(-5.127, 43.253), timeinterval = 365, nyears = 3,
  runmoist = 1, snowmodel = 1, runshade = 1, soiltype = 4, warm = 0
)
micro_warm <- extend_micro(micro_warm_in)

cat("done\n")
flush.console()

# simulations and plots
cat(sprintf("\nRunning %d simulations ...\n", length(all_pars)))

for (j in seq_along(all_pars)) {
  p <- all_pars[[j]]
  cat(sprintf("  [%d/%d] %s ... ", j, length(all_pars), p$label))
  flush.console()

  # cold site (czechia)
  assign("micro", micro_cold, .GlobalEnv)
  ecto_cold <- do.call(ectotherm, c(bio_par,
    p[!names(p) %in% c("id", "label", "sex")]))
  debout_cold <- as.data.frame(ecto_cold$debout)
  n_hours <- nrow(debout_cold)
  x_time <- seq_len(n_hours) / 24

  # warm site (spain)
  assign("micro", micro_warm, .GlobalEnv)
  ecto_warm <- do.call(ectotherm, c(bio_par,
    p[!names(p) %in% c("id", "label", "sex")]))
  debout_warm <- as.data.frame(ecto_warm$debout)

  h_hp_cold <- which(debout_cold$E_H > p$E_Hp)[1]
  h_hp_warm <- which(debout_warm$E_H > p$E_Hp)[1]
  puberty_text_cold <- if (!is.na(h_hp_cold)) {
    sprintf("puberty at day %.0f", x_time[h_hp_cold])
  } else {
    "NEVER reaches puberty"
  }
  puberty_text_warm <- if (!is.na(h_hp_warm)) {
    sprintf("puberty at day %.0f", x_time[h_hp_warm])
  } else {
    "NEVER reaches puberty"
  }

  total_fec_cold <- sum(ecto_cold$yearsout[, "Fec"], na.rm = TRUE)
  total_fec_warm <- sum(ecto_warm$yearsout[, "Fec"], na.rm = TRUE)

  is_male <- !is.null(p$sex) && p$sex == "M" || grepl("\\(male\\)", p$label, ignore.case = TRUE)
  if (is_male) {
    total_repro_energy_cold <- total_fec_cold * p$E_0
    total_repro_energy_warm <- total_fec_warm * p$E_0
  }

  pdf(file.path(fig_dir, sprintf("growth_%s.pdf", p$id)),
    width = 10, height = 7
  )

  par(mfrow = c(2, 1), mar = c(4, 4.5, 2, 1))

  # panel a: cold site
  ymax <- max(c(debout_cold$WETMASS, debout_warm$WETMASS), na.rm = TRUE)
  plot(x_time, debout_cold$WETMASS,
    type = "l",
    xlab = "Time (days)", ylab = "Mass (g)",
    col = "pink", lwd = 2, ylim = c(0, ymax)
  )
  lines(x_time, debout_cold$V, col = "dark green", lwd = 2)
  lines(x_time, debout_cold$WETMASS - debout_cold$WETGONAD, col = "brown",
    lwd = 2)
  lines(x_time, debout_cold$WETMASS - debout_cold$WETGONAD -
    debout_cold$WETGUT, col = "grey", lwd = 2)
  if (!is.na(h_hp_cold)) abline(v = x_time[h_hp_cold], lty = 2, col = "grey50")

  legend("topright",
    c("Total mass", "Structure", "Body - gonad", "Body - gonad - gut"),
    lty = 1, col = c("pink", "dark green", "brown", "grey"), bty = "n"
  )
  if (is_male) {
    mtext(
      sprintf("(a) Cold site (Czechia): %s  (%s, repro. energy = %.1f J)",
        p$label, puberty_text_cold, total_repro_energy_cold),
      side = 3, adj = 0, line = 0.5, cex = 0.85
    )
  } else {
    mtext(
      sprintf("(a) Cold site (Czechia): %s  (%s, fecundity = %.0f)",
        p$label, puberty_text_cold, total_fec_cold),
      side = 3, adj = 0, line = 0.5, cex = 0.85
    )
  }

  # panel b: warm site
  plot(x_time, debout_warm$WETMASS,
    type = "l",
    xlab = "Time (days)", ylab = "Mass (g)",
    col = "pink", lwd = 2, ylim = c(0, ymax)
  )
  lines(x_time, debout_warm$V, col = "dark green", lwd = 2)
  lines(x_time, debout_warm$WETMASS - debout_warm$WETGONAD, col = "brown",
    lwd = 2)
  lines(x_time, debout_warm$WETMASS - debout_warm$WETGONAD -
    debout_warm$WETGUT, col = "grey", lwd = 2)
  if (!is.na(h_hp_warm)) abline(v = x_time[h_hp_warm], lty = 2, col = "grey50")

  legend("topright",
    c("Total mass", "Structure", "Body - gonad", "Body - gonad - gut"),
    lty = 1, col = c("pink", "dark green", "brown", "grey"), bty = "n"
  )
  if (is_male) {
    mtext(
      sprintf("(b) Warm site (Spain): %s  (%s, repro. energy = %.1f J)",
        p$label, puberty_text_warm, total_repro_energy_warm),
      side = 3, adj = 0, line = 0.5, cex = 0.85
    )
  } else {
    mtext(
      sprintf("(b) Warm site (Spain): %s  (%s, fecundity = %.0f)",
        p$label, puberty_text_warm, total_fec_warm),
      side = 3, adj = 0, line = 0.5, cex = 0.85
    )
  }

  dev.off()
  cat("done\n")
  flush.console()
}

cat(sprintf("\nDone. %d figures saved to %s\n", length(all_pars), fig_dir))
