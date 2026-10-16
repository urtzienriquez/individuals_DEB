# ichthyosaura alpestris individual 8045: biophysical + DEB

library(NicheMapR)
library(R.matlab)
source("utils.R")
source("params_definition.R")
source("compute_r0_correct.R")

ind_id <- 8045
ind <- Filter(function(x) x$id == ind_id, ind_pars)[[1]]
if (is.null(ind)) stop(sprintf("Individual %d not found in ind_pars", ind_id))

# microclimate
loc <- c(15.0, 49.5) # czech republic
# loc <- c(0.0, 42)     # Spain

micro_in <- micro_global(
  loc = loc,
  timeinterval = 365,
  nyears = 3,
  runmoist = 1,
  runshade = 1,
  snowmodel = 1,
  soiltype = 4,
  warm = 0
)
micro <- extend_micro(micro_in)

deb_par <- ind[!names(ind) %in% c("id", "sex")]
ecto <- do.call(ectotherm, c(bio_par, deb_par))

environ <- as.data.frame(ecto$environ)
debout <- as.data.frame(ecto$debout)
yearout <- as.data.frame(ecto$yearout)
yearsout <- as.data.frame(ecto$yearsout)

cat("\n========================================\n")
cat(sprintf("  Individual %d  -  Ichthyosaura alpestris\n", ind_id))
cat("  Start: post-metamorphic juvenile (stage 1)\n")
cat("========================================\n\n")

cat("Key parameters\n")
cat(sprintf("  p_Xm = %.1f J/d.cm^2    (with NicheMapR *10 factor)\n",
  ind$p_Xm))
cat(sprintf("  L_b  = %.4f cm\n", ind$L_b))
cat(sprintf("  E_0  = %.2f J\n", ind$E_0))
cat(sprintf("  V_init = %.4f cm^3\n", ind$V_init))
cat(sprintf("  z = %.3f, kap = %.4f, p_M = %.1f J/d.cm^3\n",
  ind$z, ind$kap, ind$p_M * 24))
cat("\n")

cat("Yearout (whole-life summary)\n")
print(yearout)
cat("\nYearsout (annual summary)\n")
print(yearsout)

cat("\nKey output metrics\n")
cat(sprintf("  Clutches:    %s\n", yearout$CLUTCHES))
cat(sprintf("  Fecundity:   %s\n", yearout$FECUNDITY))
cat(sprintf("  R0:          %s\n", yearout$R0))
cat(sprintf("  Pct_Des:     %s\n", yearout$Pct_Des))
cat(sprintf("  Total drink: %s\n", sum(yearsout$Drink)))
cat(sprintf("  Total food:  %s\n", sum(yearsout$Food)))
cat(sprintf("  Max weight:  %s\n", max(yearsout$MaxWgt, na.rm = TRUE)))

# corrected R0
res <- compute_r0_correct(debout, yearsout, yearout)
cat(sprintf("  Corrected R0:  %.6f  (built-in: %.6f, ratio %.2f)\n",
    res$R0_correct, res$R0_builtin, res$R0_correct / res$R0_builtin))
cat(sprintf("  Clutches detected from E_B drops: %d\n", res$n_clutches))
cat("\n")

# daily water needs
masbal <- as.data.frame(ecto$masbal)
alive_start <- min(which(masbal$YEAR == 1 & masbal$DAY >= 227))
alive_end   <- max(which(masbal$H2OCut_g > 0))
alive_mb <- masbal[alive_start:alive_end, ]
alive_mb$date <- paste(alive_mb$YEAR, alive_mb$DAY, sep = "-")
daily_deficit <- aggregate(
  H2OBal_g ~ date,
  data = alive_mb,
  FUN = function(x) sum(x[x < 0])
)
cat("Daily drinking water needs\n")
cat(sprintf("  Days alive:            %d\n", nrow(daily_deficit)))
cat(sprintf("  Days with deficit:     %d\n", sum(daily_deficit$H2OBal_g < 0)))
cat(sprintf("  Total drinking need:   %.2f g\n", abs(sum(daily_deficit$H2OBal_g))))
cat(sprintf("  Mean daily need:       %.4f g/day\n",
  abs(mean(daily_deficit$H2OBal_g))))
cat(sprintf("  First day need:        %.4f g\n",
  abs(daily_deficit$H2OBal_g[1])))
cat(sprintf("  Last day need:         %.4f g\n",
  abs(daily_deficit$H2OBal_g[nrow(daily_deficit)])))
cat("\n")

# growth plot
n_hours <- length(debout$WETMASS)
x_time <- seq(1, n_hours) / 24
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

plot(x_time, debout$WETMASS,
  type = "l", xlab = "date (days)",
  ylab = "wet mass (g)", col = "pink", lwd = 2,
  ylim = c(0, max(debout$WETMASS))
)
points(x_time, debout$V, type = "l", col = "dark green", lwd = 2)
points(x_time, debout$WETMASS - debout$WETGONAD, type = "l", lwd = 2,
  col = "brown")
points(x_time, debout$WETMASS - debout$WETGONAD - debout$WETGUT,
  type = "l", lwd = 2, col = "grey"
)
h_hp <- which(debout$E_H > deb_par$E_Hp)[1]
if (!is.na(h_hp)) abline(v = x_time[h_hp], lty = 2, col = "grey")
legend("topleft",
  c("repro. buffer", "food in gut", "reserve", "structure"),
  lty = 1, col = c("pink", "brown", "grey", "dark green"), bty = "n"
)
mtext("(a) DEB biomass components", side = 3, adj = 0, line = 0.5)

plot(x_time, environ$TC,
  type = "l", xlab = "date (days)",
  ylab = expression("T"["b"] * " (" * degree * "C)"), col = "blue", lwd = 1
)
abline(h = c(bio_par$T_F_min, bio_par$T_F_max), lty = 3, col = "orange")
abline(h = bio_par$CT_max, lty = 2, col = "red")
abline(h = bio_par$CT_min, lty = 2, col = "blue")
legend("topleft",
  c("Tb", "activity range", "CT limits"),
  lty = c(1, 3, 2), col = c("blue", "orange", "red"), bty = "n"
)
mtext("(b) Body temperature", side = 3, adj = 0, line = 0.5)
