source("parse_ind_pars.R")
source("utils.R")

library(tikzDevice)
options(tikzLatexPackages = c(
  getOption("tikzLatexPackages"),
  "\\usepackage{underscore}"
))
options(tikzSanitizeCharacters = "%")
options(tikzReplacementCharacters = list("%" = "\\%"))
fig_folder <- "../figures/1_base/"

# population-average constants
sp_pars <- parse_species_pars("../sp/pars_init_Ichthyosaura_alpestris.m")
z_sp     <- sp_pars["z"]
z_m_sp   <- sp_pars["z_m"]
v_sp     <- sp_pars["v"]
p_M_sp   <- sp_pars["p_M"]
kap_sp   <- sp_pars["kap"]
del_M_sp <- sp_pars["del_M"]
f_sp     <- sp_pars["f"]
L_ref    <- 1

free_params <- c("z", "v", "p_M", "z_m")
varying <- result[free_params]

d <- read.table("../raw_data/Baskiera_Gvozdik_Dataset_B.txt", header = TRUE)
sex_map <- aggregate(Sex ~ ID, data = d, FUN = function(x) names(which.max(table(x))))
result$individual <- as.integer(result$individual)
result <- merge(result, sex_map, by.x = "individual", by.y = "ID", all.x = TRUE)

# derived quantities
result$p_Am      <- result$z   * L_ref * result$p_M / kap_sp
result$p_Am_m    <- result$z_m * L_ref * result$p_M / kap_sp
result$SVL_max   <- result$z   * L_ref * f_sp / del_M_sp
result$SVL_m_max <- result$z_m * L_ref * f_sp / del_M_sp

normalized <- as.data.frame(lapply(varying, function(x) x / mean(x)))

cv <- sapply(varying, function(x) sd(x) / mean(x) * 100)
cv_sorted <- sort(cv, decreasing = TRUE)
saveRDS(cv_sorted, "../results/param_cv.rds")
normalized <- normalized[names(cv_sorted)]

names(normalized) <- param_latex[names(normalized)]
names(cv_sorted) <- param_latex[names(cv_sorted)]

# figure 1: CV barplot for z, v, p_M, z_m
tex_file <- paste0(fig_folder, "param_zvpM_cv.tex")

X11(width = 6, height = 4)
par(mar = c(4, 4, 4, 2))
bp <- barplot(cv_sorted,
  las = 2, ylab = "CV (\\%)",
  main = "Parameter variability (CV)",
  col = "tomato", ylim = c(0, max(cv_sorted) * 1.15)
)
text(bp, cv_sorted, labels = sprintf("%.1f", cv_sorted), pos = 3, cex = 0.8)
dev.off()

tikz(tex_file, width = 4, height = 3, standAlone = TRUE)
par(mar = c(4, 4, 3, 2))
bp <- barplot(cv_sorted,
  las = 2, ylab = "CV (\\%)",
  col = "tomato", ylim = c(0, max(cv_sorted) * 1.15)
)
text(bp, cv_sorted, labels = sprintf("%.1f\\%%", cv_sorted), pos = 3, cex = 0.6)
dev.off()

old_wd <- getwd()
setwd(fig_folder)
tools::texi2pdf("param_zvpM_cv.tex", clean = TRUE)
setwd(old_wd)
unlink(tex_file)

# figure 2: normalized boxplot
tex_file <- paste0(fig_folder, "param_zvpM_boxplot.tex")

X11(width = 6, height = 4)
par(mar = c(6, 4, 4, 2))
boxplot(normalized,
  las = 2, ylab = "Value / mean",
  col = "steelblue", outpch = 16, outcex = 0.6
)
abline(h = 1, lty = 2, col = "grey50")
dev.off()

tikz(tex_file, width = 4, height = 3, standAlone = TRUE)
par(mar = c(4, 4, 3, 2))
boxplot(normalized,
  las = 2, ylab = "Value / mean",
  col = "steelblue", outpch = 16, outcex = 0.6
)
abline(h = 1, lty = 2, col = "grey50")
dev.off()

old_wd <- getwd()
setwd(fig_folder)
tools::texi2pdf("param_zvpM_boxplot.tex", clean = TRUE)
setwd(old_wd)
unlink(tex_file)

# figure 3: boxplot of derived p_Am by sex
tex_file <- paste0(fig_folder, "explore_zvpM_pAm_boxplot.tex")

X11(width = 8, height = 5)
par(mfrow = c(1, 2))
boxplot(p_Am ~ Sex, data = result,
        ylab = expression(paste({p[Am]}, " (J/d/cm"^2, ")")),
        main = expression(paste({p[Am]}, " by sex")),
        col = c("#440154", "#FDE725"), outpch = 16, outcex = 0.8)
boxplot(p_Am_m ~ Sex, data = result,
        ylab = expression(paste({p[Am]}[m], " (J/d/cm"^2, ")")),
        main = expression(paste({p[Am]}[m], " by sex")),
        col = c("#440154", "#FDE725"), outpch = 16, outcex = 0.8)
dev.off()

tikz(tex_file, width = 5, height = 3, standAlone = TRUE)
par(mfrow = c(1, 2), mar = c(2, 4, 2, 1))
boxplot(p_Am ~ Sex, data = result,
        ylab = "$\\{p_{Am}\\}$ (J\\,d$^{-1}$\\,cm$^{-2}$)",
        col = c("#440154", "#FDE725"), outpch = 16, outcex = 0.5)
boxplot(p_Am_m ~ Sex, data = result,
        ylab = "$\\{p_{Am,m}\\}$ (J\\,d$^{-1}$\\,cm$^{-2}$)",
        col = c("#440154", "#FDE725"), outpch = 16, outcex = 0.5)
dev.off()

old_wd <- getwd()
setwd(fig_folder)
tools::texi2pdf("explore_zvpM_pAm_boxplot.tex", clean = TRUE)
setwd(old_wd)
unlink(tex_file)

# summary statistics
cat("\n")
cat("  CV analysis for z, v, p_M, z_m\n")
cat("\n")

for (nm in names(cv_sorted)) {
  cat(sprintf("  %s: CV = %.2f%%\n", nm, cv_sorted[nm]))
}

cat("\n\nSex-specific summary\n\n")
for (sex in c("F", "M")) {
  sub <- result[result$Sex == sex, ]
  cat(sprintf("--- %s (n = %d) ---\n", ifelse(sex == "F", "Females", "Males"), nrow(sub)))
  for (p in free_params) {
    cat(sprintf("  %s: mean = %.5f  sd = %.5f  CV = %.2f%%\n",
        p, mean(sub[[p]]), sd(sub[[p]]), sd(sub[[p]])/mean(sub[[p]])*100))
  }
  cat(sprintf("  p_Am:   mean = %.1f  sd = %.2f\n", mean(sub$p_Am), sd(sub$p_Am)))
  cat(sprintf("  p_Am_m: mean = %.1f  sd = %.2f\n", mean(sub$p_Am_m), sd(sub$p_Am_m)))
  cat(sprintf("  max SVL (F): %.2f cm\n", mean(sub$SVL_max)))
  cat(sprintf("  max SVL (M): %.2f cm\n", mean(sub$SVL_m_max)))
  cat("\n")
}

cat("Wilcoxon rank-sum tests\n\n")
for (p in free_params) {
  wt <- wilcox.test(result[[p]] ~ Sex, data = result, exact = FALSE)
  cat(sprintf("  %s: W = %.1f, p = %.4f\n", p, wt$statistic, wt$p.value))
}

cat(sprintf("\nComparison to population average\n\n"))
cat(sprintf("Individual z   mean deviates from species by %+.4f (%.1f%%)\n",
    mean(result$z) - z_sp, (mean(result$z) - z_sp) / z_sp * 100))
cat(sprintf("Individual z_m mean deviates from species by %+.4f (%.1f%%)\n",
    mean(result$z_m) - z_m_sp, (mean(result$z_m) - z_m_sp) / z_m_sp * 100))
cat(sprintf("Individual v   mean deviates from species by %+.5f (%.1f%%)\n",
    mean(result$v) - v_sp, (mean(result$v) - v_sp) / v_sp * 100))
cat(sprintf("Individual p_M mean deviates from species by %+.2f (%.1f%%)\n",
    mean(result$p_M) - p_M_sp, (mean(result$p_M) - p_M_sp) / p_M_sp * 100))
