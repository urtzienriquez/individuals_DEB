source("parse_ind_pars.R")
source("utils.R")

library(tikzDevice)
options(tikzLatexPackages = c(
  getOption("tikzLatexPackages"),
  "\\usepackage{underscore}"
))
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

d <- read.table("../raw_data/Baskiera_Gvozdik_Dataset_B.txt", header = TRUE)
sex_map <- aggregate(Sex ~ ID, data = d, FUN = function(x) names(which.max(table(x))))
result$individual <- as.integer(result$individual)
result <- merge(result, sex_map, by.x = "individual", by.y = "ID", all.x = TRUE)

# derived quantities
result$p_Am      <- result$z   * L_ref * result$p_M / kap_sp
result$p_Am_m    <- result$z_m * L_ref * result$p_M / kap_sp
result$SVL_max   <- result$z   * L_ref * f_sp / del_M_sp
result$SVL_m_max <- result$z_m * L_ref * f_sp / del_M_sp

col_sex <- ifelse(result$Sex == "F", "#440154", "#FDE725")
sex_vec <- result$Sex
free_params <- c("z", "v", "p_M", "z_m")

# figure 1: distribution of z, z_m, v, p_M (2x2)
tex_file <- paste0(fig_folder, "explore_zvpM_distribution.tex")

X11(width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(result$z, col = "grey90", border = "grey50", breaks = 10,
     main = expression(z ~ "distribution"),
     xlab = "z (zoom factor)", ylab = "Frequency")
abline(v = z_sp, lty = 2, col = "tomato", lwd = 2)
legend("topleft", legend = sprintf("species = %.4f", z_sp),
       lty = 2, col = "tomato", bty = "n", cex = 0.8)

hist(result$z_m, col = "grey90", border = "grey50", breaks = 8,
     main = expression(z[m] ~ "distribution"),
     xlab = expression(z[m]), ylab = "Frequency")
abline(v = z_m_sp, lty = 2, col = "tomato", lwd = 2)
legend("topleft", legend = sprintf("species = %.4f", z_m_sp),
       lty = 2, col = "tomato", bty = "n", cex = 0.8)

hist(result$v, col = "grey90", border = "grey50", breaks = 10,
     main = expression(dot(v) ~ "distribution"),
     xlab = expression(dot(v) ~ "(cm/d)"), ylab = "Frequency")
abline(v = v_sp, lty = 2, col = "tomato", lwd = 2)
legend("topleft", legend = sprintf("species = %.4f", v_sp),
       lty = 2, col = "tomato", bty = "n", cex = 0.8)

hist(result$p_M, col = "grey90", border = "grey50", breaks = 10,
     main = expression("[" * p[M] * "] distribution"),
     xlab = expression("[" * p[M] * "] (J/d/cm"^3 * ")"), ylab = "Frequency")
abline(v = p_M_sp, lty = 2, col = "tomato", lwd = 2)
legend("topleft", legend = sprintf("species = %.1f", p_M_sp),
       lty = 2, col = "tomato", bty = "n", cex = 0.8)

dev.off()

tikz(tex_file, width = 7, height = 5.5, standAlone = TRUE)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

hist(result$z, col = "grey90", border = "grey50", breaks = 10,
     main = "", xlab = param_latex["z"], ylab = "Frequency")
abline(v = z_sp, lty = 2, col = "tomato", lwd = 1.2)
legend("topleft", legend = sprintf("population average = %.4f", z_sp),
       lty = 2, col = "tomato", bty = "n", cex = 0.9)

hist(result$z_m, col = "grey90", border = "grey50", breaks = 8,
     main = "", xlab = param_latex["z_m"], ylab = "Frequency")
abline(v = z_m_sp, lty = 2, col = "tomato", lwd = 1.2)
legend("topleft", legend = sprintf("population average = %.4f", z_m_sp),
       lty = 2, col = "tomato", bty = "n", cex = 0.9)

hist(result$v, col = "grey90", border = "grey50", breaks = 10,
     main = "", xlab = param_latex["v"], ylab = "Frequency")
abline(v = v_sp, lty = 2, col = "tomato", lwd = 1.2)
legend("topleft", legend = sprintf("population average = %.4f", v_sp),
       lty = 2, col = "tomato", bty = "n", cex = 0.9)

hist(result$p_M, col = "grey90", border = "grey50", breaks = 10,
     main = "", xlab = param_latex["p_M"], ylab = "Frequency")
abline(v = p_M_sp, lty = 2, col = "tomato", lwd = 1.2)
legend("topleft", legend = sprintf("population average = %.1f", p_M_sp),
       lty = 2, col = "tomato", bty = "n", cex = 0.9)

dev.off()

old_wd <- getwd()
setwd(fig_folder)
tools::texi2pdf("explore_zvpM_distribution.tex", clean = TRUE)
setwd(old_wd)
unlink(tex_file)

# figure 2: boxplots of z, z_m, v, p_M by sex (2x2)
tex_file <- paste0(fig_folder, "explore_zvpM_boxplot_sex.tex")

X11(width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (p in free_params) {
  boxplot(result[[p]] ~ Sex, data = result,
          ylab = p, main = paste(p, "by sex"),
          col = c("#440154", "#FDE725"), outpch = 16, outcex = 0.8)
  stripchart(result[[p]] ~ Sex, data = result, vertical = TRUE, add = TRUE,
             method = "jitter", pch = 16, col = "#00000060")
}
dev.off()

tikz(tex_file, width = 5.5, height = 5, standAlone = TRUE)
par(mfrow = c(2, 2), mar = c(2, 4, 2, 1))
for (p in free_params) {
  boxplot(result[[p]] ~ Sex, data = result,
          ylab = param_latex[p],
          col = c("#440154", "#FDE725"), outpch = 16, outcex = 0.5)
  stripchart(result[[p]] ~ Sex, data = result, vertical = TRUE, add = TRUE,
             method = "jitter", pch = 16, col = "#00000060", cex = 0.5)
}
dev.off()

old_wd <- getwd()
setwd(fig_folder)
tools::texi2pdf("explore_zvpM_boxplot_sex.tex", clean = TRUE)
setwd(old_wd)
unlink(tex_file)

tex_file <- paste0(fig_folder, "explore_boxplot_sex.tex")

tikz(tex_file, width = 5.5, height = 5, standAlone = TRUE)
par(mfrow = c(2, 2), mar = c(2, 4, 2, 1))
for (p in free_params) {
  boxplot(result[[p]] ~ Sex, data = result,
          ylab = param_latex[p],
          col = c("#440154", "#FDE725"), outpch = 16, outcex = 0.5)
  stripchart(result[[p]] ~ Sex, data = result, vertical = TRUE, add = TRUE,
             method = "jitter", pch = 16, col = "#00000060", cex = 0.5)
}
dev.off()

old_wd <- getwd()
setwd(fig_folder)
tools::texi2pdf("explore_boxplot_sex.tex", clean = TRUE)
setwd(old_wd)
unlink(tex_file)

# figure 3: scatterplot matrix (z, z_m, v, p_M)
tex_file <- paste0(fig_folder, "explore_zvpM_scatter_matrix.tex")

X11(width = 8, height = 8)
pairs(~ z + v + p_M + z_m, data = result,
      col = col_sex, pch = 19, cex = 0.8,
      labels = c("z", "v", "p[M]", "z[m]"),
      upper.panel = function(x, y, ...) {
        points(x, y, col = col_sex, pch = 19, cex = 0.8)
        usr <- par("usr")
        text(usr[2] - 0.15 * diff(usr[1:2]),
             usr[3] + 0.15 * diff(usr[3:4]),
             sprintf("r = %.3f", cor(x, y)), adj = c(1, 0), cex = 0.7)
      },
      diag.panel = NULL,
      lower.panel = function(x, y, ...) {
        points(x, y, col = col_sex, pch = 19, cex = 0.8)
      })
dev.off()

tikz(tex_file, width = 5.5, height = 5.5, standAlone = TRUE)
panel.cor <- function(x, y, ...) {
  usr <- par("usr")
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- sprintf("$r = %.3f$", r)
  text(0.5, 0.5, txt, cex = 0.8)
}
panel.pts <- function(x, y, ...) {
  points(x, y, col = col_sex, pch = 19, cex = 0.6)
}
pairs(~ z + v + p_M + z_m, data = result,
      diag.panel = NULL,
      upper.panel = panel.cor,
      lower.panel = panel.pts,
      labels = c(param_latex["z"], param_latex["v"],
                 param_latex["p_M"], param_latex["z_m"]),
      gap = 0.4)
dev.off()

old_wd <- getwd()
setwd(fig_folder)
tools::texi2pdf("explore_zvpM_scatter_matrix.tex", clean = TRUE)
setwd(old_wd)
unlink(tex_file)

# figure 4: z (females) and z_m (males) per individual
tex_file <- paste0(fig_folder, "explore_zvpM_barplot_z.tex")

females <- result[result$Sex == "F", ]
males   <- result[result$Sex == "M", ]
ord_f <- order(females$z)
ord_m <- order(males$z_m)
fem_ord <- females[ord_f, ]
mal_ord <- males[ord_m, ]

ymax_z  <- max(result$z)  * 1.15
ymax_zm <- max(result$z_m) * 1.15
ymax <- max(ymax_z, ymax_zm)

X11(width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))

bp_f <- barplot(fem_ord$z, names.arg = fem_ord$individual,
                las = 2, col = "#440154", border = NA, cex.names = 0.7,
                main = expression("a) Female zoom factor " * z),
                ylab = "Zoom factor", ylim = c(0, ymax))
abline(h = z_sp, lty = 2, col = "tomato", lwd = 1.5)
text(mean(bp_f), z_sp, labels = sprintf("species = %.4f", z_sp),
     pos = 3, col = "tomato", cex = 0.7)

bp_m <- barplot(mal_ord$z_m, names.arg = mal_ord$individual,
                las = 2, col = "#FDE725", border = NA, cex.names = 0.7,
                main = expression("b) Male zoom factor " * z[m]),
                ylab = "Zoom factor", ylim = c(0, ymax))
abline(h = z_m_sp, lty = 2, col = "tomato", lwd = 1.5)
text(mean(bp_m), z_m_sp, labels = sprintf("species = %.4f", z_m_sp),
     pos = 3, col = "tomato", cex = 0.7)
dev.off()

tikz(tex_file, width = 7, height = 3, standAlone = TRUE)
par(mfrow = c(1, 2), mar = c(2, 4, 2, 1))

bp_f <- barplot(fem_ord$z, names.arg = fem_ord$individual,
                las = 2, col = "#440154", border = NA, cex.names = 0.45,
                ylab = "Zoom factor", ylim = c(0, ymax))
abline(h = z_sp, lty = 2, col = "tomato", lwd = 1)
text(mean(bp_f), z_sp, labels = sprintf("species = %.4f", z_sp),
     pos = 3, col = "tomato", cex = 0.5)

bp_m <- barplot(mal_ord$z_m, names.arg = mal_ord$individual,
                las = 2, col = "#FDE725", border = NA, cex.names = 0.45,
                ylab = "Zoom factor", ylim = c(0, ymax))
abline(h = z_m_sp, lty = 2, col = "tomato", lwd = 1)
text(mean(bp_m), z_m_sp, labels = sprintf("species = %.4f", z_m_sp),
     pos = 3, col = "tomato", cex = 0.5)
dev.off()

old_wd <- getwd()
setwd(fig_folder)
tools::texi2pdf("explore_zvpM_barplot_z.tex", clean = TRUE)
setwd(old_wd)
unlink(tex_file)

# summary statistics
cat("\n")
cat("  z, v, p_M, z_m parameter summary\n")
cat("\n")

for (p in free_params) {
  cat(sprintf("%s:\n", p))
  cat(sprintf("  mean = %.5f   sd = %.5f   CV = %.2f%%\n",
      mean(result[[p]]), sd(result[[p]]), sd(result[[p]]) / mean(result[[p]]) * 100))
  cat(sprintf("  range = [%.5f, %.5f]\n", min(result[[p]]), max(result[[p]])))
  cat(sprintf("  species = %.5f\n\n", get(paste0(p, "_sp"))))
}

cat("Derived quantities:\n")
cat(sprintf("  p_Am (F): range [%.1f, %.1f] J/d/cm^2\n",
    min(result$p_Am), max(result$p_Am)))
cat(sprintf("  p_Am (M): range [%.1f, %.1f] J/d/cm^2\n",
    min(result$p_Am_m), max(result$p_Am_m)))
cat(sprintf("  max SVL (F): range [%.2f, %.2f] cm\n",
    min(result$SVL_max), max(result$SVL_max)))
cat(sprintf("  max SVL (M): range [%.2f, %.2f] cm\n\n",
    min(result$SVL_m_max), max(result$SVL_m_max)))

cat("Wilcoxon rank-sum tests (sex differences):\n")
for (p in free_params) {
  wt <- wilcox.test(result[[p]] ~ Sex, data = result, exact = FALSE)
  cat(sprintf("  %s: W = %.1f, p = %.4f\n", p, wt$statistic, wt$p.value))
}

cat(sprintf("\nCorrelation matrix:\n"))
cor_m <- cor(result[free_params])
print(round(cor_m, 4))
