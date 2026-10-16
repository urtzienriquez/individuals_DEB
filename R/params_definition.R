# all parameters for ichthyosaura alpestris DEB model

library(NicheMapR)
library(R.matlab)
source("utils.R")

# constants
CT_min      <- -2
CT_max      <- 36
clutch_ab   <- c(110, 0)
minclutch   <- 300
maxclutch   <- 600
batch       <- 1
cond_clutch <- 1
photostart  <- 4   # spring equinox (doy 80)
photofinish <- 1   # summer solstice (doy 173)
startday    <- 227
X           <- 10  # food density
raindrink   <- 0
z.mult      <- 1

# behav_stages matrix
stages <- 8

behav_stages <- matrix(c(
  1, 1, 1, 1, 0, 2, 10, 1, 0, 1, 0, 0.1, 0, 0, 0, # stage 0
  1, 1, 1, 1, 0, 2, 10, 1, 0, 1, 0, 0.1, 0, 0, 0, # stage 1
  1, 1, 1, 1, 0, 2, 10, 1, 0, 1, 0, 0.1, 0, 0, 0, # stage 2
  1, 1, 1, 1, 0, 2, 10, 1, 0, 1, 0, 0.1, 0, 0, 0, # stage 3
  0, 1, 1, 1, 0, 2, 10, 1, 0, 1, 0, 0.1, 0, 0, 0, # stage 4
  0, 1, 1, 1, 0, 2, 10, 1, 0, 1, 0, 0.1, 0, 0, 0, # stage 5
  0, 1, 1, 1, 0, 2, 10, 1, 0, 1, 0, 0.1, 0, 0, 0, # stage 6
  0, 1, 1, 1, 0, 2, 10, 1, 0, 1, 0, 0.1, 0, 0, 0  # stage 7
), nrow = 8, ncol = 15, byrow = TRUE)

# biophysical parameters
bio_par <- list(
  DEB          = 1,
  CT_kill = 1,
  CT_minthresh = 5,
  behav_stages = behav_stages,
  pct_H_R      = 50,
  pct_H_death  = 100,
  viviparous   = 0,
  shape        = 3,
  alpha_max    = 0.85,
  alpha_min    = 0.85,
  T_F_min      = 4,
  T_F_max      = 28,
  T_B_min      = 4,
  T_RB_min     = 4,
  T_pref       = 17,
  CT_max       = CT_max,
  CT_min       = CT_min,
  diurn        = 1,
  nocturn      = 1,
  crepus       = 1,
  shade_seek   = 1,
  shdburrow    = 1,
  burrow       = 1,
  climb        = 0,
  mindepth     = 2,
  maxdepth     = 10,
  pct_wet      = 90,
  photostart   = photostart,
  photofinish  = photofinish,
  startday     = startday,
  X            = X,
  raindrink    = raindrink,
  z.mult       = z.mult
)

# sex lookup
dataset <- read.table("../raw_data/Baskiera_Gvozdik_Dataset_B.txt",
  header = TRUE, stringsAsFactors = FALSE
)
sex_lookup <- aggregate(Sex ~ ID,
  data = dataset,
  FUN = function(x) names(sort(table(x), decreasing = TRUE))[1]
)
sex_map <- setNames(sex_lookup$Sex, as.character(sex_lookup$ID))

# population-average, females
cat("Building population-average parameters ...\n")
flush.console()

sp_par_f <- build_par_list(
  mat_path         = "../sp/results_Ichthyosaura_alpestris.mat",
  res_html_path    = "../sp/Ichthyosaura_alpestris_res.html",
  report_html_path = "../sp/report_Ichthyosaura_alpestris.html",
  sex = "F", CT_min = CT_min, CT_max = CT_max
)
sp_par_f$clutch_ab    <- clutch_ab
sp_par_f$minclutch    <- minclutch
sp_par_f$maxclutch    <- maxclutch
sp_par_f$batch        <- batch
sp_par_f$cond_clutch  <- cond_clutch

sp_par_m <- build_par_list(
  mat_path         = "../sp/results_Ichthyosaura_alpestris.mat",
  res_html_path    = "../sp/Ichthyosaura_alpestris_res.html",
  report_html_path = "../sp/report_Ichthyosaura_alpestris.html",
  sex = "M", CT_min = CT_min, CT_max = CT_max
)
# males: continuous reproduction
sp_par_m$clutch_ab    <- c(0, 0)
sp_par_m$minclutch    <- 1
sp_par_m$maxclutch    <- 1
sp_par_m$batch        <- 0
sp_par_m$cond_clutch  <- 0

# individuals
ind_ids <- list.dirs("../ind", recursive = FALSE, full.names = FALSE)
ind_ids <- sort(as.numeric(ind_ids))
ind_ids <- ind_ids[!is.na(ind_ids)]

ind_pars <- list()
for (id in ind_ids) {
  sex <- sex_map[as.character(id)]
  if (is.na(sex)) {
    cat(sprintf("  WARNING: ID %d has unknown sex, defaulting to F\n", id))
    sex <- "F"
  }
  p <- build_par_list(
    mat_path         = sprintf("../ind/%d/results_Ichthyosaura_alpestris.mat", id),
    res_html_path    = sprintf("../ind/%d/Ichthyosaura_alpestris_res.html", id),
    report_html_path = sprintf("../ind/%d/report_Ichthyosaura_alpestris.html", id),
    sex = sex, CT_min = CT_min, CT_max = CT_max
  )
  p$id  <- id
  p$sex <- sex
  if (sex == "M") {
    p$clutch_ab    <- c(0, 0)
    p$minclutch    <- 1
    p$maxclutch    <- 1
    p$batch        <- 0
    p$cond_clutch  <- 0
  } else {
    p$clutch_ab    <- clutch_ab
    p$minclutch    <- minclutch
    p$maxclutch    <- maxclutch
    p$batch        <- batch
    p$cond_clutch  <- cond_clutch
  }
  ind_pars[[length(ind_pars) + 1]] <- p
}
cat(sprintf("  Built %d individual parameter sets\n", length(ind_pars)))
flush.console()
