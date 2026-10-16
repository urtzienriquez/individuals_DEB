#!/usr/bin/env Rscript
# european grid simulation (chunked for cluster)

TIMEINTERVAL  <- 365
NYEARS        <- 3
NYEARS_EXT    <- 20
N_CORES       <- 10
EU_EXT        <- c(-10, 40, 35, 70)

args <- commandArgs(trailingOnly = TRUE)
CHUNK <- NULL
N_CHUNKS <- NULL
WARM <- NULL
SCENARIO <- NULL

for (a in args) {
  if (grepl("^--chunk=", a)) {
    val <- sub("^--chunk=", "", a)
    parts <- strsplit(val, "/")[[1]]
    CHUNK     <- as.integer(parts[1])
    N_CHUNKS  <- as.integer(parts[2])
  } else if (grepl("^--warm=", a)) {
    WARM <- as.numeric(sub("^--warm=", "", a))
  } else if (grepl("^--scenario=", a)) {
    SCENARIO <- sub("^--scenario=", "", a)
  }
}

if (is.null(CHUNK) || is.null(WARM) || is.null(SCENARIO)) {
  stop("Usage: Rscript 3.1_run_europe_chunk.R --chunk=M/N --warm=0|2|4 --scenario=baseline|w2|w4")
}

OUTPUT_FILE <- sprintf("../data/europe_results_%s_chunk%dof%d.rds", SCENARIO, CHUNK, N_CHUNKS)

library(NicheMapR)
library(terra)
library(future)
library(future.apply)

source("utils.R")
source("params_definition.R")

# simulate one grid cell
run_cell <- function(lon, lat, elev, sp_par_f, sp_par_m, ind_pars, bio_par,
                     timeinterval, nyears, nyears_ext, warm) {
  micro <- micro_global(
    loc = c(lon, lat), elev = elev,
    timeinterval = timeinterval, nyears = nyears,
    runmoist = 1, snowmodel = 1,
    run.gads = 0, warm = warm
  )
  micro <- extend_micro(micro, nyears_out = nyears_ext)
  assign("micro", micro, .GlobalEnv)

  ecto_sp_f <- do.call(ectotherm, c(bio_par, sp_par_f))
  sp_res_f <- extract_metrics(ecto_sp_f)
  out_sp_f <- data.frame(
    type = "species", id = "sp", sex = "F",
    lon = lon, lat = lat, elev = elev, sp_res_f
  )

  ecto_sp_m <- do.call(ectotherm, c(bio_par, sp_par_m))
  sp_res_m <- extract_metrics(ecto_sp_m)
  out_sp_m <- data.frame(
    type = "species", id = "sp", sex = "M",
    lon = lon, lat = lat, elev = elev, sp_res_m
  )

  ind_rows <- lapply(seq_along(ind_pars), function(j) {
    p <- ind_pars[[j]]
    ecto_args <- p[!names(p) %in% c("id", "sex")]
    ecto <- do.call(ectotherm, c(bio_par, ecto_args))
    met <- extract_metrics(ecto)
    data.frame(
      type = "individual", id = p$id, sex = p$sex,
      lon = lon, lat = lat, elev = elev, met
    )
  })
  out_ind <- do.call(rbind, ind_rows)

  rbind(out_sp_f, out_sp_m, out_ind)
}

# build europe grid
cat(sprintf("Building Europe grid (chunk %d / %d, scenario: %s, warm: %d)\n",
            CHUNK, N_CHUNKS, SCENARIO, WARM))

gcfolder_file <- paste0(.libPaths()[1], "/gcfolder.rda")
if (!file.exists(gcfolder_file)) {
  stop("Global climate data not found. Run get.global.climate() first.")
}
load(gcfolder_file)

global_rast <- terra::rast(paste0(folder, "/global_climate.nc"))
eu_ext <- ext(EU_EXT[1], EU_EXT[2], EU_EXT[3], EU_EXT[4])
eu_rast <- crop(global_rast[[1]], eu_ext)

lons <- xFromCol(eu_rast)
lats <- yFromRow(eu_rast)
elev <- values(eu_rast)

grid <- expand.grid(lon = lons, lat = lats)
grid$elev <- as.vector(elev)
grid <- grid[!is.na(grid$elev), ]
rownames(grid) <- NULL

total <- nrow(grid)
chunk_size <- ceiling(total / N_CHUNKS)
start_idx <- (CHUNK - 1L) * chunk_size + 1L
end_idx   <- min(CHUNK * chunk_size, total)
grid_chunk <- grid[start_idx:end_idx, ]

cat(sprintf("  Total grid cells: %d\n", total))
cat(sprintf("  This chunk: cells %d to %d (%d cells)\n", start_idx, end_idx, nrow(grid_chunk)))

# run
cat(sprintf("  Cores: %d\n", N_CORES))
cat(sprintf("  Output: %s\n\n", OUTPUT_FILE))

plan(multisession, workers = N_CORES)

results_list <- future_lapply(seq_len(nrow(grid_chunk)), function(i) {
  gi <- start_idx + i - 1L
  tryCatch(
    {
      r <- run_cell(
        lon          = grid_chunk$lon[i],
        lat          = grid_chunk$lat[i],
        elev         = grid_chunk$elev[i],
        sp_par_f     = sp_par_f,
        sp_par_m     = sp_par_m,
        ind_pars     = ind_pars,
        bio_par      = bio_par,
        timeinterval = TIMEINTERVAL,
        nyears       = NYEARS,
        nyears_ext   = NYEARS_EXT,
        warm         = WARM
      )
      if (i %% 500 == 0) cat(sprintf("    cell %d / %d\n", i, nrow(grid_chunk)))
      r
    },
    error = function(e) {
      cat(sprintf(
        "    CELL %d FAILED (lon=%.1f lat=%.1f): %s\n",
        gi, grid_chunk$lon[i], grid_chunk$lat[i], conditionMessage(e)
      ))
      data.frame(
        type = "error", id = NA, sex = NA,
        lon = grid_chunk$lon[i], lat = grid_chunk$lat[i], elev = grid_chunk$elev[i],
        MaxWgt = NA, MaxLen = NA, Tmax = NA, Tmin = NA, MinRes = NA,
        Bsk = NA, Forage = NA, Food = NA, Fecundity = NA, Clutches = NA,
        Pct_Des = NA, LifeSpan = NA, R0 = NA, CauseDeath = NA,
        TotDeficit = NA, DaysDeficit = NA, MeanDailyDef = NA
      )
    }
  )
},
future.packages = "NicheMapR",
future.seed = TRUE,
future.chunk.size = 5
)

# save
cat("  Saving chunk results ...\n")
results <- do.call(rbind, results_list)
rownames(results) <- NULL

dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)
saveRDS(results, OUTPUT_FILE)

n_ok <- sum(results$type != "error" | is.na(results$type))
cat(sprintf("\nDone. %d rows written to %s\n", nrow(results), OUTPUT_FILE))
cat(sprintf("  Successful cells: %d / %d\n", n_ok, nrow(grid_chunk)))
