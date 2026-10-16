repo_root <- normalizePath(file.path(getwd(), ".."))
database <- read.table(file.path(repo_root, "raw_data", "Baskiera_Gvozdik_Dataset_B.txt"), header = TRUE)
database$date <- as.Date(database$date, tryFormats = c("%d.%m.%Y"))
str(database)

ids <- unique(database$ID)

# fix how julian date increases with year
for (i in ids) {
  db <- database[database$ID == i, ]
  db$date.julian <- julian(db$date, origin = as.Date("2014/01/01"))

  ind_dir <- file.path(repo_root, "ind", i)
  dir.create(ind_dir)
  setwd(ind_dir)

  sink(paste0("mydata_Ichthyosaura_alpestris.m"))
  cat(readLines(file.path(repo_root, "ind", "temp_1.m")), sep = "\n")
  write.table(unname(as.matrix(db[, c("date.julian", "BM")])), row.names = FALSE, col.names = FALSE)
  cat(readLines(file.path(repo_root, "ind", "temp_2.m")), sep = "\n")
  write.table(unname(as.matrix(db[, c("date.julian", "VO2")])), row.names = FALSE, col.names = FALSE)
  cat(readLines(file.path(repo_root, "ind", "temp_3.m")), sep = "\n")
  sink()

  if (length(unique(db$Sex)) > 1) {
    sx <- names(which.max(table(db$Sex)))
    db$Sex <- sx
  }

  if (unique(db$Sex) == "F") {
    file.copy(file.path(repo_root, "ind", "predict_Ichthyosaura_alpestris_female.m"), "./predict_Ichthyosaura_alpestris.m")
  } else if (unique(db$Sex) == "M") {
    file.copy(file.path(repo_root, "ind", "predict_Ichthyosaura_alpestris_male.m"), "./predict_Ichthyosaura_alpestris.m")
  }

  file.copy(file.path(repo_root, "sp", "run_Ichthyosaura_alpestris.m"), "./run_Ichthyosaura_alpestris.m")
  file.copy(file.path(repo_root, "sp", "pars_init_Ichthyosaura_alpestris.m"), "./pars_init_Ichthyosaura_alpestris.m")

  setwd(repo_root)
}
