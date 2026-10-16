repo_root <- normalizePath(file.path(getwd(), ".."))
ind_dir <- file.path(repo_root, "ind")
sp_dir  <- file.path(repo_root, "sp")

unlink(ind_dir, recursive = TRUE)
dir.create(ind_dir)

# predict templates
sp_p <- readLines(file.path(sp_dir, "predict_Ichthyosaura_alpestris.m"))

# female template
p_f <- sp_p

idx <- grep("  tvel = get_tp(pars_tpm, f_tW, get_lp(pars_tjm, f_tW), tWw_m", p_f, fixed = TRUE)
end_mw <- grep("  EWw_m = (L_mm * tvel(:,4)).^3 .* (1 + ome_m * tvel(:,3));", p_f, fixed = TRUE)
end_mw <- end_mw[end_mw > idx[1]][1]
p_f <- p_f[-( (idx[1] - 1) : end_mw )]

idx <- grep("  p_refm = p_Am_m * L_mm^2;", p_f, fixed = TRUE)
end_mj <- grep("  EJO_m = - 1e3 * J_M(3,:)'/ 24/ X_gas;", p_f, fixed = TRUE)
end_mj <- end_mj[end_mj > idx[1]][1]
p_f <- p_f[-( (idx[1] - 1) : end_mj )]

p_f <- p_f[!grepl("prdData.tWw_m", p_f, fixed = TRUE)]
p_f <- p_f[!grepl("prdData.tJO_m", p_f, fixed = TRUE)]

p_f <- gsub("temp.tWw_f", "temp.tWw", p_f, fixed = TRUE)
p_f <- gsub("tWw_f", "tWw", p_f, fixed = TRUE)
p_f <- gsub("EWw_f", "EWw", p_f, fixed = TRUE)
p_f <- gsub("tJO_f", "tJO", p_f, fixed = TRUE)
p_f <- gsub("EJO_f", "EJO", p_f, fixed = TRUE)

writeLines(p_f, file.path(ind_dir, "predict_Ichthyosaura_alpestris_female.m"))

# male template
p_m <- sp_p

idx <- grep("  tvel = get_tp(pars_tp, f_tW, get_lp(pars_tj, f_tW), tWw_f", p_m, fixed = TRUE)
end_fw <- grep("  EWw_f = (L_m * tvel(:,4)).^3 .* (1 + ome * tvel(:,3));", p_m, fixed = TRUE)
end_fw <- end_fw[end_fw > idx[1]][1]
p_m <- p_m[-( (idx[1] - 1) : end_fw )]

idx <- grep("  p_ref = p_Am * L_m^2;", p_m, fixed = TRUE)
end_fj <- grep("  EJO_f = - 1e3 * J_M(3,:)'/ 24/ X_gas;", p_m, fixed = TRUE)
end_fj <- end_fj[end_fj > idx[1]][1]
p_m <- p_m[-( (idx[1] - 1) : end_fj )]

p_m <- p_m[!grepl("prdData.tWw_f", p_m, fixed = TRUE)]
p_m <- p_m[!grepl("prdData.tJO_f", p_m, fixed = TRUE)]

p_m <- gsub("temp.tWw_f", "temp.tWw", p_m, fixed = TRUE)
p_m <- gsub("tWw_m", "tWw", p_m, fixed = TRUE)
p_m <- gsub("EWw_m", "EWw", p_m, fixed = TRUE)
p_m <- gsub("tJO_m", "tJO", p_m, fixed = TRUE)
p_m <- gsub("EJO_m", "EJO", p_m, fixed = TRUE)

writeLines(p_m, file.path(ind_dir, "predict_Ichthyosaura_alpestris_male.m"))

# mydata template
sp_m <- readLines(file.path(sp_dir, "mydata_Ichthyosaura_alpestris.m"))

tL_end <- grep("comment.tL_m = 'Data for males'", sp_m, fixed = TRUE)
pre_tWw_end <- tL_end + 1

weights_start <- grep("%% set weights for all real data", sp_m, fixed = TRUE)
links_start <- grep("%% Links", sp_m, fixed = TRUE)

header <- sp_m[1:pre_tWw_end]
post_data <- sp_m[weights_start:(links_start - 1)]
links_refs <- sp_m[links_start:length(sp_m)]

w_idx <- grep("weights.tWw_f", post_data, fixed = TRUE)
wm_idx <- grep("weights.tWw_m", post_data, fixed = TRUE)
post_data[w_idx] <- "weights.tWw = 5 * weights.tWw;"
post_data <- post_data[-wm_idx]

s2_idx <- grep("set2 = ", post_data, fixed = TRUE)
post_data[s2_idx] <- "set2 = {'tWw'}; subtitle2 = {'Wet weight'};"
s3_idx <- grep("set3 = ", post_data, fixed = TRUE)
post_data[s3_idx] <- "set3 = {'tJO'}; subtitle3 = {'O2 consumption'};"

post_data <- gsub("tWw_f", "tWw", post_data, fixed = TRUE)
post_data <- gsub("tWw_m", "tWw", post_data, fixed = TRUE)
post_data <- gsub("tJO_f", "tJO", post_data, fixed = TRUE)
post_data <- gsub("tJO_m", "tJO", post_data, fixed = TRUE)

empty_tWw <- c(
  "% time-weight",
  "data.tWw = [ ... % time since metam (d), wet weight (g)",
  "];",
  "units.tWw = {'d', 'g'}; label.tWw = {'time since metam', 'Wet weight'};",
  "temp.tWw = C2K(12); units.temp.tWw = 'K'; label.temp.tWw = 'temperature';",
  "comment.tWw = 'Individual data'; bibkey.tWw = 'BaskGvoz2021';"
)

empty_tJO <- c(
  "% time-respiration",
  "data.tJO = [ ... % time since birth (d), O2 consumption (mL/h)",
  "];",
  "units.tJO  = {'d', 'mL/h'};  label.tJO = {'time since birth', 'O2 consumption'};",
  "temp.tJO   = C2K(12);  units.temp.tJO = 'K'; label.temp.tJO = 'temperature';",
  "bibkey.tJO = 'BaskGvoz2021';",
  "comment.tJO = 'Individual data';"
)

template <- c(header, empty_tWw, "", empty_tJO, "", post_data, links_refs)
template_path <- file.path(ind_dir, "mydata_Ichthyosaura_alpestris_template.m")
writeLines(template, template_path)

# split into temp_1, temp_2, temp_3
tw_line <- grep("data.tWw = \\[ \\.\\.\\. % time since metam \\(d\\), wet weight \\(g\\)", template)
temp1 <- template[1:tw_line]
writeLines(temp1, file.path(ind_dir, "temp_1.m"))

jo_line <- grep("data.tJO = \\[ \\.\\.\\. % time since birth \\(d\\), O2 consumption \\(mL/h\\)", template)
temp2 <- template[(tw_line + 2):jo_line]
temp2 <- template[(tw_line + 1):jo_line]
writeLines(temp2, file.path(ind_dir, "temp_2.m"))

temp3 <- template[(jo_line + 1):length(template)]
writeLines(temp3, file.path(ind_dir, "temp_3.m"))

# copy run.m and pars_init.m
file.copy(file.path(sp_dir, "run_Ichthyosaura_alpestris.m"), ind_dir)
file.copy(file.path(sp_dir, "pars_init_Ichthyosaura_alpestris.m"), ind_dir)

cat("Done. Files created in:", ind_dir, "\n")
cat("  predict_Ichthyosaura_alpestris_female.m\n")
cat("  predict_Ichthyosaura_alpestris_male.m\n")
cat("  mydata_Ichthyosaura_alpestris_template.m\n")
cat("  temp_1.m / temp_2.m / temp_3.m\n")
cat("  run_Ichthyosaura_alpestris.m\n")
cat("  pars_init_Ichthyosaura_alpestris.m\n")
