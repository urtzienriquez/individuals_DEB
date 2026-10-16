files <- list.files("../ind", pattern = "pars_init_Ichthyosaura_alpestris\\.m$",
                    recursive = TRUE, full.names = TRUE)
files <- files[grepl("/\\d+/", files)]

parse_m_file <- function(file) {
  ind_id <- basename(dirname(file))
  lines <- readLines(file)
  matches <- regmatches(lines, gregexpr("par\\.(\\w+)\\s*=\\s*([+-]?\\d*\\.?\\d+(?:[eE][+-]?\\d+)?)\\s*;", lines))
  params <- lapply(matches, function(m) {
    if (length(m) == 0) return(NULL)
    name  <- gsub("par\\.(\\w+)\\s*=\\s*[^;]+;", "\\1", m)
    value <- gsub("par\\.\\w+\\s*=\\s*([^;]+);", "\\1", m)
    setNames(as.numeric(value), name)
  })
  c(individual = ind_id, unlist(params))
}

parsed <- lapply(files, parse_m_file)
result <- as.data.frame(do.call(rbind, parsed), stringsAsFactors = FALSE)
result[-1] <- lapply(result[-1], as.numeric)

param_latex <- c(
  z     = "$z$",
  v     = "$\\dot{v}$",
  kap   = "$\\kappa$",
  p_M   = "$[p_M]$",
  E_G   = "$[E_G]$",
  E_Hb  = "$E_H^b$",
  E_Hp  = "$E_H^p$",
  h_a   = "$\\ddot{h}_a$",
  E_Hj  = "$E_H^j$",
  E_Hpm = "$E_{H,m}^{p}$",
  del_M = "$\\delta_M$",
  z_m   = "$z_m$"
)
