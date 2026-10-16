# shared utility functions for ichthyosaura alpestris DEB modelling

# read parameters from an AmP html report
read_amp_html <- function(html_file) {
  require(xml2)
  doc <- read_html(html_file)
  out <- list()

  for (tr in xml_find_all(doc, "//tr")) {
    id <- xml_attr(tr, "id")
    tds <- xml_find_all(tr, "./td")
    n <- length(tds)
    if (n == 0) next

    # full report
    if (!is.na(id)) {
      if (n >= 3) {
        out[[id]] <- suppressWarnings(type.convert(xml_text(tds[[3]]), as.is = TRUE))
      }
      next
    }

    first <- trimws(xml_text(tds[[1]]))
    if (first %in% c("symbol", "data", "")) next

    if (n >= 3) {
      # results html, parameter row
      if (grepl("^[a-zA-Z_.][a-zA-Z0-9_.]*$", first) && !first %in% c("TRUE", "FALSE", "NA", "NaN")) {
        out[[first]] <- suppressWarnings(type.convert(xml_text(tds[[3]]), as.is = TRUE))
      }
    }

    if (n >= 4) {
      # results html, prediction row
      sym4 <- trimws(xml_text(tds[[4]]))
      if (grepl("^[a-zA-Z_.][a-zA-Z0-9_.]*$", sym4) && !sym4 %in% c("symbol", "data", "prd")) {
        if (is.null(out[[sym4]])) {
          out[[sym4]] <- suppressWarnings(type.convert(xml_text(tds[[2]]), as.is = TRUE))
        }
      }
    }
  }
  out
}

# repeat the last microclimate year nyears_out times
extend_micro <- function(micro_in, nyears_out = 20, nyears_in = NULL) {
  if (is.null(nyears_in)) nyears_in <- micro_in$nyears
  N_HR  <- 365 * 24
  N_DAY <- 365
  last_yr_hr  <- ((nyears_in - 1) * N_HR + 1):(nyears_in * N_HR)
  last_yr_day <- ((nyears_in - 1) * N_DAY + 1):(nyears_in * N_DAY)

  m <- list()
  for (n in names(micro_in)) {
    x <- micro_in[[n]]
    if (identical(n, "dates")) {
      m[[n]] <- seq(0, nyears_out * 365 - 1 / 24, length.out = nyears_out * N_HR)
    } else if (identical(n, "dates2")) {
      m[[n]] <- seq_len(nyears_out * N_DAY)
    } else if (is.matrix(x) && nrow(x) == nyears_in * N_HR) {
      m[[n]] <- do.call(rbind, replicate(nyears_out, x[last_yr_hr, , drop = FALSE], simplify = FALSE))
    } else if (is.vector(x) && length(x) == nyears_in * N_DAY && !is.list(x)) {
      m[[n]] <- rep(x[last_yr_day], nyears_out)
    } else if (is.vector(x) && length(x) == nyears_in * N_HR && !is.list(x)) {
      m[[n]] <- rep(x[last_yr_hr], nyears_out)
    } else {
      m[[n]] <- x
    }
  }
  m$nyears <- nyears_out
  m$ndays <- nyears_out * N_DAY
  m$timeinterval <- 365
  m
}

# DEB parameter list from .mat and html
build_par_list <- function(mat_path, res_html_path, report_html_path,
                           sex = "F", CT_min = -2.5, CT_max = 36) {
  library(R.matlab)
  d <- readMat(mat_path)
  pnames <- dimnames(d$par)[[1]]
  par <- list()
  for (nm in pnames) {
    val <- d$par[which(pnames == nm), , 1][[1]]
    if (length(val) == 1) par[[nm]] <- as.numeric(val)
  }

  z_i  <- par$z
  zm_i <- par$z.m
  v_i <- par$v
  kap_i <- par$kap
  pM_i <- par$p.M
  EG_i <- par$E.G
  kapX_i <- par$kap.X
  kapR_i <- par$kap.R
  kJ_i <- par$k.J
  sG_i <- par$s.G
  ha_i <- par$h.a
  Tref_i <- par$T.ref
  TA_i <- par$T.A
  f_i <- par$f
  dV_i <- par$d.V
  dE_i <- par$d.E
  delM_i <- par$del.M
  EHb_i <- par$E.Hb
  EHp_i <- par$E.Hp
  EHpm_i <- par$E.Hpm
  Ehj_i <- par$E.Hj

  if (is.null(Ehj_i)) Ehj_i <- EHb_i
  if (is.null(zm_i)) zm_i <- z_i
  z_use <- if (sex == "F") z_i else zm_i

  n_CX <- par$n.CX; n_HX <- par$n.HX; n_OX <- par$n.OX; n_NX <- par$n.NX
  n_CE <- par$n.CE; n_HE <- par$n.HE; n_OE <- par$n.OE; n_NE <- par$n.NE
  n_CV <- par$n.CV; n_HV <- par$n.HV; n_OV <- par$n.OV; n_NV <- par$n.NV
  n_CP <- par$n.CP; n_HP <- par$n.HP; n_OP <- par$n.OP; n_NP <- par$n.NP
  n_CN <- par$n.CN; n_HN <- par$n.HN; n_ON <- par$n.ON; n_NN <- par$n.NN

  muX_i <- par$mu.X; muE_i <- par$mu.E; muV_i <- par$mu.V
  muP_i <- par$mu.P; kapP_i <- par$kap.P

  amp_res <- read_amp_html(res_html_path)
  Lj_svl <- amp_res$Lj
  delM_html <- amp_res$del_M
  L.j <- Lj_svl * delM_html

  amp <- read_amp_html(report_html_path)

  z_ratio <- z_use / z_i
  p_Xm_val <- amp$p_Xm * 10 * z_ratio
  p_Am_val <- amp$p_Am * z_ratio
  E_m_val  <- amp$E_m * z_ratio

  w <- (E_m_val / EG_i) * dE_i / dV_i
  Lb_i <- (amp_res$Wwb / (1 + f_i * w))^(1 / 3)

  V_b <- Lb_i^3
  E0_i <- (EHb_i + EG_i * V_b + E_m_val * V_b) / kapR_i

  EHp_use <- if (sex == "F") EHp_i else EHpm_i

  T_L_val  <- if (!is.null(par$T.L))  par$T.L  else (CT_min + 273.15)
  T_H_val  <- if (!is.null(par$T.H))  par$T.H  else (CT_max + 273.15)
  T_AL_val <- if (!is.null(par$T.AL)) par$T.AL else 5e04
  T_AH_val <- if (!is.null(par$T.AH)) par$T.AH else 9e04

  list(
    z         = z_use,
    del_M     = delM_i,
    p_Xm      = p_Xm_val,
    kap_X     = kapX_i,
    v         = v_i / 24,
    kap       = kap_i,
    p_M       = pM_i / 24,
    E_G       = EG_i,
    kap_R     = kapR_i,
    k_J       = kJ_i / 24,
    E_Hb      = EHb_i,
    E_Hj      = Ehj_i,
    E_Hp      = EHp_use,
    E_He      = EHp_use,
    h_a       = ha_i / (24^2),
    s_G       = sG_i,
    T_REF     = Tref_i,
    T_A       = TA_i,
    T_AL      = T_AL_val,
    T_AH      = T_AH_val,
    T_L       = T_L_val,
    T_H       = T_H_val,
    E_0       = E0_i,
    f         = f_i,
    d_V       = dV_i,
    d_E       = dE_i,
    d_Egg     = dE_i,
    mu_X      = muX_i,
    mu_E      = muE_i,
    mu_V      = muV_i,
    mu_P      = muP_i,
    kap_X_P   = kapP_i,
    n_X       = c(n_CX, n_HX, n_OX, n_NX),
    n_E       = c(n_CE, n_HE, n_OE, n_NE),
    n_V       = c(n_CV, n_HV, n_OV, n_NV),
    n_P       = c(n_CP, n_HP, n_OP, n_NP),
    n_M_nitro = c(n_CN, n_HN, n_ON, n_NN),
    L_b       = Lb_i,
    V_init    = L.j^3,
    E_init    = E_m_val,
    E_H_init  = Ehj_i,
    stage     = 1
  )
}

# per-individual E_0 (J per egg)
build_E0_lookup <- function(ind_ids, sp_par_m = NULL, base_dir = "../ind") {
  library(R.matlab)
  rows <- list()
  for (id in ind_ids) {
    mat_path  <- file.path(base_dir, id, "results_Ichthyosaura_alpestris.mat")
    res_html  <- file.path(base_dir, id, "Ichthyosaura_alpestris_res.html")
    rep_html  <- file.path(base_dir, id, "report_Ichthyosaura_alpestris.html")
    if (!all(file.exists(mat_path, res_html, rep_html))) next
    tryCatch({
      d <- readMat(mat_path)
      pnames <- dimnames(d$par)[[1]]
      par <- list()
      for (nm in pnames) {
        val <- d$par[which(pnames == nm), , 1][[1]]
        if (length(val) == 1) par[[nm]] <- as.numeric(val)
      }
      amp_res <- read_amp_html(res_html)
      amp     <- read_amp_html(rep_html)
      E_G     <- par$E.G
      kap_R   <- par$kap.R
      E_Hb    <- par$E.Hb
      dE_dV   <- par$d.E / par$d.V
      E_m     <- amp$p_Am / par$v * par$kap
      E_m_val <- amp$E_m
      f_val   <- par$f
      w       <- (E_m_val / E_G) * dE_dV
      Lb      <- (amp_res$Wwb / (1 + f_val * w))^(1/3)
      V_b     <- Lb^3
      E0      <- (E_Hb + E_G * V_b + E_m_val * V_b) / kap_R
      rows[[length(rows)+1]] <- data.frame(id = as.integer(id), E_0 = E0,
                                           stringsAsFactors = FALSE)
    }, error = function(e) {
      warning(sprintf("build_E0_lookup: failed for id=%s: %s", id, e$message))
    })
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame(id = integer(0), E_0 = numeric(0))
  # fall back to population-average E_0
  if (!is.null(sp_par_m) && nrow(out) > 0) {
    missing_ids <- setdiff(as.character(ind_ids), as.character(out$id))
    if (length(missing_ids) > 0) {
      fb <- data.frame(id = as.integer(missing_ids), E_0 = sp_par_m$E_0,
                       stringsAsFactors = FALSE)
      out <- rbind(out, fb)
    }
  }
  out[order(out$id), ]
}

# colour-scale limits
compute_limits <- function(vals, robust = FALSE, type = c("range", "symmax", "posmax")) {
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2) return(NULL)
  type <- match.arg(type)

  if (type == "range") {
    lim <- if (robust) quantile(vals, probs = c(0.01, 0.99), names = FALSE, na.rm = TRUE)
           else range(vals, na.rm = TRUE, finite = TRUE)
    if (diff(lim) == 0) lim <- lim + c(-1, 1)
    return(lim)
  }

  if (type == "symmax") {
    m <- if (robust) quantile(abs(vals), probs = 0.98, names = FALSE, na.rm = TRUE)
         else max(abs(vals), na.rm = TRUE, finite = TRUE)
    if (m == 0) m <- 1
    return(c(-m, m))
  }

  if (type == "posmax") {
    vals <- vals[is.finite(vals)]
    if (length(vals) < 1) return(NULL)
    m <- if (robust) quantile(vals, probs = 0.98, names = FALSE, na.rm = TRUE)
         else max(vals, na.rm = TRUE, finite = TRUE)
    if (m == 0) m <- 1
    return(m)
  }
}

# metrics from ectotherm output
extract_metrics <- function(ecto) {
  yo <- ecto$yearout
  ys <- ecto$yearsout

  # water deficit
  mb <- as.data.frame(ecto$masbal)
  alive <- which(mb$H2OCut_g > 0)
  if (length(alive) > 1) {
    alive <- alive[alive >= min(which(mb$YEAR == 1 & mb$DAY >= 227))]
    mb_alive <- mb[alive, ]
    mb_alive$date <- paste(mb_alive$YEAR, mb_alive$DAY, sep = "-")
    daily_def <- aggregate(H2OBal_g ~ date, data = mb_alive,
                           FUN = function(x) sum(x[x < 0]))
    TotDeficit   <- abs(sum(daily_def$H2OBal_g, na.rm = TRUE))
    DaysDeficit  <- sum(daily_def$H2OBal_g < 0, na.rm = TRUE)
    MeanDailyDef <- if (DaysDeficit > 0) TotDeficit / DaysDeficit else 0
  } else {
    TotDeficit   <- 0
    DaysDeficit  <- 0
    MeanDailyDef <- 0
  }

  data.frame(
    MaxWgt       = max(ys[, 3], na.rm = TRUE),
    MaxLen       = max(ys[, 4], na.rm = TRUE),
    Tmax         = ys[1, 5],
    Tmin         = ys[1, 6],
    MinRes       = ys[1, 7],
    Bsk          = sum(ys[, 13], na.rm = TRUE),
    Forage       = sum(ys[, 14], na.rm = TRUE),
    Food         = sum(ys[, 16], na.rm = TRUE),
    Fecundity    = yo[1, 7],
    Clutches     = yo[1, 8],
    Pct_Des      = yo[1, 15],
    LifeSpan     = yo[1, 16],
    R0           = yo[1, 18],
    CauseDeath   = ys[1, 23],
    TotDeficit   = TotDeficit,
    DaysDeficit  = DaysDeficit,
    MeanDailyDef = MeanDailyDef
  )
}

# numeric par.xxx values from pars_init.m
parse_species_pars <- function(filepath) {
  lines <- readLines(filepath)
  pattern <- "par\\.(\\w+)\\s*=\\s*([^;]+);"
  matches <- regmatches(lines, gregexpr(pattern, lines))
  items <- lapply(matches, function(m) {
    if (length(m) == 0) return(NULL)
    name  <- sub("par\\.(\\w+)\\s*=\\s*[^;]+;", "\\1", m)
    value <- sub("par\\.\\w+\\s*=\\s*([^;]+);", "\\1", m)
    setNames(trimws(value), name)
  })
  all_items <- unlist(items)
  numeric_mask <- !is.na(suppressWarnings(as.numeric(all_items)))
  setNames(as.numeric(all_items[numeric_mask]), names(all_items)[numeric_mask])
}

# label.xxx strings from pars_init.m
parse_species_labels <- function(filepath) {
  lines <- readLines(filepath)
  pattern <- "label\\.(\\w+)\\s*=\\s*'([^']+)';"
  matches <- regmatches(lines, gregexpr(pattern, lines))
  items <- lapply(matches, function(m) {
    if (length(m) == 0) return(NULL)
    name  <- sub("label\\.(\\w+)\\s*=\\s*'[^']+';", "\\1", m)
    value <- sub("label\\.\\w+\\s*=\\s*'([^']+)';", "\\1", m)
    setNames(value, name)
  })
  unlist(items)
}
