# =============================================================================
# Quick screening: a few responses x a few structures, settings through environment variables.
# Same models and scores as analysis/model_comparison.qmd, but sequential and without figures.
# Input : parcels_with_soil_complete.gpkg  (one row per field x year, soil columns NA when not sampled)
# Output: <OUTDIR>/results.csv        one row per variable x fixed set x structure (criteria + hold-out scores)
#         <OUTDIR>/holdout_pred.csv   hold-out predictions (for observed-vs-predicted plots)
#         <OUTDIR>/hyper.csv          hyperparameters of the full-data fits
# Usage : Rscript R/quick_screen.R          (run from the repository root; settings below)
# Example: VARS=pH,OC MODELS="S02|S08|S14" Rscript R/quick_screen.R
# =============================================================================
suppressPackageStartupMessages({ library(sf); library(Matrix); library(INLA) })

# ---- 1. SETTINGS ------------------------------------------------------------------------------
env <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
GPKG    <- env("GPKG", "data/parcels_with_soil_complete.gpkg")
OUTDIR  <- env("OUTDIR", "outputs/quick_screen")
VARS    <- strsplit(env("VARS", "pH,OC,log_P2O5,Ca,ILR3_1,ILR3_2"), ",")[[1]]
FIXED   <- strsplit(env("FIXED", "1;Bloc + irrg_ty"), ";")[[1]]
MODELS  <- env("MODELS", "all")
SCHEMES <- strsplit(env("SCHEMES", "row,field"), ",")[[1]]
N_REP   <- as.integer(env("N_REP", 1)); HOLD <- 0.10; SEED <- 123
DIST    <- as.numeric(env("DIST", 110))         # neighbourhood: boundary-to-boundary distance (m)
STRAT   <- env("INT", "eb")                     # "eb" (fast screening) or "ccd" (full integration)
inla.setOption(num.threads = as.integer(env("NTHREADS", 1)))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
TLIM <- as.integer(env("INLA_TIMEOUT", 0))        # optional time limit per fit (seconds, Linux `timeout`); 0 = none
if (TLIM > 0) {
  wr <- file.path(normalizePath(OUTDIR), "inla_timeout.sh")
  writeLines(c("#!/bin/sh", sprintf('exec timeout %d "%s" "$@"', TLIM, inla.getOption("inla.call"))), wr); Sys.chmod(wr, "755")
  inla.setOption(inla.call = wr)
}

# ---- 2. DATA + NEIGHBOURHOOD GRAPH ---------------------------------------------------------------
g <- st_read(GPKG, quiet = TRUE)
fld <- g[!duplicated(g$ID.Field), c("ID.Field")]; fld <- fld[order(fld$ID.Field), ]
stopifnot(all(fld$ID.Field == seq_len(nrow(fld))))
nb <- st_is_within_distance(fld, fld, dist = DIST)
W <- sparseMatrix(i = rep(seq_along(nb), lengths(nb)), j = unlist(nb), x = 1, dims = rep(nrow(fld), 2)); diag(W) <- 0; W <- drop0(W)
cat("fields:", nrow(W), " mean neighbours:", round(mean(rowSums(W)), 2), " isolated:", sum(rowSums(W) == 0), "\n")
d <- st_drop_geometry(g); d <- d[order(d$ID.Field, d$Year), ]
d$ID.Year <- d$Year - min(d$Year) + 1L
d$Bloc <- factor(d$Bloc); d$irrg_ty <- factor(d$irrg_ty)
# derived responses: logs for right-skewed variables, ILR for sand / silt / clay (closed to 100 %)
d$log_EC <- log(d$Electrical_Conductivity_mmho_cm); d$log_P2O5 <- log(d$Available_P2O5_mg_kg)
d$log_K <- log(d$Exchangeable_Potassium_meq_100g);  d$log_Na <- log(d$Exchangeable_Sodium_meq_100g)
d$Ca <- d$Exchangeable_Calcium_meq_100g; d$Mg <- d$Exchangeable_Magnesium_meq_100g
d$OC <- d$Organic_Carbon_per_mille; d$N <- d$Total_Nitrogen_per_mille
# pivot ILR coordinates of sand / silt / clay (same definition as in the .qmd)
sa <- pmax(d$Total_Sand_percent, 0.1); si <- pmax(d$Total_Silt_percent, 0.1); cl <- pmax(d$Clay_percent, 0.1)
d$ILR3_1 <- sqrt(2 / 3) * log(sa / sqrt(si * cl))    # sand | silt, clay
d$ILR3_2 <- sqrt(1 / 2) * log(si / cl)               # silt | clay
d$ID.Field2 <- d$ID.Field; d$ID.Field3 <- d$ID.Field; d$ID.Year2 <- d$ID.Year
d$ID.Obs <- seq_len(nrow(d)); d$YearC <- as.numeric(scale(d$ID.Year))

# ---- 3. PRIORS (PC priors; the response is standardised, so the same prior is used for every variable) ---
pcp   <- function(u = 1, a = 0.01) list(prec = list(prior = "pc.prec", param = c(u, a)))
rho_  <- list(rho = list(prior = "pc.cor1", param = c(0.9, 0.9)))
bym2h <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), phi = list(prior = "pc", param = c(0.5, 0.5)))

# ---- 4. CANDIDATE STRUCTURES ---------------------------------------------------------------------
bes  <- 'f(ID.Field, model = "besag", graph = W, scale.model = TRUE, constr = TRUE, hyper = pcp())'
ar1t <- 'f(ID.Year, model = "ar1", hyper = c(pcp(), rho_))'
kron <- function(idx, gm) paste0('f(', idx, ', model = "besag", graph = W, scale.model = TRUE, constr = TRUE, group = ID.Year, control.group = list(model = "', gm, '"',
                                 if (gm == "ar1") ", hyper = rho_" else "", '), hyper = pcp())')
STRUCT <- list(
  S00_none                 = "",
  S01_iidField             = 'f(ID.Field, model = "iid", hyper = pcp())',
  S02_besag                = bes,
  S03_bym2                 = 'f(ID.Field, model = "bym2", graph = W, scale.model = TRUE, constr = TRUE, hyper = bym2h)',
  S04_ar1Time              = ar1t,
  S05_rw1Time              = 'f(ID.Year, model = "rw1", scale.model = TRUE, constr = TRUE, hyper = pcp())',
  S06_rw2Time              = 'f(ID.Year, model = "rw2", scale.model = TRUE, constr = TRUE, hyper = pcp())',
  S07_iidField_ar1Time     = paste('f(ID.Field, model = "iid", hyper = pcp())', ar1t, sep = " + "),
  S08_besag_ar1Time        = paste(bes, ar1t, sep = " + "),
  S09_bym2_ar1Time         = paste('f(ID.Field, model = "bym2", graph = W, scale.model = TRUE, constr = TRUE, hyper = bym2h)', ar1t, sep = " + "),
  S10_KH1_iidObs           = paste(bes, ar1t, 'f(ID.Obs, model = "iid", hyper = pcp())', sep = " + "),
  S11_KH2_iidSp_x_ar1      = paste(bes, ar1t, 'f(ID.Field2, model = "iid", group = ID.Year, control.group = list(model = "ar1", hyper = rho_), hyper = pcp())', sep = " + "),
  S12_KH3_besag_per_year   = paste(bes, ar1t, 'f(ID.Field2, model = "besag", graph = W, scale.model = TRUE, constr = TRUE, replicate = ID.Year, hyper = pcp())', sep = " + "),
  S13_KH4_besag_x_ar1      = paste(bes, ar1t, kron("ID.Field2", "ar1"), sep = " + "),
  S14_CURRENT_besag_x_ar1  = kron("ID.Field", "ar1"),                                  # the model used so far
  S15_besag_x_rw1          = kron("ID.Field", "rw1"),
  S16_besag_x_exch         = kron("ID.Field", "exchangeable"),
  S17_iidField_besag_x_ar1 = paste('f(ID.Field2, model = "iid", hyper = pcp())', kron("ID.Field", "ar1"), sep = " + "),
  S18_ar1PerField          = 'f(ID.Year, model = "ar1", replicate = ID.Field, hyper = c(pcp(), rho_))',
  S19_iid_ar1PerField      = 'f(ID.Field2, model = "iid", hyper = pcp()) + f(ID.Year, model = "ar1", replicate = ID.Field, hyper = c(pcp(), rho_))',
  S20_besag_ar1PerField    = paste(bes, 'f(ID.Year, model = "ar1", replicate = ID.Field3, hyper = c(pcp(), rho_))', sep = " + "),
  S21_besag_slope          = paste(bes, 'f(ID.Field2, YearC, model = "iid", hyper = pcp())', sep = " + "),
  S22_besag_ar1_slope      = paste(bes, ar1t, 'f(ID.Field2, YearC, model = "iid", hyper = pcp())', sep = " + "),
  S23_bym2_ar1_iidObs      = paste('f(ID.Field, model = "bym2", graph = W, scale.model = TRUE, constr = TRUE, hyper = bym2h)', ar1t, 'f(ID.Obs, model = "iid", hyper = pcp())', sep = " + "),
  S24_bym2_besag_x_ar1     = paste('f(ID.Field2, model = "bym2", graph = W, scale.model = TRUE, constr = TRUE, hyper = bym2h)', kron("ID.Field", "ar1"), sep = " + ")
)
if (MODELS != "all") STRUCT <- STRUCT[grepl(MODELS, names(STRUCT))]

# ---- 5. HELPERS ----------------------------------------------------------------------------------
make_mask <- function(y, field, scheme, seed) {          # TRUE = held out
  set.seed(seed); obs <- which(!is.na(y)); m <- rep(FALSE, length(y))
  if (scheme == "row") m[sample(obs, ceiling(HOLD * length(obs)))] <- TRUE
  else { f <- unique(field[obs]); h <- sample(f, ceiling(HOLD * length(f))); m[field %in% h & !is.na(y)] <- TRUE }
  m
}
fit_one <- function(v, st, fs, held, theta = NULL) {
  dd <- d; y <- dd[[v]]; y[!is.finite(y)] <- NA; dd$y_raw <- y
  dd$Y <- ifelse(held, NA, y); mu <- mean(dd$Y, na.rm = TRUE); sdv <- sd(dd$Y, na.rm = TRUE); dd$Y <- (dd$Y - mu) / sdv
  rhs <- paste(c(fs, if (nzchar(STRUCT[[st]])) STRUCT[[st]]), collapse = " + ")
  fml <- as.formula(paste("Y ~", rhs))
  t0 <- Sys.time()
  fit <- tryCatch(inla(fml, data = dd, family = "gaussian", control.family = list(hyper = pcp()),
                       control.predictor = list(compute = TRUE, link = 1),
                       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                       control.inla = list(strategy = "gaussian", int.strategy = STRAT),       # Gaussian likelihood: the Gaussian approximation is exact
                       control.mode = if (is.null(theta)) list() else list(theta = theta, restart = TRUE), verbose = FALSE),
                  error = function(e) { cat("   ! ", st, ": ", substr(conditionMessage(e), 1, 150), "\n"); NULL })
  if (is.null(fit) || is.null(fit$summary.fitted.values)) return(NULL)
  list(fit = fit, d = dd, mu = mu, sdv = sdv, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
app <- function(x, file) write.table(x, file, sep = ",", row.names = FALSE, col.names = !file.exists(file), append = file.exists(file), qmethod = "double")

# ---- 6. MAIN LOOP (resumable: finished variable x fixed x structure combinations are skipped) ----
RF <- file.path(OUTDIR, "results.csv"); PF <- file.path(OUTDIR, "holdout_pred.csv"); HF <- file.path(OUTDIR, "hyper.csv")
done <- if (file.exists(RF)) with(read.csv(RF), paste(Variable, Fixed, Model, sep = "|")) else character(0)
for (v in VARS) for (fs in FIXED) for (st in names(STRUCT)) {
  key <- paste(v, fs, st, sep = "|"); if (key %in% done) next
  cat(format(Sys.time(), "%H:%M:%S"), key, "\n")
  full <- fit_one(v, st, fs, rep(FALSE, nrow(d)))
  if (is.null(full)) { app(data.frame(Variable = v, Fixed = fs, Model = st, ok = 0), file.path(OUTDIR, "failed.csv")); next }
  f <- full$fit; cpo <- f$cpo$cpo; okc <- !is.na(cpo) & cpo > 0 & (is.na(f$cpo$failure) | f$cpo$failure < 1)
  out <- data.frame(Variable = v, Fixed = fs, Model = st, n_obs = sum(!is.na(full$d$Y)), minutes = full$minutes,
                    DIC = f$dic$dic, pD = f$dic$p.eff, WAIC = f$waic$waic, pW = f$waic$p.eff,
                    LCPO = -mean(log(cpo[okc])), cpo_fail = sum(!okc & !is.na(full$d$Y)), MLIK = f$mlik[1, 1])
  h <- f$summary.hyperpar; app(data.frame(Variable = v, Fixed = fs, Model = st, par = rownames(h), mean = h$mean, sd = h$sd, q025 = h$`0.025quant`, q975 = h$`0.975quant`), HF)
  for (sc in SCHEMES) {
    S <- NULL
    for (rp in seq_len(N_REP)) {
      held <- make_mask(full$d$y_raw, d$ID.Field, sc, SEED + rp); r <- fit_one(v, st, fs, held, theta = f$mode$theta); if (is.null(r)) next
      i <- which(held); s <- r$fit$summary.fitted.values
      hh <- r$fit$summary.hyperpar; s2e <- 1 / hh[grep("Gaussian observations", rownames(hh)), "mean"]      # observation-noise variance
      p <- s$mean[i]; ps <- sqrt(s$sd[i]^2 + s2e); lo <- p - 1.96 * ps; hi <- p + 1.96 * ps                  # predictive distribution of a new measurement
      y <- (r$d$y_raw[i] - r$mu) / r$sdv; z <- (y - p) / ps                                                 # standardised scale (training mean / sd)
      S <- rbind(S, data.frame(RMSE = sqrt(mean((p - y)^2)), MAE = mean(abs(p - y)), R2 = 1 - sum((p - y)^2) / sum((y - mean(y))^2),
                               COV95 = mean(y >= lo & y <= hi), WIDTH = mean(hi - lo),
                               CRPS = mean(ps * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))),
                               LOGS = -mean(dnorm(y, p, ps, log = TRUE)), minutes = r$minutes))
      app(data.frame(Variable = v, Fixed = fs, Model = st, scheme = sc, rep = rp, row = i, field = d$ID.Field[i], year = d$Year[i],
                     y = y, pred = p, lo = lo, hi = hi, y_raw = r$d$y_raw[i], mu = r$mu, sd = r$sdv), PF)
    }
    if (!is.null(S)) { m <- colMeans(S); names(m) <- paste0(names(m), "_", sc); out <- cbind(out, as.data.frame(t(m))) }
  }
  app(out, RF)
}
cat("DONE\n")
