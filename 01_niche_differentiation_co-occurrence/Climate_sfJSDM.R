#!/usr/bin/env Rscript
# Analyze climatic responses, PCA-env niche overlap, and co-occurrence coverage
# Author: Ryosuke Ito
# Qmon-Qser climatic segregation with spOccupancy::sfJSDM(), followed by
# availability-corrected PCA-env niche overlap and co-occurrence coverage.
#
# Standalone script: no other local R scripts are required.
#
# Smoke test (pipeline check only):
#   QMON_MODE=smoke Rscript Climate_sfJSDM.R Quercus_climate.csv
#
# Full analysis:
#   QMON_MODE=full Rscript Climate_sfJSDM.R Quercus_climate.csv

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260818)

required_packages <- c(
  "spOccupancy", "sf", "readr", "dplyr", "ggplot2", "coda", "RANN", "tibble",
  "ade4", "ecospat", "terra", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install missing packages first:\ninstall.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args)) args[1] else Sys.getenv(
  "QMON_DATA",
  "Quercus_climate.csv"
)
mode <- tolower(Sys.getenv("QMON_MODE", "smoke"))
if (!mode %in% c("smoke", "full")) {
  stop("QMON_MODE must be smoke or full.", call. = FALSE)
}
if (!file.exists(input_file)) stop("Input file not found: ", input_file, call. = FALSE)

out_dir <- file.path("spoccupancy_results", paste0("bradfer_", mode))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

climate_vars <- c("AMT", "APCP", "SD", "RH")
species_names <- c("Qmon", "Qser")
species_colours <- c(Qmon = "#313D32", Qser = "#D2C0A7")
vif_limit <- 2.5

# -----------------------------------------------------------------------------
# 1. Read, validate and document the analysis data
# -----------------------------------------------------------------------------
dat0 <- readr::read_csv(
  input_file,
  col_types = readr::cols(plot_id = readr::col_character()),
  show_col_types = FALSE,
  na = c("", "NA", "NaN", "NULL")
)

needed <- c("plot_id", species_names, "longitude", "latitude", climate_vars)
missing_columns <- setdiff(needed, names(dat0))
if (length(missing_columns)) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

dat0 <- dat0 |>
  dplyr::mutate(
    Qmon = as.integer(.data$Qmon),
    Qser = as.integer(.data$Qser),
    longitude = as.numeric(.data$longitude),
    latitude = as.numeric(.data$latitude),
    dplyr::across(dplyr::all_of(climate_vars), as.numeric)
  )

if (any(!is.na(dat0$Qmon) & !dat0$Qmon %in% 0:1) ||
    any(!is.na(dat0$Qser) & !dat0$Qser %in% 0:1)) {
  stop("Qmon and Qser must contain only 0, 1, or NA.", call. = FALSE)
}
if (anyDuplicated(dat0$plot_id)) stop("plot_id is duplicated.", call. = FALSE)

complete <- stats::complete.cases(dat0[, needed])
audit <- data.frame(
  input_rows = nrow(dat0),
  complete_rows = sum(complete),
  excluded_incomplete_rows = sum(!complete),
  Qmon_present = sum(dat0$Qmon[complete] == 1L),
  Qser_present = sum(dat0$Qser[complete] == 1L),
  observed_cooccurrence = sum(dat0$Qmon[complete] == 1L & dat0$Qser[complete] == 1L),
  SSD_included = FALSE,
  SSD_exclusion_reason = "A priori exclusion: strong national correlation with SD"
)
readr::write_csv(audit, file.path(out_dir, "data_audit.csv"))
dat <- dat0[complete, , drop = FALSE]

# VIF is calculated for the untransformed focal predictors, as in the
# predictor-screening stage of Bradfer-Lawrence et al. (2025).
calculate_vif <- function(x) {
  out <- vapply(names(x), function(response) {
    others <- setdiff(names(x), response)
    fit_vif <- stats::lm(stats::reformulate(others, response), data = x)
    r2 <- summary(fit_vif)$r.squared
    if (!is.finite(r2) || r2 >= 1) Inf else 1 / (1 - r2)
  }, numeric(1))
  data.frame(variable = names(out), VIF = unname(out), row.names = NULL)
}

vif_table <- calculate_vif(dat[, climate_vars, drop = FALSE]) |>
  dplyr::arrange(dplyr::desc(.data$VIF))
vif_table$threshold <- vif_limit
vif_table$passes <- vif_table$VIF < vif_limit
readr::write_csv(vif_table, file.path(out_dir, "vif_focal_predictors.csv"))
if (any(!vif_table$passes)) {
  offenders <- vif_table$variable[!vif_table$passes]
  stop(
    "Bradfer-Lawrence screening criterion not met: VIF >= ", vif_limit,
    " for ", paste(offenders, collapse = ", "),
    ". See ", file.path(out_dir, "vif_focal_predictors.csv"),
    ". Do not remove a variable automatically; revise the biological model first.",
    call. = FALSE
  )
}

cor_matrix <- stats::cor(dat[, climate_vars], method = "pearson")
readr::write_csv(
  tibble::rownames_to_column(as.data.frame(cor_matrix), "variable"),
  file.path(out_dir, "climate_correlations.csv")
)

# Centre and scale using all complete forest plots, not only presence sites.
scaling <- data.frame(
  variable = climate_vars,
  mean = vapply(dat[climate_vars], mean, numeric(1)),
  sd = vapply(dat[climate_vars], stats::sd, numeric(1)),
  row.names = NULL
)
if (any(!is.finite(scaling$sd) | scaling$sd <= 0)) {
  stop("A focal predictor has zero or invalid standard deviation.", call. = FALSE)
}
readr::write_csv(scaling, file.path(out_dir, "climate_scaling.csv"))
for (v in climate_vars) {
  dat[[paste0(v, "_z")]] <-
    (dat[[v]] - scaling$mean[scaling$variable == v]) /
    scaling$sd[scaling$variable == v]
}

# -----------------------------------------------------------------------------
# 2. Project coordinates and resolve exact ties caused by rounded source data
# -----------------------------------------------------------------------------
japan_lcc <- paste(
  "+proj=lcc +lat_1=30 +lat_2=46 +lat_0=38 +lon_0=138",
  "+x_0=0 +y_0=0 +ellps=GRS80 +units=km +no_defs"
)
points_ll <- sf::st_as_sf(
  dat, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
)
points_lcc <- sf::st_transform(points_ll, crs = japan_lcc)
coords <- sf::st_coordinates(points_lcc)
colnames(coords) <- c("x_km", "y_km")

coord_key <- paste(round(coords[, 1], 6), round(coords[, 2], 6), sep = "_")
tied <- duplicated(coord_key) | duplicated(coord_key, fromLast = TRUE)
if (any(tied)) {
  tied_groups <- split(which(tied), coord_key[tied])
  for (idx in tied_groups) {
    angle <- seq(0, 2 * pi, length.out = length(idx) + 1L)[-(length(idx) + 1L)]
    coords[idx, 1] <- coords[idx, 1] + 0.001 * cos(angle)
    coords[idx, 2] <- coords[idx, 2] + 0.001 * sin(angle)
  }
  warning(sum(tied), " rows had tied rounded coordinates and were separated by 1 m.")
}

nn <- RANN::nn2(coords, coords, k = 2)
nearest_dist <- nn$nn.dists[, 2]
positive_dist <- nearest_dist[is.finite(nearest_dist) & nearest_dist > 0]
if (!length(positive_dist)) stop("Could not determine positive nearest-neighbour distances.")
min_dist <- min(positive_dist)
median_dist <- stats::median(positive_dist)
max_dist <- sqrt(diff(range(coords[, 1]))^2 + diff(range(coords[, 2]))^2)
readr::write_csv(
  data.frame(
    duplicate_coordinate_rows_jittered = sum(tied),
    jitter_radius_km = ifelse(any(tied), 0.001, 0),
    min_nearest_km = min_dist,
    median_nearest_km = median_dist,
    bounding_box_diagonal_km = max_dist
  ),
  file.path(out_dir, "spatial_distance_summary.csv")
)

# -----------------------------------------------------------------------------
# 3. Spatial JSDM
# -----------------------------------------------------------------------------
y <- rbind(Qmon = dat$Qmon, Qser = dat$Qser)
covs <- as.data.frame(dat[, paste0(climate_vars, "_z"), drop = FALSE])
jsdm_data <- list(y = y, covs = covs, coords = coords)

# Quadratic terms allow climatic optima. The linear coefficient is the local
# slope at the nationwide mean of that predictor (z = 0).
jsdm_formula <- ~ AMT_z + I(AMT_z^2) +
  APCP_z + I(APCP_z^2) +
  SD_z + I(SD_z^2) +
  RH_z + I(RH_z^2)

n_factors <- 1L
lambda_inits <- matrix(c(1, stats::rnorm(1)), nrow = 2, ncol = 1,
                       dimnames = list(species_names, NULL))
inits <- list(
  beta.comm = 0,
  beta = 0,
  tau.sq.beta = 1,
  lambda = lambda_inits,
  phi = 3 / median_dist
)
priors <- list(
  beta.comm.normal = list(mean = 0, var = 2.72),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  phi.unif = list(3 / max_dist, 3 / min_dist)
)
tuning <- list(phi = 1)

settings <- if (mode == "full") {
  # 100,000 iterations/chain - 50,000 burn-in, then thin by 10 =
  # 5,000 retained/chain and 20,000 retained across four chains.
  list(n.batch = 4000L, batch.length = 25L, n.burn = 50000L,
       n.thin = 10L, n.chains = 4L, n.neighbors = 15L)
} else {
  list(n.batch = 50L, batch.length = 20L, n.burn = 500L,
       n.thin = 5L, n.chains = 2L, n.neighbors = 10L)
}
settings$n.iter_per_chain <- settings$n.batch * settings$batch.length
settings$n.retained_per_chain <-
  (settings$n.iter_per_chain - settings$n.burn) / settings$n.thin
settings$n.retained_total <- settings$n.retained_per_chain * settings$n.chains
readr::write_csv(as.data.frame(settings), file.path(out_dir, "mcmc_settings.csv"))

fit_path <- file.path(out_dir, "sfJSDM_fit.rds")
reuse_fit <- identical(Sys.getenv("QMON_REUSE_FIT", "0"), "1")
if (reuse_fit) {
  if (!file.exists(fit_path)) {
    stop("QMON_REUSE_FIT=1 but no saved model was found: ", fit_path, call. = FALSE)
  }
  message("Reusing saved model and regenerating summaries/figures: ", fit_path)
  fit <- readRDS(fit_path)
} else {
  message("Fitting Bradfer-style sfJSDM in '", mode, "' mode to ", ncol(y), " sites...")
  fit <- spOccupancy::sfJSDM(
    formula = jsdm_formula,
    data = jsdm_data,
    inits = inits,
    priors = priors,
    tuning = tuning,
    cov.model = "exponential",
    NNGP = TRUE,
    n.neighbors = settings$n.neighbors,
    search.type = "cb",
    n.factors = n_factors,
    n.batch = settings$n.batch,
    batch.length = settings$batch.length,
    accept.rate = 0.43,
    n.omp.threads = 1,
    verbose = TRUE,
    n.report = if (mode == "full") 100L else 10L,
    n.burn = settings$n.burn,
    n.thin = settings$n.thin,
    n.chains = settings$n.chains,
    # Avoid retaining multi-gigabyte site-level psi, w and likelihood arrays.
    # The requested coefficient violins and marginal curves require beta only.
    monitors = c("beta.comm", "tau.sq.beta", "beta", "lambda", "theta")
  )
  saveRDS(fit, fit_path, compress = FALSE)
}

capture.output(summary(fit), file = file.path(out_dir, "model_summary.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

# -----------------------------------------------------------------------------
# 4. Posterior summaries and diagnostics
# -----------------------------------------------------------------------------
beta_draws <- as.matrix(fit$beta.samples)
beta_summary <- data.frame(
  parameter = colnames(beta_draws),
  mean = colMeans(beta_draws),
  sd = apply(beta_draws, 2, stats::sd),
  q2.5 = apply(beta_draws, 2, stats::quantile, probs = 0.025),
  median = apply(beta_draws, 2, stats::median),
  q97.5 = apply(beta_draws, 2, stats::quantile, probs = 0.975),
  row.names = NULL
)
if (length(fit$rhat$beta) == nrow(beta_summary)) {
  beta_summary$Rhat <- as.numeric(fit$rhat$beta)
}
if (length(fit$ESS$beta) == nrow(beta_summary)) {
  beta_summary$ESS <- as.numeric(fit$ESS$beta)
}
beta_summary$credible_nonzero <- with(
  beta_summary, (q2.5 > 0 & q97.5 > 0) | (q2.5 < 0 & q97.5 < 0)
)
readr::write_csv(beta_summary, file.path(out_dir, "species_climate_coefficients.csv"))

safe_max <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}
safe_min_positive <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x) & x > 0]
  if (length(x)) min(x) else NA_real_
}
diagnostics <- data.frame(
  metric = c(
    "max_Rhat_species_beta", "min_ESS_species_beta",
    "max_Rhat_spatial_theta", "min_ESS_spatial_theta",
    "max_Rhat_factor_loading", "min_ESS_factor_loading",
    "max_Rhat_community_beta", "min_ESS_community_beta"
  ),
  value = c(
    safe_max(fit$rhat$beta), safe_min_positive(fit$ESS$beta),
    safe_max(fit$rhat$theta), safe_min_positive(fit$ESS$theta),
    safe_max(fit$rhat$lambda.lower.tri), safe_min_positive(fit$ESS$lambda),
    safe_max(fit$rhat$beta.comm), safe_min_positive(fit$ESS$beta.comm)
  )
)
diagnostics$bradfer_Rhat_pass <- ifelse(
  grepl("Rhat", diagnostics$metric), diagnostics$value < 1.1, NA
)
readr::write_csv(diagnostics, file.path(out_dir, "mcmc_diagnostics.csv"))

grDevices::pdf(file.path(out_dir, "traceplots.pdf"), width = 10, height = 7)
plot(fit$beta.samples, main = "Species-level climate coefficients")
plot(fit$theta.samples, main = "Spatial range parameter")
plot(fit$lambda.samples, main = "Spatial factor loadings")
grDevices::dev.off()

# Paired posterior contrasts directly test whether climatic responses differ.
contrast_rows <- list()
for (v in climate_vars) {
  for (term_type in c("linear", "quadratic")) {
    term <- if (term_type == "linear") paste0(v, "_z") else paste0("I(", v, "_z^2)")
    qmon_name <- paste0(term, "-Qmon")
    qser_name <- paste0(term, "-Qser")
    if (!all(c(qmon_name, qser_name) %in% colnames(beta_draws))) {
      stop("Could not find coefficient columns: ", qmon_name, ", ", qser_name)
    }
    delta <- beta_draws[, qmon_name] - beta_draws[, qser_name]
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      variable = v,
      term = term_type,
      contrast = "Qmon - Qser",
      mean = mean(delta),
      sd = stats::sd(delta),
      median = stats::median(delta),
      q2.5 = stats::quantile(delta, 0.025),
      q97.5 = stats::quantile(delta, 0.975),
      Pr_Qmon_gt_Qser = mean(delta > 0),
      credible_difference = stats::quantile(delta, 0.025) > 0 |
        stats::quantile(delta, 0.975) < 0
    )
  }
}
contrast_table <- do.call(rbind, contrast_rows)
readr::write_csv(contrast_table, file.path(out_dir, "Qmon_Qser_coefficient_contrasts.csv"))

# Portable posterior files preserve chains and saved iteration numbers. These
# files reproduce coefficient violins and convergence summaries without relying
# on the package-specific fitted-model object.
mcmc_draws_table <- function(samples) {
  if (inherits(samples, "mcmc.list")) {
    return(dplyr::bind_rows(lapply(seq_along(samples), function(chain) {
      values <- as.matrix(samples[[chain]])
      data.frame(
        chain = chain,
        iteration = as.numeric(stats::time(samples[[chain]])),
        values,
        check.names = FALSE
      )
    })))
  }
  if (inherits(samples, "mcmc")) {
    values <- as.matrix(samples)
    return(data.frame(
      chain = 1L,
      iteration = as.numeric(stats::time(samples)),
      values,
      check.names = FALSE
    ))
  }
  values <- as.matrix(samples)
  data.frame(
    chain = 1L,
    iteration = seq_len(nrow(values)),
    values,
    check.names = FALSE
  )
}

posterior_objects <- list(
  beta = fit$beta.samples,
  beta_community = fit$beta.comm.samples,
  beta_variance = fit$tau.sq.beta.samples,
  spatial_theta = fit$theta.samples,
  factor_loading = fit$lambda.samples
)
for (posterior_name in names(posterior_objects)) {
  samples <- posterior_objects[[posterior_name]]
  if (!is.null(samples)) {
    readr::write_csv(
      mcmc_draws_table(samples),
      file.path(out_dir, paste0("posterior_", posterior_name, "_draws.csv.gz"))
    )
  }
}

# -----------------------------------------------------------------------------
# 5. Split violin plots (left Qmon, right Qser)
# -----------------------------------------------------------------------------
half_violin_polygon <- function(values, x_center, side, width = 0.42, n = 256L) {
  values <- values[is.finite(values)]
  if (length(values) < 2L || stats::sd(values) == 0) {
    return(data.frame(x = numeric(), y = numeric()))
  }
  den <- stats::density(values, n = n, adjust = 1.1)
  scaled <- den$y / max(den$y) * width
  edge <- if (side == "left") x_center - scaled else x_center + scaled
  data.frame(
    x = c(x_center, edge, x_center),
    y = c(den$x[1], den$x, den$x[length(den$x)])
  )
}

# Main inferential figure: posterior distributions of standardized linear
# effects. With quadratic terms included, these are local slopes at z = 0.
posterior_poly <- list()
posterior_intervals <- list()
for (i in seq_along(climate_vars)) {
  v <- climate_vars[i]
  for (sp in species_names) {
    nm <- paste0(v, "_z-", sp)
    vals <- beta_draws[, nm]
    poly <- half_violin_polygon(
      vals, x_center = i, side = if (sp == "Qmon") "left" else "right"
    )
    poly$variable <- v
    poly$species <- sp
    posterior_poly[[length(posterior_poly) + 1L]] <- poly
    posterior_intervals[[length(posterior_intervals) + 1L]] <- data.frame(
      variable = v,
      species = sp,
      x = i + ifelse(sp == "Qmon", -0.055, 0.055),
      q2.5 = stats::quantile(vals, 0.025),
      median = stats::median(vals),
      q97.5 = stats::quantile(vals, 0.975)
    )
  }
}
posterior_poly <- do.call(rbind, posterior_poly)
posterior_intervals <- do.call(rbind, posterior_intervals)

p_beta <- ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey45") +
  ggplot2::geom_polygon(
    data = posterior_poly,
    ggplot2::aes(x = x, y = y, group = interaction(variable, species), fill = species),
    alpha = 0.85, colour = "grey20", linewidth = 0.25
  ) +
  ggplot2::geom_linerange(
    data = posterior_intervals,
    ggplot2::aes(x = x, ymin = q2.5, ymax = q97.5, colour = species),
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    data = posterior_intervals,
    ggplot2::aes(x = x, y = median, colour = species), size = 1.7
  ) +
  ggplot2::scale_x_continuous(breaks = seq_along(climate_vars), labels = climate_vars) +
  ggplot2::scale_fill_manual(values = species_colours) +
  ggplot2::scale_colour_manual(values = species_colours) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::labs(
    x = NULL,
    y = "Posterior linear effect on logit occurrence\n(per 1 SD; local slope at mean climate)",
    fill = NULL,
    colour = NULL
  )
ggplot2::ggsave(
  file.path(out_dir, "posterior_climate_effect_split_violin.pdf"),
  p_beta, width = 8, height = 5.4
)
ggplot2::ggsave(
  file.path(out_dir, "posterior_climate_effect_split_violin.png"),
  p_beta, width = 8, height = 5.4, dpi = 400
)

# Descriptive supplement: realized climates at observed presence sites.
# Co-occurrence sites correctly contribute once to each species distribution.
observed_poly <- list()
observed_medians <- list()
for (v in climate_vars) {
  for (sp in species_names) {
    vals <- dat[[v]][dat[[sp]] == 1L]
    poly <- half_violin_polygon(
      vals, x_center = 1, side = if (sp == "Qmon") "left" else "right"
    )
    poly$variable <- v
    poly$species <- sp
    observed_poly[[length(observed_poly) + 1L]] <- poly
    observed_medians[[length(observed_medians) + 1L]] <- data.frame(
      variable = v,
      species = sp,
      x = 1 + ifelse(sp == "Qmon", -0.04, 0.04),
      median = stats::median(vals)
    )
  }
}
observed_poly <- do.call(rbind, observed_poly)
observed_medians <- do.call(rbind, observed_medians)

p_observed <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = observed_poly,
    ggplot2::aes(x = x, y = y, group = interaction(variable, species), fill = species),
    alpha = 0.85, colour = "grey20", linewidth = 0.25
  ) +
  ggplot2::geom_point(
    data = observed_medians,
    ggplot2::aes(x = x, y = median, colour = species), size = 1.8
  ) +
  ggplot2::facet_wrap(~ variable, scales = "free_y", nrow = 1) +
  ggplot2::scale_x_continuous(limits = c(0.5, 1.5), breaks = NULL) +
  ggplot2::scale_fill_manual(values = species_colours) +
  ggplot2::scale_colour_manual(values = species_colours) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::labs(x = NULL, y = "Observed climate value at presence sites", fill = NULL, colour = NULL)
ggplot2::ggsave(
  file.path(out_dir, "observed_presence_climate_split_violin.pdf"),
  p_observed, width = 10, height = 4.5
)

# -----------------------------------------------------------------------------
# 6. Bradfer-style marginal response curves
# -----------------------------------------------------------------------------
response_rows <- list()
for (v in climate_vars) {
  raw_grid <- seq(
    stats::quantile(dat[[v]], 0.01),
    stats::quantile(dat[[v]], 0.99),
    length.out = 120L
  )
  mu <- scaling$mean[scaling$variable == v]
  sig <- scaling$sd[scaling$variable == v]
  z_grid <- (raw_grid - mu) / sig
  for (sp in species_names) {
    intercept <- beta_draws[, paste0("(Intercept)-", sp)]
    linear <- beta_draws[, paste0(v, "_z-", sp)]
    quadratic <- beta_draws[, paste0("I(", v, "_z^2)-", sp)]
    eta <- outer(intercept, rep(1, length(z_grid))) +
      outer(linear, z_grid) + outer(quadratic, z_grid^2)
    prob <- stats::plogis(eta)
    response_rows[[length(response_rows) + 1L]] <- data.frame(
      variable = v,
      species = sp,
      climate_value = raw_grid,
      mean = colMeans(prob),
      q2.5 = apply(prob, 2, stats::quantile, probs = 0.025),
      q97.5 = apply(prob, 2, stats::quantile, probs = 0.975)
    )
  }
}
response_table <- do.call(rbind, response_rows)
readr::write_csv(response_table, file.path(out_dir, "marginal_climate_response_curves.csv"))

p_response <- ggplot2::ggplot(
  response_table,
  ggplot2::aes(x = climate_value, y = mean, colour = species, fill = species)
) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = q2.5, ymax = q97.5),
                       alpha = 0.18, colour = NA) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::facet_wrap(~ variable, scales = "free_x", ncol = 2) +
  ggplot2::scale_colour_manual(values = species_colours) +
  ggplot2::scale_fill_manual(values = species_colours) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::labs(
    x = "Climate value",
    y = "Occurrence probability\n(other standardized predictors fixed at 0)",
    colour = NULL,
    fill = NULL
  )
ggplot2::ggsave(
  file.path(out_dir, "marginal_climate_response_curves.pdf"),
  p_response, width = 8.5, height = 7
)

# -----------------------------------------------------------------------------
# 7. Availability-corrected PCA-env niche overlap and co-occurrence coverage
# -----------------------------------------------------------------------------
pca_dir <- file.path(out_dir, "pca_env")
dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)

pca_permutations_default <- if (mode == "full") 100L else 10L
pca_permutations <- suppressWarnings(as.integer(Sys.getenv(
  "QMON_PCA_PERMUTATIONS", as.character(pca_permutations_default)
)))
if (!is.finite(pca_permutations) || pca_permutations < 2L) {
  stop("QMON_PCA_PERMUTATIONS must be an integer >= 2.", call. = FALSE)
}

pca_audit <- data.frame(
  input_rows = nrow(dat0),
  analysed_rows = nrow(dat),
  excluded_incomplete_rows = nrow(dat0) - nrow(dat),
  Qmon_presences = sum(dat$Qmon == 1L),
  Qser_presences = sum(dat$Qser == 1L),
  observed_both = sum(dat$Qmon == 1L & dat$Qser == 1L),
  permutations = pca_permutations
)
readr::write_csv(pca_audit, file.path(pca_dir, "data_audit.csv"))
readr::write_csv(vif_table, file.path(pca_dir, "vif_selected_variables.csv"))

pca_environment <- as.data.frame(dat[, climate_vars, drop = FALSE])
pca_fit <- ade4::dudi.pca(
  pca_environment, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2
)
pca_background <- as.data.frame(pca_fit$li[, 1:2, drop = FALSE])
names(pca_background) <- c("PC1", "PC2")
pca_Qmon <- pca_background[dat$Qmon == 1L, , drop = FALSE]
pca_Qser <- pca_background[dat$Qser == 1L, , drop = FALSE]
pca_both <- pca_background[dat$Qmon == 1L & dat$Qser == 1L, , drop = FALSE]

pca_percent_variance <- 100 * pca_fit$eig[1:2] / sum(pca_fit$eig)
pca_loadings <- data.frame(
  variable = rownames(pca_fit$co),
  PC1 = pca_fit$co[, 1],
  PC2 = pca_fit$co[, 2],
  row.names = NULL
)
readr::write_csv(pca_loadings, file.path(pca_dir, "pca_loadings.csv"))
readr::write_csv(
  data.frame(axis = c("PC1", "PC2"), percent_variance = pca_percent_variance),
  file.path(pca_dir, "pca_variance_explained.csv")
)

pca_grid_resolution <- 100L
pca_z_Qmon <- ecospat::ecospat.grid.clim.dyn(
  pca_background, pca_background, pca_Qmon, R = pca_grid_resolution
)
pca_z_Qser <- ecospat::ecospat.grid.clim.dyn(
  pca_background, pca_background, pca_Qser, R = pca_grid_resolution
)
pca_z_both <- ecospat::ecospat.grid.clim.dyn(
  pca_background, pca_background, pca_both, R = pca_grid_resolution
)

pca_overlap <- ecospat::ecospat.niche.overlap(pca_z_Qmon, pca_z_Qser, cor = TRUE)
pca_overlap_table <- data.frame(
  comparison = "Qmon_vs_Qser",
  Schoener_D = as.numeric(pca_overlap$D),
  Hellinger_I = as.numeric(pca_overlap$I),
  environment_availability_correction = TRUE
)
readr::write_csv(pca_overlap_table, file.path(pca_dir, "niche_overlap.csv"))

message("Running PCA-env equivalency tests with ", pca_permutations, " permutations...")
set.seed(20260818)
pca_equiv_higher <- ecospat::ecospat.niche.equivalency.test(
  pca_z_Qmon, pca_z_Qser, rep = pca_permutations, intersection = 0,
  overlap.alternative = "higher", ncores = 1
)
pca_equiv_lower <- ecospat::ecospat.niche.equivalency.test(
  pca_z_Qmon, pca_z_Qser, rep = pca_permutations, intersection = 0,
  overlap.alternative = "lower", ncores = 1
)
pca_extract_named <- function(x, name) {
  if (is.list(x)) as.numeric(x[[name]]) else as.numeric(x[name])
}
pca_test_row <- function(result, alternative) {
  data.frame(
    test = "niche_equivalency",
    alternative = alternative,
    D_observed = pca_extract_named(result$obs, "D"),
    I_observed = pca_extract_named(result$obs, "I"),
    p_D = as.numeric(result$p.D),
    p_I = as.numeric(result$p.I),
    permutations = pca_permutations
  )
}
pca_equiv_table <- dplyr::bind_rows(
  pca_test_row(pca_equiv_higher, "higher (conservatism)"),
  pca_test_row(pca_equiv_lower, "lower (divergence)")
)
readr::write_csv(pca_equiv_table, file.path(pca_dir, "niche_equivalency_tests.csv"))

pca_raster_density <- function(raster, group) {
  x <- terra::as.data.frame(raster, xy = TRUE, na.rm = FALSE)
  if (ncol(x) < 3L) stop("Unexpected PCA-env raster for ", group, call. = FALSE)
  names(x)[3] <- "density"
  x$density[!is.finite(x$density)] <- 0
  if (any(x$density < 0)) stop("Negative PCA-env density for ", group, call. = FALSE)
  total <- sum(x$density)
  if (!is.finite(total) || total <= 0) stop("Zero PCA-env density for ", group, call. = FALSE)
  x$density <- x$density / total
  x$group <- group
  x
}

# Formal grids are used for D/I and coverage statistics.
pca_formal_Qmon <- pca_raster_density(pca_z_Qmon$z.cor, "Qmon")
pca_formal_Qser <- pca_raster_density(pca_z_Qser$z.cor, "Qser")
pca_formal_both <- pca_raster_density(pca_z_both$z.cor, "Both")
pca_same_grid <-
  identical(pca_formal_Qmon$x, pca_formal_Qser$x) &&
  identical(pca_formal_Qmon$y, pca_formal_Qser$y) &&
  identical(pca_formal_Qmon$x, pca_formal_both$x) &&
  identical(pca_formal_Qmon$y, pca_formal_both$y)
if (!pca_same_grid) stop("PCA-env density grids are not aligned.", call. = FALSE)
readr::write_csv(
  dplyr::bind_rows(pca_formal_Qmon, pca_formal_Qser, pca_formal_both),
  file.path(pca_dir, "pca_env_density_grid.csv.gz")
)

pca_hdr_level <- function(density, probability = 0.95) {
  values <- sort(density[is.finite(density) & density > 0], decreasing = TRUE)
  if (!length(values)) return(NA_real_)
  values[which(cumsum(values) >= probability)[1L]]
}
pca_level_Qmon <- pca_hdr_level(pca_formal_Qmon$density)
pca_level_Qser <- pca_hdr_level(pca_formal_Qser$density)
pca_level_both <- pca_hdr_level(pca_formal_both$density)
pca_in_Qmon <- pca_formal_Qmon$density >= pca_level_Qmon & pca_formal_Qmon$density > 0
pca_in_Qser <- pca_formal_Qser$density >= pca_level_Qser & pca_formal_Qser$density > 0
pca_in_both <- pca_formal_both$density >= pca_level_both & pca_formal_both$density > 0
pca_in_shared <- pca_in_Qmon & pca_in_Qser
pca_in_intersection <- pca_in_shared & pca_in_both
if (!sum(pca_in_shared) || !sum(pca_in_both)) {
  stop("A PCA-env 95% highest-density region contains no cells.", call. = FALSE)
}

pca_shared_density <- pmin(pca_formal_Qmon$density, pca_formal_Qser$density)
pca_shared_density_total <- sum(pca_shared_density)
pca_coverage <- data.frame(
  metric = c(
    "shared_95_area_covered_by_cooccurrence_95",
    "cooccurrence_95_area_inside_shared_95",
    "shared_density_mass_inside_cooccurrence_95",
    "cooccurrence_density_mass_inside_shared_95"
  ),
  percent = c(
    100 * sum(pca_in_intersection) / sum(pca_in_shared),
    100 * sum(pca_in_intersection) / sum(pca_in_both),
    100 * sum(pca_shared_density[pca_in_both]) / pca_shared_density_total,
    100 * sum(pca_formal_both$density[pca_in_shared]) / sum(pca_formal_both$density)
  )
)
readr::write_csv(pca_coverage, file.path(pca_dir, "cooccurrence_niche_coverage.csv"))
readr::write_csv(
  data.frame(
    grid_cells = nrow(pca_formal_Qmon),
    Qmon_95_cells = sum(pca_in_Qmon),
    Qser_95_cells = sum(pca_in_Qser),
    shared_95_cells = sum(pca_in_shared),
    cooccurrence_95_cells = sum(pca_in_both),
    intersection_cells = sum(pca_in_intersection),
    Qmon_95_density_threshold = pca_level_Qmon,
    Qser_95_density_threshold = pca_level_Qser,
    cooccurrence_95_density_threshold = pca_level_both,
    shared_density_total = pca_shared_density_total
  ),
  file.path(pca_dir, "cooccurrence_niche_coverage_audit.csv")
)

# Plot-only grids receive a 10% margin so contours are not clipped.
pca_pc1_range <- range(pca_background$PC1)
pca_pc2_range <- range(pca_background$PC2)
pca_plot_extension <- c(
  -0.10 * diff(pca_pc1_range), 0.10 * diff(pca_pc1_range),
  -0.10 * diff(pca_pc2_range), 0.10 * diff(pca_pc2_range)
)
pca_z_Qmon_plot <- ecospat::ecospat.grid.clim.dyn(
  pca_background, pca_background, pca_Qmon,
  R = pca_grid_resolution, extend.extent = pca_plot_extension
)
pca_z_Qser_plot <- ecospat::ecospat.grid.clim.dyn(
  pca_background, pca_background, pca_Qser,
  R = pca_grid_resolution, extend.extent = pca_plot_extension
)
pca_z_both_plot <- ecospat::ecospat.grid.clim.dyn(
  pca_background, pca_background, pca_both,
  R = pca_grid_resolution, extend.extent = pca_plot_extension
)
pca_available_plot <- pca_raster_density(pca_z_Qmon_plot$Z, "Available climate")
pca_Qmon_plot <- pca_raster_density(pca_z_Qmon_plot$z.cor, "Qmon")
pca_Qser_plot <- pca_raster_density(pca_z_Qser_plot$z.cor, "Qser")
pca_both_plot <- pca_raster_density(pca_z_both_plot$z.cor, "Both")
pca_plot_grids <- dplyr::bind_rows(
  pca_available_plot, pca_Qmon_plot, pca_Qser_plot, pca_both_plot
)
readr::write_csv(
  pca_plot_grids,
  file.path(pca_dir, "pca_env_plot_density_grid.csv.gz")
)

pca_contour_spec <- dplyr::bind_rows(lapply(
  list(Qmon = pca_Qmon_plot, Qser = pca_Qser_plot, Both = pca_both_plot),
  function(grid) data.frame(
    group = unique(grid$group), probability = 0.95,
    density_level = pca_hdr_level(grid$density)
  )
))
readr::write_csv(pca_contour_spec, file.path(pca_dir, "contour_95_percent_levels.csv"))

pca_make_contours <- function(grid, level, group) {
  x_values <- sort(unique(grid$x))
  y_values <- sort(unique(grid$y))
  z_matrix <- matrix(0, nrow = length(x_values), ncol = length(y_values))
  z_matrix[cbind(match(grid$x, x_values), match(grid$y, y_values))] <- grid$density
  lines <- grDevices::contourLines(x_values, y_values, z_matrix, levels = level)
  if (!length(lines)) return(NULL)
  dplyr::bind_rows(lapply(seq_along(lines), function(index) {
    data.frame(
      PC1 = lines[[index]]$x, PC2 = lines[[index]]$y,
      group = group, path = paste(group, index, sep = "_")
    )
  }))
}
pca_contours <- dplyr::bind_rows(lapply(seq_len(nrow(pca_contour_spec)), function(i) {
  group <- pca_contour_spec$group[i]
  grid <- pca_plot_grids[pca_plot_grids$group == group, , drop = FALSE]
  pca_make_contours(grid, pca_contour_spec$density_level[i], group)
}))
readr::write_csv(pca_contours, file.path(pca_dir, "pca_env_95_contours.csv"))

pca_density_centroid <- function(grid) {
  data.frame(
    group = unique(grid$group),
    PC1 = stats::weighted.mean(grid$x, grid$density),
    PC2 = stats::weighted.mean(grid$y, grid$density)
  )
}
pca_centroids <- dplyr::bind_rows(lapply(
  list(Qmon = pca_Qmon_plot, Qser = pca_Qser_plot, Both = pca_both_plot),
  pca_density_centroid
))
readr::write_csv(pca_centroids, file.path(pca_dir, "pca_env_centroids.csv"))
pca_centroid_Qmon <- pca_centroids[pca_centroids$group == "Qmon", ]
pca_centroid_Qser <- pca_centroids[pca_centroids$group == "Qser", ]
pca_centroid_distance <- sqrt(
  (pca_centroid_Qmon$PC1 - pca_centroid_Qser$PC1)^2 +
    (pca_centroid_Qmon$PC2 - pca_centroid_Qser$PC2)^2
)
readr::write_csv(
  data.frame(
    comparison = "Qmon_vs_Qser",
    centroid_distance_PCA_units = pca_centroid_distance
  ),
  file.path(pca_dir, "pca_centroid_distance.csv")
)

pca_arrow_multiplier <- 0.25 * min(
  diff(range(pca_available_plot$x)), diff(range(pca_available_plot$y))
) / max(sqrt(pca_loadings$PC1^2 + pca_loadings$PC2^2))
pca_arrow_data <- pca_loadings |>
  dplyr::mutate(
    PC1_plot = .data$PC1 * pca_arrow_multiplier,
    PC2_plot = .data$PC2 * pca_arrow_multiplier
  )
pca_required_extent <- max(
  abs(pca_contours$PC1), abs(pca_contours$PC2),
  abs(1.15 * pca_arrow_data$PC1_plot),
  abs(1.15 * pca_arrow_data$PC2_plot), na.rm = TRUE
)
pca_axis_limit <- max(4, ceiling(pca_required_extent * 1.08 * 2) / 2)
pca_axis_limits <- c(-pca_axis_limit, pca_axis_limit)
readr::write_csv(
  data.frame(axis_min = -pca_axis_limit, axis_max = pca_axis_limit),
  file.path(pca_dir, "figure_axis_limits.csv")
)

pca_colours <- c(
  Qmon = species_colours[["Qmon"]],
  Qser = species_colours[["Qser"]],
  Both = "#889D35"
)
pca_figure <- ggplot2::ggplot() +
  ggplot2::geom_raster(
    data = pca_available_plot,
    ggplot2::aes(.data$x, .data$y, fill = .data$density), alpha = 0.72
  ) +
  ggplot2::scale_fill_gradient(
    low = "#FFFFFF", high = "#3F3F3F", trans = "sqrt",
    name = "Climatic availability"
  ) +
  ggplot2::geom_path(
    data = pca_contours,
    ggplot2::aes(.data$PC1, .data$PC2, colour = .data$group,
                 group = .data$path), linewidth = 1.15
  ) +
  ggplot2::geom_segment(
    data = pca_arrow_data,
    ggplot2::aes(x = 0, y = 0, xend = .data$PC1_plot, yend = .data$PC2_plot),
    inherit.aes = FALSE, colour = "grey10", linewidth = 0.48,
    arrow = grid::arrow(length = grid::unit(2.1, "mm"), type = "closed")
  ) +
  ggplot2::geom_label(
    data = pca_arrow_data,
    ggplot2::aes(
      x = 1.12 * .data$PC1_plot,
      y = 1.12 * .data$PC2_plot,
      label = .data$variable
    ),
    inherit.aes = FALSE, size = 2.8, linewidth = 0,
    fill = scales::alpha("white", 0.78)
  ) +
  ggplot2::geom_point(
    data = pca_centroids,
    ggplot2::aes(.data$PC1, .data$PC2, colour = .data$group),
    shape = 4, stroke = 1.2, size = 3.1
  ) +
  ggplot2::scale_colour_manual(values = pca_colours, name = "95% HDR") +
  ggplot2::annotate(
    "label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.15,
    label = sprintf(
      "Schoener's D = %.3f\nHellinger's I = %.3f\nBoth n = %d",
      pca_overlap_table$Schoener_D, pca_overlap_table$Hellinger_I,
      pca_audit$observed_both
    ),
    size = 3.0, linewidth = 0.2
  ) +
  ggplot2::coord_cartesian(
    xlim = pca_axis_limits, ylim = pca_axis_limits,
    expand = FALSE, clip = "off"
  ) +
  ggplot2::labs(
    x = sprintf("PC1 (%.1f%%)", pca_percent_variance[1]),
    y = sprintf("PC2 (%.1f%%)", pca_percent_variance[2])
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(8, 18, 8, 10)
  )
ggplot2::ggsave(
  file.path(pca_dir, "figure_pca_env_niche_overlap.pdf"),
  pca_figure, width = 8.2, height = 6.7
)
ggplot2::ggsave(
  file.path(pca_dir, "figure_pca_env_niche_overlap.png"),
  pca_figure, width = 8.2, height = 6.7, dpi = 600
)

pca_reported_results <- dplyr::bind_rows(
  data.frame(metric = "Schoener_D", value = pca_overlap_table$Schoener_D),
  data.frame(metric = "Hellinger_I", value = pca_overlap_table$Hellinger_I),
  data.frame(
    metric = "Qmon_Qser_centroid_distance_PCA_units",
    value = pca_centroid_distance
  ),
  data.frame(metric = pca_coverage$metric, value = pca_coverage$percent)
)
readr::write_csv(
  pca_reported_results,
  file.path(pca_dir, "pca_env_reported_results.csv")
)
capture.output(sessionInfo(), file = file.path(pca_dir, "sessionInfo.txt"))

message("Finished. Results: ", normalizePath(out_dir))
message("PCA-env results: ", normalizePath(pca_dir))
if (mode == "smoke") {
  message("Important: smoke-mode MCMC is only a pipeline check, not an inferential result.")
}

