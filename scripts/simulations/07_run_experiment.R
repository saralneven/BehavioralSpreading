library(here)

# 07_run_experiment.R

library(dplyr)
library(parallel)
library(data.table)

# Ensure parallel-safe RNG streams
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# Source all modules
src_dir <- here::here("scripts", "simulations")
source(file.path(src_dir, "00_config.R"))
source(file.path(src_dir, "01_io_and_helpers.R"))
source(file.path(src_dir, "02_load_model_fits.R"))
source(file.path(src_dir, "03_fr_speed.R"))
source(file.path(src_dir, "04_gmm_loader.R"))
source(file.path(src_dir, "04_empirical_loader.R"))
source(file.path(src_dir, "05_geometry.R"))
source(file.path(src_dir, "06_simulation_core.R"))

# ------------------------------------------------
# If binary FR-speed: derive weak/strong from DATA
# ------------------------------------------------
if (EXP$FR_speed_mode == "binary") {
  ws <- get_weak_strong_FR_speeds(
    path  = FR_information,
    probs = c(0.25, 0.75)
  )
  print(ws)

  EXP$FR_speed_binary_values <- ws$FR_speed_m_per_s
  EXP$FR_speed_labels        <- ws$label
}

# ------------------------------------------------
# Sanity: if grid mode, require FR_speed_grid_values
# ------------------------------------------------
if (EXP$FR_speed_mode == "grid") {
  if (is.null(EXP$FR_speed_grid_values)) {
    stop("EXP$FR_speed_mode == 'grid' but EXP$FR_speed_grid_values is NULL. ",
         "Add it in 00_config.R, e.g. seq(0.4, 3.0, by=0.1).")
  }
  EXP$FR_speed_grid_values <- as.numeric(EXP$FR_speed_grid_values)
  EXP$FR_speed_grid_values <- EXP$FR_speed_grid_values[is.finite(EXP$FR_speed_grid_values)]
  if (length(EXP$FR_speed_grid_values) == 0) {
    stop("EXP$FR_speed_grid_values is empty after filtering non-finite values.")
  }
}

# Build run tag + output folder
run_tag   <- make_run_tag(EXP)
safe_path <- make_safe_folder(out_base, run_tag)

log_msg("Run tag: %s", run_tag)
log_msg("Output folder: %s", safe_path)
log_msg("Network source: %s", EXP$network_source)

# Load fitted model parameters
mf <- load_model_fits(fits_csv, EXP)
get_model_spec <- mf$get_model_spec

# Keep only models available in fitted params
available_models <- unique(mf$params_by_model$model)
if (!all(EXP$models_to_run %in% available_models)) {
  missing <- setdiff(EXP$models_to_run, available_models)
  warning(sprintf("Skipping models not found in CSV: %s", paste(missing, collapse = ", ")))
  EXP$models_to_run <- intersect(EXP$models_to_run, available_models)
}
stopifnot(length(EXP$models_to_run) > 0)

# Fit FR speed distribution (used in lognormal mode; also supplies mu/sigma but ignored in binary/fixed/grid)
speed_dist <- fit_fr_speed_lognormal(FR_information)

# ------------------------------------------------
# FR speed conditions to run
# - lognormal: one condition (no label/value)
# - fixed: one condition
# - binary: labels in EXP$FR_speed_labels
# - grid: numeric values in EXP$FR_speed_grid_values
# ------------------------------------------------
fr_mode <- EXP$FR_speed_mode

if (fr_mode == "grid") {
  fr_cond_type   <- "value"
  fr_cond_values <- EXP$FR_speed_grid_values
} else if (fr_mode == "binary") {
  fr_cond_type   <- "label"
  fr_cond_values <- EXP$FR_speed_labels
} else {
  fr_cond_type   <- "none"
  fr_cond_values <- ""  # single iteration
}

# Ensure N_FRs iterable
N_FRs_vec <- as.integer(EXP$N_FRs)

# Prepare network samplers depending on source
get_gmm_for_size <- NULL
sample_from_gmm_3d <- NULL
sample_real_positions_3d <- NULL
available_real_sizes <- NULL

if (EXP$network_source == "gmm") {
  gmm_obj <- load_gmms_by_size(gmm_csv)
  get_gmm_for_size   <- gmm_obj$get_gmm_for_size
  sample_from_gmm_3d <- gmm_obj$sample_from_gmm_3d

} else if (EXP$network_source == "empirical") {
  emp <- load_empirical_positions_3d(real_positions_csv, EXP)
  sample_real_positions_3d <- emp$sample_real_positions_3d
  available_real_sizes <- emp$available_real_sizes

  log_msg("Empirical available sizes: %s", paste(available_real_sizes, collapse = ", "))

} else {
  stop("Unknown EXP$network_source: ", EXP$network_source)
}

# ------------------------------------------------
# Helper: construct a per-condition filename WITHOUT overwriting
# (We don't rely on make_partial_filename for grid mode.)
# ------------------------------------------------
make_partial_file_safe <- function(folder, model_name, A, N, nfr, EXP, run_tag,
                                  fr_speed_label = NA_character_,
                                  fr_speed_value = NA_real_) {
  parts <- c(
    paste0("Model", model_name),
    paste0("A", fmt_num(A, 1)),
    paste0("N", N),
    paste0("NFR", nfr),
    paste0("FRspd", EXP$FR_speed_mode)
  )

  if (EXP$FR_speed_mode == "fixed") {
    parts <- c(parts, paste0("FRv", fmt_num(EXP$FR_speed_fixed_value, 3)))
  }
  if (EXP$FR_speed_mode == "binary") {
    parts <- c(parts, paste0("FR", fr_speed_label))
  }
  if (EXP$FR_speed_mode == "grid") {
    parts <- c(parts, paste0("FRv", fmt_num(fr_speed_value, 3)))
  }

  parts <- c(
    parts,
    paste0("dec", fmt_num(EXP$speed_decay, 4)),
    run_tag
  )

  file.path(folder, paste0(paste(parts, collapse = "_"), ".csv"))
}

all_results <- list()
idx <- 1L

for (A in EXP$areas) {
  log_msg("Area label: %.1f", A)

  for (N in EXP$group_sizes) {
    log_msg("  Group size: %d", N)

    for (m in EXP$models_to_run) {

      for (nfr in N_FRs_vec) {

        for (cond in fr_cond_values) {

          use_lab   <- NA_character_
          use_value <- NA_real_

          if (fr_cond_type == "label") {
            use_lab <- as.character(cond)
          } else if (fr_cond_type == "value") {
            use_value <- as.numeric(cond)
          }

          if (EXP$network_source == "gmm") {
            sim_results <- simulate_condition_on_geometry(
              model_name = m,
              Area = A,
              group_size = N,
              N_FRs = nfr,
              num_reps = EXP$num_reps,
              speed_dist = speed_dist,
              EXP = EXP,
              get_model_spec = get_model_spec,
              get_gmm_for_size = get_gmm_for_size,
              sample_from_gmm_3d = sample_from_gmm_3d,
              mc_cores = mc_cores,
              fr_speed_label = use_lab,
              fr_speed_value = use_value
            )
          } else {
            sim_results <- simulate_condition_on_empirical_geometry(
              model_name = m,
              Area = A,
              group_size = N,   # will map to nearest available size if EXP$empirical_use_nearest_size=TRUE
              N_FRs = nfr,
              num_reps = EXP$num_reps,
              speed_dist = speed_dist,
              EXP = EXP,
              get_model_spec = get_model_spec,
              sample_real_positions_3d = sample_real_positions_3d,
              mc_cores = mc_cores,
              fr_speed_label = use_lab,
              fr_speed_value = use_value
            )
          }

          partial_file <- make_partial_file_safe(
            folder = safe_path,
            model_name = m,
            A = A,
            N = N,
            nfr = nfr,
            EXP = EXP,
            run_tag = run_tag,
            fr_speed_label = use_lab,
            fr_speed_value = use_value
          )

          data.table::fwrite(sim_results, partial_file)
          log_msg("    Saved: %s (rows=%d)", partial_file, nrow(sim_results))

          all_results[[idx]] <- sim_results
          idx <- idx + 1L
        }
      }
    }

    gc(FALSE)
  }
}

results_df <- dplyr::bind_rows(all_results)

final_file <- file.path(
  out_base,
  if (EXP$network_source == "gmm") {
    sprintf("Experimental_Simulations_3D_GMM_Results_%s.csv", run_tag)
  } else {
    sprintf("Experimental_Simulations_3D_REAL_Results_%s.csv", run_tag)
  }
)
data.table::fwrite(results_df, final_file)
log_msg("Saved final %s (rows=%d)", final_file, nrow(results_df))

err_df <- subset(results_df, !is.na(Error))
if (nrow(err_df) > 0) {
  log_msg("WARNING: %d replicate(s) returned errors. Showing first 10:", nrow(err_df))
  print(utils::head(err_df[, c("Model","Density","Video_ID","Group_Size_Target","Group_Size_Used","Error")], 10))
}
