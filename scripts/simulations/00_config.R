# 00_config.R

EXP <- list(
  # -------------------------
  # What to simulate
  # -------------------------
  models_to_run = c("Di", "DiT", "DiTDe", "DiTDeS"),  # Choose the model types to simulate
  areas         = 1,                                  # Area size (Area x Area in meters)
  group_sizes   = 5:40,                               # Group sizes to simulate
  num_reps      = 200,                                # Number of repetitions per condition (model x group size) 
  N_FRs         = c(1),                               # Number of initiators (scalar or vector)

  # -------------------------
  # Network source
  # -------------------------
  # "gmm" | "empirical"
  network_source = "gmm",                             # Source of geometric networks for simulations

  # Only used if network_source == "empirical"
  empirical_snapshot_strategy = "unique_fish_first_row",  # robust default
  empirical_use_nearest_size  = FALSE,   # if requested size not available, map to nearest available
  empirical_max_tries         = 25L,     # resample video if unusable
  empirical_min_fish          = 2L,      # must have >=2 fish to simulate

  # -------------------------
  # FR speed specification
  # -------------------------
  FR_speed_mode          = "binary",  # "lognormal" | "fixed" | "binary" | "grid"
  FR_speed_fixed_value   = 1.2,
  FR_speed_grid_values.  = seq(0.4, 3.0, by = 0.1),
  FR_speed_binary_values = c(0.8, 2.0),
  FR_speed_labels        = c("weak", "strong"),

  # -------------------------
  # Simulation dynamics
  # -------------------------
  max_interaction_dist = 5,      # max distance for any influence (in meters)
  detection_delay      = 5,      # sensory-motor delay (in frames)
  frame_rate           = 240,    # frames per second
  speed_decay          = 0, #0.0551, # fitted value: 0.0551 is fitted to median decay (with baseline) - was 0.0333 before, without baseline 
  t_max                = 400,    # max time to simulate (in frames)
  stall_limit          = 40      # stop the simulation if nothing happens for these many frames (based on After_FR for SRNL)
)

# Parallel settings (Mac/Linux for mclapply)
mc_cores <- max(1L, parallel::detectCores() - 2L)

# Paths (relative to project root via here package)
library(here)
data_path         <- here::here("data")
modelfits_path    <- here::here("outputs", "model_fitting")
dataanalysis_path <- here::here("outputs", "data_analysis")

# All fish observations — FR speeds are extracted from here (Response_Category == "FR")
FR_information     <- file.path(data_path, "all_observations.csv")
# Fitted model parameters (output of scripts/04_model_fitting.R)
fits_csv           <- file.path(modelfits_path, "model_comparisons.csv")

# GMM parameters by group size (output of notebooks/03_data_analysis.ipynb)
gmm_csv            <- file.path(dataanalysis_path, "GMM_values_by_size.csv")

# 3D positions from video processing (output of notebooks/01_video_processing.ipynb)
real_positions_csv <- file.path(data_path, "derived", "all_positions.csv")

stopifnot(file.exists(FR_information), file.exists(fits_csv))
if (EXP$network_source == "gmm") {
  stopifnot(file.exists(gmm_csv))
}
if (EXP$network_source == "empirical") {
  stopifnot(file.exists(real_positions_csv))
}

# Output base (separate folders so you don't mix results)
out_base <- if (EXP$network_source == "gmm") {
  here::here("outputs", "simulations", "GMM_3D")
} else {
  here::here("outputs", "simulations", "Real_3D")
}
if (!dir.exists(out_base)) dir.create(out_base, recursive = TRUE)
