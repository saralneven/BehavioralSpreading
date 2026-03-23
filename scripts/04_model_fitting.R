# 04_model_fitting.R
# Fits all model variants (combinations of time decay, neighbour proximity,
# speed, and pooling) and saves a comparison table (AIC, R2, parameters).
#
# Input:  data/derived/model_input.csv       (from 02_prepare_model_input.ipynb)
# Output: outputs/model_fitting/model_comparisons.csv

library(here)

suppressPackageStartupMessages({
  library(dplyr)
})

# ================================
# PATHS
# ================================
IN_PATH  <- here::here("data", "derived")
OUT_PATH <- here::here("outputs", "model_fitting")

# ================================
# HELPERS
# ================================
clean_numeric_column <- function(column) as.numeric(gsub(",", "", column))

# ================================
# LOAD + CLEAN
# ================================
df <- read.csv(file.path(IN_PATH, "model_input.csv"),
               sep = ",", dec = ".", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::filter(Focal_Response_Category %in% c("SRNL", "NRNL_r")) %>%
  dplyr::mutate(
    Input_Responseframe_LF = clean_numeric_column(Input_Responseframe_cam3),
    focal_Response_frame   = clean_numeric_column(Focal_Responseframe_cam3),
    Input_Speed            = clean_numeric_column(Input_Speed),
    Input_Distance         = clean_numeric_column(Input_Distance),
    Focal_Speed            = clean_numeric_column(Focal_Speed),
    Focal_Distance_Coral   = clean_numeric_column(Focal_Distance_Coral),
    Focal_Neighbour_Proximity    = clean_numeric_column(Focal_Neighbour_Proximity)
  )

# ================================
# FRAME WINDOW PER VIDEO
# ================================
s_times <- df %>%
  dplyr::group_by(Video_ID) %>%
  dplyr::summarize(
    min = min(Input_Responseframe_LF, na.rm = TRUE),
    max = max(Input_Responseframe_LF, na.rm = TRUE) + 30,
    .groups = "drop"
  )

# ================================
# PER-FRAME BUILDER
#  Distance is ALWAYS included as 1/d^2 (p_dist = 2).
#  Poolings here: sum, mean, max, nearest (NO 'rescale' here).
# ================================
calculate_decision_and_response <- function(df, s_times, tau = 5, b0, b1, c,
                                            detection_delay = 5,
                                            include_speed = 1,
                                            include_density = 1,
                                            pooling = "sum") {
  unique_fish <- unique(df[, c('Video_ID', 'Focal_Fish_ID')])
  D_i_ts <- responses <- local_densities <- focal_speeds <- focal_coral_ms <- t_rel_all <- list()
  segment_ids <- list()
  seg_counter <- 0

  for (i in seq_len(nrow(unique_fish))) {
    video_id      <- unique_fish$Video_ID[i]
    focal_fish_id <- unique_fish$Focal_Fish_ID[i]

    focal_rows <- df$Video_ID == video_id & df$Focal_Fish_ID == focal_fish_id
    focal_cat  <- df$Focal_Response_Category[focal_rows][1]

    start_frame <- s_times$min[s_times$Video_ID == video_id]
    end_frame   <- if (focal_cat %in% c("SRNL", "SRL")) {
      df$focal_Response_frame[focal_rows][1]
    } else {
      s_times$max[s_times$Video_ID == video_id]
    }

    if (!is.finite(start_frame) || !is.finite(end_frame) || start_frame >= end_frame) next

    input_data <- df[focal_rows, ]
    input_data <- input_data[input_data$Input_Response_Category %in% c("SRNL", "SRL", "FR"), ]

    speeds          <- input_data$Input_Speed
    distances       <- input_data$Input_Distance
    response_frames <- input_data$Input_Responseframe_LF + detection_delay

    local_density   <- input_data$Focal_Neighbour_Proximity[1]
    focal_speed     <- input_data$Focal_Speed[1]
    focal_coral_m   <- input_data$Focal_Distance_Coral[1] / 1000  # meters

    frame_list    <- start_frame:end_frame
    nT            <- length(frame_list)
    D_i_t_list    <- numeric(nT)
    response_list <- rep(0, nT)
    t_rel_list    <- seq_len(nT) - 1

    if (focal_cat %in% c("SRNL", "SRL")) {
      response_list[nT] <- 1
    }

    for (t in seq_len(nT)) {
      current_frame <- frame_list[t]
      valid <- which(current_frame >= response_frames)

      if (length(valid) > 0) {
        times <- (current_frame - response_frames[valid])
        d_m   <- pmax(distances[valid] / 1000, 1e-9)

        # distance always: 1/d^2
        dist_fac <- 1 / (d_m ^ 2)

        # optional speed
        if (include_speed == 1) {
          s_mps   <- pmax(speeds[valid] / 1000, 0)
          spd_fac <- (s_mps ^ c)
        } else {
          spd_fac <- 1
        }

        w <- exp(-times / tau) * (spd_fac * dist_fac)

        pooled <- switch(pooling,
                         "sum"     = sum(w),
                         "max"     = max(w),
                         "nearest" = { j <- which.min(d_m); w[j] },
                         "mean"    = mean(w),
                         stop("Unknown pooling: ", pooling))

        D_i_t_list[t] <- ifelse(is.finite(pooled), pooled, 0)
      } else {
        D_i_t_list[t] <- 0
      }
    }

    seg_counter <- seg_counter + 1
    seg_ids_vec <- rep(seg_counter, nT)

    D_i_ts[[i]]          <- D_i_t_list
    responses[[i]]       <- response_list
    t_rel_all[[i]]       <- t_rel_list
    local_densities[[i]] <- rep(local_density, nT)
    focal_speeds[[i]]    <- rep(focal_speed, nT)
    focal_coral_ms[[i]]  <- rep(focal_coral_m, nT)
    segment_ids[[i]]     <- seg_ids_vec
  }

  list(
    D_i_t          = unlist(D_i_ts),
    response       = unlist(responses),
    t_rel          = unlist(t_rel_all),
    local_density  = unlist(local_densities),
    focal_speed    = unlist(focal_speeds),
    focal_coral_m  = unlist(focal_coral_ms),
    segment_id     = unlist(segment_ids)
  )
}

# ================================
# OBJECTIVE (HAZARD ONLY; no rescale here)
#  b0 is ONLY used/estimated when include_density == 1.
# ================================
objective_function <- function(params, df, s_times, settings) {
  idx <- 1
  if (settings$include_time == 1) { tau <- max(1e-6, params[idx]); idx <- idx + 1 } else { tau <- 1e6 }

  # b0 present ONLY for density models (DiTDe / DiTDeS)
  if (settings$include_density == 1) { b0 <- params[idx]; idx <- idx + 1 } else { b0 <- 0 }

  # Gain on social input
  b1 <- params[idx]; idx <- idx + 1

  # Speed exponent
  if (settings$include_speed == 1) { c <- params[idx]; idx <- idx + 1 } else { c <- 1 }

  res <- calculate_decision_and_response(
    df, s_times, tau = tau, b0 = b0, b1 = b1, c = c,
    detection_delay   = settings$detection_delay,
    include_speed     = settings$include_speed,
    include_density   = settings$include_density,
    pooling           = settings$pooling
  )

  y  <- res$response
  ld <- res$local_density

  # ---- HAZARD LINK ----
  D_in <- pmax(res$D_i_t, 0)
  lambda <- 1e-5 + b1 * D_in
  if (settings$include_density == 1) {
    lambda <- lambda * exp(b0 * ld)
  }
  lambda <- pmax(lambda, 0)

  p <- 1 - exp(-lambda)

  eps <- 1e-12
  p <- pmin(pmax(p, eps), 1 - eps)
  -2 * sum(y * log(p) + (1 - y) * log(1 - p))
}

# ================================
# MODEL GRID (NO 'rescale')
# ================================
pooling_types <- c("sum","max","nearest","mean")

core_structs <- list(
  Di      = list(include_time=0, include_density=0, include_speed=0),
  DiT     = list(include_time=1, include_density=0, include_speed=0),
  DiTDe   = list(include_time=1, include_density=1, include_speed=0),
  DiTDeS  = list(include_time=1, include_density=1, include_speed=1),
  DiS     = list(include_time=0, include_density=0, include_speed=1),
  DiDeS   = list(include_time=0, include_density=1, include_speed=1),
  DiDe    = list(include_time=0, include_density=1, include_speed=0),
  DiTS    = list(include_time=1, include_density=0, include_speed=1)
)

grid <- do.call(
  rbind,
  lapply(names(core_structs), function(core) {
    base <- core_structs[[core]]
    expand.grid(
      model_core        = core,
      include_time      = base$include_time,
      include_density   = base$include_density,
      include_speed     = base$include_speed,
      pooling           = pooling_types,
      KEEP.OUT.ATTRS    = FALSE,
      stringsAsFactors  = FALSE
    )
  })
)

name_model <- function(row){
  suffix <- switch(row$pooling,
                   sum="p1", max="sel-max", nearest="sel-near", mean="avg")
  paste0(row$model_core, "_", suffix)
}

model_variants <- setNames(
  lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ])),
  sapply(split(grid, seq_len(nrow(grid))), name_model)
)

# ================================
# NULL DEVIANCE (for R2)
# ================================
tmp <- calculate_decision_and_response(
  df, s_times,
  tau = 10, b0 = 0, b1 = 0, c = 0,
  detection_delay = 5,
  include_speed   = 1,
  include_density = 1,
  pooling         = "sum"
)
resp_series   <- tmp$response
p0            <- min(max(mean(resp_series), 1e-12), 1 - 1e-12)
null_deviance <- -2 * sum(resp_series * log(p0) + (1 - resp_series) * log(1 - p0))

# ================================
# FIT LOOP
# ================================
results <- data.frame()

for (model_name in names(model_variants)) {
  cat("\n--- Running model:", model_name, "---\n")
  settings <- model_variants[[model_name]]
  settings$detection_delay <- 5

  # Build init ONLY with params that are actually used
  init <- numeric(0)
  if (settings$include_time == 1)  init <- c(init, tau = 12.0)
  if (settings$include_density == 1) {
    init <- c(init, b0 = -5.0e-3)            # density coefficient (hazard)
  }
  init <- c(init, b1 = 3.0e-3)
  if (settings$include_speed == 1) init <- c(init, c = 1.8)

  names(init) <- c(
    if (settings$include_time == 1)     "tau",
    if (settings$include_density == 1)  "b0",
    "b1",
    if (settings$include_speed == 1)    "c"
  )

  fit <- tryCatch({
    # Build bounds matching init order
    lower <- numeric(0)
    upper <- numeric(0)

    if (settings$include_time == 1) {
      lower <- c(lower, tau = 1e-6)
      upper <- c(upper, tau = 1e6)
    }
    if (settings$include_density == 1) {
      # wide bounds; adjust if you have a better biological range
      lower <- c(lower, b0 = -1)
      upper <- c(upper, b0 =  1)
    }

    # b1 must be >= 0
    lower <- c(lower, b1 = 0)
    upper <- c(upper, b1 = 1e6)

    if (settings$include_speed == 1) {
      # c >= 0; you can cap upper if you want
      lower <- c(lower, c = 0)
      upper <- c(upper, c = 10)
    }

    # Make sure names align with init
    stopifnot(identical(names(init), names(lower)))
    stopifnot(identical(names(init), names(upper)))

    optim_result <- optim(
      par = init,
      fn = objective_function,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      df = df,
      s_times = s_times,
      settings = settings,
      control = list(maxit = 2000)
    )


    # parameter count k (for AIC) — count only included params
    k <- (settings$include_time == 1) +          # tau
         (settings$include_density == 1) +       # b0 (only if density)
         1 +                                     # b1
         (settings$include_speed == 1)           # c

    dev <- optim_result$value
    aic <- dev + 2 * k
    r2  <- 1 - (dev / null_deviance)

    pars <- optim_result$par
    
    if (any(abs(pars - lower) < 1e-10) || any(abs(pars - upper) < 1e-10)) {
      message("  Note: parameters at bounds in ", model_name)
    }

    getp <- function(nm) if (nm %in% names(pars)) as.numeric(pars[[nm]]) else NA_real_

    data.frame(
      model   = model_name,
      core    = settings$model_core,
      pooling = settings$pooling,
      tau     = getp("tau"),
      b0      = getp("b0"),      # NA for DiT
      b1      = getp("b1"),
      c       = getp("c"),
      deviance = dev,
      AIC      = aic,
      R2       = r2,
      Include_Time    = settings$include_time,
      Include_Density = settings$include_density,
      Include_Speed   = settings$include_speed,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning(paste("Model", model_name, "failed:", e$message))
    data.frame(
      model = model_name, core = settings$model_core, pooling = settings$pooling,
      tau = NA, b0 = NA, b1 = NA, c = NA,
      deviance = NA, AIC = NA, R2 = NA,
      Include_Time      = settings$include_time,
      Include_Density   = settings$include_density,
      Include_Speed     = settings$include_speed,
      stringsAsFactors = FALSE
    )
  })

  results <- dplyr::bind_rows(results, fit)
}

# ================================
# SAVE
# ================================
outfile <- file.path(OUT_PATH, "model_comparisons.csv")
write.csv(results, outfile, row.names = FALSE)
message("Wrote ModelComparisons: ", outfile)

print(results)
