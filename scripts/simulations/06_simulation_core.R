# 06_simulation_core.R

get_inputs_fast <- function(t, startle_times, W, include_speed, tau, c,
                            speed0, detection_delay, speed_decay) {
  n     <- length(startle_times)
  valid <- which(!is.na(startle_times) & startle_times < (t - detection_delay))
  if (length(valid) == 0) return(numeric(n))

  s_times <- startle_times[valid]

  speed_fac <- if (include_speed == 1) {
    (speed0[valid] * exp(-speed_decay * s_times))^c
  } else {
    1.0
  }

  time_kernel <- exp(-((t - (s_times + detection_delay)) / tau)) * speed_fac
  as.numeric(W[, valid, drop = FALSE] %*% time_kernel)
}

# ------------------------------------------------------------
# Shared core: simulate on GIVEN positions
# ------------------------------------------------------------
simulate_condition_on_positions <- function(
  model_name, Area, group_size_requested, positions, meta,
  N_FRs, speed_dist, EXP, get_model_spec,
  fr_speed_label = NA_character_,
  fr_speed_value = NA_real_
) {
  spec <- get_model_spec(model_name)
  p    <- spec$params

  N <- nrow(positions)
  if (N < 2) stop("Too few fish in positions (N < 2).")

  G <- make_distance_weights(positions, EXP$max_interaction_dist)
  W <- G$W; LD <- G$LD; D <- G$D

  # --- seed FR(s) ---
  startle_times <- rep(NA_real_, N)
  FR_primary <- sample.int(N, 1)
  startle_times[FR_primary] <- 0
  FR_IDs <- FR_primary

  if (N_FRs > 1) {
    cand <- setdiff(seq_len(N), FR_primary)
    nn <- cand[order(D[cand, FR_primary], na.last = TRUE)][1:(N_FRs - 1)]
    nn <- nn[is.finite(nn)]
    if (length(nn)) {
      startle_times[nn] <- 0
      FR_IDs <- c(FR_primary, nn)
    }
  }

  # --- speeds ---
  spd <- make_speed0(
    N = N,
    fr_index = FR_primary,
    mu = speed_dist$mu,
    sigma = speed_dist$sigma,
    EXP = EXP,
    fr_speed_label = fr_speed_label,
    fr_speed_value = fr_speed_value
  )
  speed0    <- spd$speed0
  fr_speed0 <- spd$fr_speed0

  # --- FR LD stats ---
  FR_LD_primary    <- LD[FR_primary]
  FR_LD_percentile <- mean(LD <= FR_LD_primary)
  FR_LD_z          <- (FR_LD_primary - mean(LD)) / (sd(LD) + 1e-12)
  FR_LD_meanFR     <- mean(LD[FR_IDs])

  stalled <- 0L
  prev_n  <- sum(!is.na(startle_times))

  for (t in seq_len(EXP$t_max)) {
    D_i_t <- get_inputs_fast(
      t, startle_times, W,
      spec$include_speed, p$tau, p$c,
      speed0, p$detection_delay, p$speed_decay
    )

    lambda <- 1e-5 + p$b1 * pmax(D_i_t, 0)
    if (spec$include_density == 1) lambda <- lambda * exp(p$b0 * LD)

    p_vec <- 1 - exp(-lambda)
    p_vec <- pmin(pmax(p_vec, 1e-12), 1 - 1e-12)

    new_ids <- which(is.na(startle_times) & (runif(N) < p_vec))
    if (length(new_ids)) startle_times[new_ids] <- t

    now_n <- sum(!is.na(startle_times))
    if (now_n == prev_n) {
      stalled <- stalled + 1L
      if (stalled >= EXP$stall_limit) break
    } else {
      stalled <- 0L
      prev_n  <- now_n
    }
  }

  started <- which(!is.na(startle_times))

  # Build output row (keep your GMM columns + add empirical metadata columns as NA when absent)
  data.frame(
    Model            = model_name,
    Area             = Area,
    Density          = as.integer(group_size_requested),  # keep name "Density" for compatibility

    Include_Speed    = spec$include_speed,
    Include_Density  = spec$include_density,
    Include_Time     = spec$include_time,

    N_FRs            = as.integer(N_FRs),
    FR_speed         = as.numeric(fr_speed0),
    FR_speed_mode    = as.character(EXP$FR_speed_mode),
    FR_speed_label   = as.character(fr_speed_label),

    N_Fish           = as.integer(N),
    N_Startles       = as.integer(length(started)),

    FR_LD_primary    = as.numeric(FR_LD_primary),
    FR_LD_percentile = as.numeric(FR_LD_percentile),
    FR_LD_meanFR     = as.numeric(FR_LD_meanFR),
    FR_LD_z          = as.numeric(FR_LD_z),

    # empirical-only metadata (NA for GMM unless provided)
    Video_ID             = if (!is.null(meta$Video_ID)) as.integer(meta$Video_ID) else NA_integer_,
    Group_Size_Used      = if (!is.null(meta$Group_Size_Used)) as.integer(meta$Group_Size_Used) else NA_integer_,
    Group_Size_Target    = if (!is.null(meta$Group_Size_Target)) as.integer(meta$Group_Size_Target) else NA_integer_,
    Group_Size_Requested = if (!is.null(meta$Group_Size_Requested)) as.integer(meta$Group_Size_Requested) else as.integer(group_size_requested),

    Error = NA_character_,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Wrapper A: GMM networks (your existing entry point)
# ------------------------------------------------------------
simulate_condition_on_geometry <- function(
  model_name, Area, group_size, N_FRs, num_reps,
  speed_dist, EXP,
  get_model_spec,
  get_gmm_for_size,
  sample_from_gmm_3d,
  mc_cores,
  fr_speed_label = NA_character_,
  fr_speed_value = NA_real_
) {
  worker <- function(rep_id) {
    tryCatch({
      N <- as.integer(group_size)
      set.seed(seed_for(Area, N, rep_id))

      gmm <- get_gmm_for_size(N)
      positions <- sample_from_gmm_3d(N, gmm)

      simulate_condition_on_positions(
        model_name = model_name,
        Area = Area,
        group_size_requested = N,
        positions = positions,
        meta = list(),
        N_FRs = N_FRs,
        speed_dist = speed_dist,
        EXP = EXP,
        get_model_spec = get_model_spec,
        fr_speed_label = fr_speed_label,
        fr_speed_value = fr_speed_value
      )
    }, error = function(e) {
      spec <- get_model_spec(model_name)
      data.frame(
        Model = model_name, Area = Area, Density = as.integer(group_size),
        Include_Speed = spec$include_speed, Include_Density = spec$include_density, Include_Time = spec$include_time,
        N_FRs = as.integer(N_FRs), FR_speed = NA_real_, FR_speed_mode = as.character(EXP$FR_speed_mode),
        FR_speed_label = as.character(fr_speed_label),
        N_Fish = NA_integer_, N_Startles = NA_integer_,
        FR_LD_primary = NA_real_, FR_LD_percentile = NA_real_, FR_LD_meanFR = NA_real_, FR_LD_z = NA_real_,
        Video_ID = NA_integer_, Group_Size_Used = NA_integer_, Group_Size_Target = NA_integer_, Group_Size_Requested = as.integer(group_size),
        Error = paste0("rep_id=", rep_id, " : ", conditionMessage(e)),
        stringsAsFactors = FALSE
      )
    })
  }

  parallel::mclapply(seq_len(num_reps), worker, mc.cores = mc_cores) |>
    dplyr::bind_rows()
}

# ------------------------------------------------------------
# Wrapper B: Empirical networks
# ------------------------------------------------------------
simulate_condition_on_empirical_geometry <- function(
  model_name, Area, group_size, N_FRs, num_reps,
  speed_dist, EXP,
  get_model_spec,
  sample_real_positions_3d,   # from 04_empirical_loader.R
  mc_cores,
  fr_speed_label = NA_character_,
  fr_speed_value = NA_real_
) {
  worker <- function(rep_id) {
    tryCatch({
      Nreq <- as.integer(group_size)
      set.seed(seed_for(Area, Nreq, rep_id))

      geom <- sample_real_positions_3d(group_size = Nreq, Area = Area, rep_id = rep_id)
      positions <- geom$positions

      simulate_condition_on_positions(
        model_name = model_name,
        Area = Area,
        group_size_requested = Nreq,
        positions = positions,
        meta = geom,
        N_FRs = N_FRs,
        speed_dist = speed_dist,
        EXP = EXP,
        get_model_spec = get_model_spec,
        fr_speed_label = fr_speed_label,
        fr_speed_value = fr_speed_value
      )
    }, error = function(e) {
      spec <- get_model_spec(model_name)
      data.frame(
        Model = model_name, Area = Area, Density = as.integer(group_size),
        Include_Speed = spec$include_speed, Include_Density = spec$include_density, Include_Time = spec$include_time,
        N_FRs = as.integer(N_FRs), FR_speed = NA_real_, FR_speed_mode = as.character(EXP$FR_speed_mode),
        FR_speed_label = as.character(fr_speed_label),
        N_Fish = NA_integer_, N_Startles = NA_integer_,
        FR_LD_primary = NA_real_, FR_LD_percentile = NA_real_, FR_LD_meanFR = NA_real_, FR_LD_z = NA_real_,
        Video_ID = NA_integer_, Group_Size_Used = NA_integer_, Group_Size_Target = NA_integer_, Group_Size_Requested = as.integer(group_size),
        Error = paste0("rep_id=", rep_id, " : ", conditionMessage(e)),
        stringsAsFactors = FALSE
      )
    })
  }

  parallel::mclapply(seq_len(num_reps), worker, mc.cores = mc_cores) |>
    dplyr::bind_rows()
}
