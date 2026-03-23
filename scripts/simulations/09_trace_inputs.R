library(here)

# 09_trace_inputs_DiTDeS.R
# Trace social input over time to the k-nearest neighbours of the primary FR
# under the DiTDeS model.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(parallel)
  library(Matrix)
})

# -------------------------
# 0) USER SETTINGS
# -------------------------
# IMPORTANT: set this to your simulations code folder
src_dir <- here::here("scripts", "simulations")

# Where to save output traces (a new folder will be created)
out_base_trace <- here::here("outputs", "simulations", "GMM_3D", "diagnostics")

# ---- extra mechanistic output (Top-K contributors at response time) ----
STORE_CONTRIB_TOPK <- TRUE
TOPK_CONTRIB <- 40   # store top 40 contributing responders per responding target fish (at its response time)

# Trace settings (keep small-ish: tracing is heavy!)
TRACE <- list(
  model_name   = "DiTDeS",
  areas        = c(1),
  group_sizes  = c(8, 30),      # c(5, 8, 10, 12, 15, 20, 25, 30, 40),
  N_FRs        = 1L,
  num_reps     = 500,          # reps per condition
  k_neighbours = NA_integer_,  # <-- set below (cannot refer to TRACE inside its own list)
  record_every = 1L,           # record every X frames (1 = every frame)
  t_max        = 400L,         # override EXP$t_max if you want
  stall_limit  = 40L,          # override EXP$stall_limit if you want
  mc_cores     = max(1L, detectCores() - 2L)
)

# Set k-neighbours AFTER TRACE exists
TRACE$k_neighbours <- max(TRACE$group_sizes)

# Choose network source: "gmm" or "empirical"
FORCE_NETWORK_SOURCE <- NA_character_  # e.g. "gmm" or "empirical" or NA to keep EXP$network_source

# Use binary weak/strong initiator speeds 
FORCE_FR_SPEED_MODE <- "binary"  # "lognormal" | "fixed" | "binary" | "grid"

# -------------------------
# 1) SOURCE MODULES
# -------------------------
source(file.path(src_dir, "00_config.R"))
source(file.path(src_dir, "01_io_and_helpers.R"))
source(file.path(src_dir, "02_load_model_fits.R"))
source(file.path(src_dir, "03_fr_speed.R"))
source(file.path(src_dir, "04_gmm_loader.R"))
source(file.path(src_dir, "04_empirical_loader.R"))
source(file.path(src_dir, "05_geometry.R"))
source(file.path(src_dir, "06_simulation_core.R"))

# -------------------------
# 2) APPLY TRACE OVERRIDES
# -------------------------
EXP$t_max       <- as.integer(TRACE$t_max)
EXP$stall_limit <- as.integer(TRACE$stall_limit)

if (!is.na(FORCE_NETWORK_SOURCE)) {
  EXP$network_source <- FORCE_NETWORK_SOURCE
}

# Apply requested mode (binary)
EXP$FR_speed_mode <- FORCE_FR_SPEED_MODE

# If binary FR-speed: derive weak/strong from DATA (same as 07_run_experiment)
if (EXP$FR_speed_mode == "binary") {
  ws <- get_weak_strong_FR_speeds(
    path  = FR_information,
    probs = c(0.25, 0.75)
  )
  print(ws)

  EXP$FR_speed_binary_values <- ws$FR_speed_m_per_s
  EXP$FR_speed_labels        <- ws$label
}

cat("Binary FR speeds used (m/s):\n")
print(data.frame(label = EXP$FR_speed_labels,
                 speed = EXP$FR_speed_binary_values))

# Output folder
run_tag <- paste0(
  "TRACEINPUTS_",
  format(Sys.time(), "%Y%m%d-%H%M%S"),
  "_", TRACE$model_name,
  "_NFR", TRACE$N_FRs,
  "_", EXP$network_source,
  "_dec", EXP$speed_decay
)
safe_path <- file.path(out_base_trace, run_tag)
if (!dir.exists(safe_path)) dir.create(safe_path, recursive = TRUE)
log_msg("Trace output folder: %s", safe_path)

# -------------------------
# 3) LOAD FITTED MODEL PARAMS (DiTDeS)
# -------------------------
mf <- load_model_fits(fits_csv, EXP)
get_model_spec <- mf$get_model_spec
spec <- get_model_spec(TRACE$model_name)
p <- spec$params

# speed dist (only used for non-FR fish; FR is forced by label/value)
speed_dist <- fit_fr_speed_lognormal(FR_information)

# -------------------------
# 4) PREP NETWORK SAMPLERS
# -------------------------
get_gmm_for_size <- NULL
sample_from_gmm_3d <- NULL
sample_real_positions_3d <- NULL

if (EXP$network_source == "gmm") {
  gmm_obj <- load_gmms_by_size(gmm_csv)
  get_gmm_for_size   <- gmm_obj$get_gmm_for_size
  sample_from_gmm_3d <- gmm_obj$sample_from_gmm_3d
} else if (EXP$network_source == "empirical") {
  emp <- load_empirical_positions_3d(real_positions_csv, EXP)
  sample_real_positions_3d <- emp$sample_real_positions_3d
} else {
  stop("Unknown EXP$network_source: ", EXP$network_source)
}

# -------------------------
# 5) CORE: simulate + trace inputs for k-nearest neighbours
# -------------------------
simulate_trace_on_positions <- function(
  positions, Area, group_size_requested,
  N_FRs, rep_id,
  fr_speed_label = NA_character_,
  fr_speed_value = NA_real_
) {
  N <- nrow(positions)
  G <- make_distance_weights(positions, EXP$max_interaction_dist)
  W <- G$W; LD <- G$LD; D <- G$D

  # --- seed FR(s): match 06_simulation_core exactly ---
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

  # --- speeds: use your existing helper, with forced FR speed ---
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

  # --- choose k-nearest neighbours of the primary FR at t=0 ---
  cand <- setdiff(seq_len(N), FR_primary)
  cand <- cand[is.finite(D[cand, FR_primary])]
  ord  <- cand[order(D[cand, FR_primary])]
  k    <- min(as.integer(TRACE$k_neighbours), length(ord))
  target_ids <- ord[seq_len(k)]
  target_k   <- seq_len(k)
  target_d   <- D[target_ids, FR_primary]

  # storage: time traces
  trace_rows <- list()
  ridx <- 1L

  # storage: top-k contributors at response time
  contrib_rows <- list()
  cidx <- 1L

  # Helper: compute top-K contributor decomposition at a focal response time
  # Uses:
  #   - responders with dt >= detection_delay
  #   - time decay exp(-dt_eff / tau)
  #   - speed decay exp(-speed_decay * s_times)^c
  #   - distance weight from W[focal, j]
  compute_topk_contrib <- function(focal_i, t_now, startle_times, W, speed0, FR_IDs, p, include_speed, TOPK=10L) {

    # responders that are detectable at t_now
    active <- which(!is.na(startle_times) & startle_times < (t_now - p$detection_delay))
    if (!length(active)) return(NULL)

    s_times <- startle_times[active]
    dt_eff  <- t_now - (s_times + p$detection_delay)

    time_w <- exp(-dt_eff / p$tau)

    speed_fac <- if (include_speed == 1) {
      (speed0[active] * exp(-p$speed_decay * s_times))^p$c
    } else {
      rep(1.0, length(active))
    }

    w_dist <- as.numeric(W[focal_i, active])

    contrib <- w_dist * time_w * speed_fac
    contrib[!is.finite(contrib)] <- NA_real_

    ord <- order(contrib, decreasing = TRUE, na.last = NA)
    if (!length(ord)) return(NULL)

    keep <- ord[seq_len(min(TOPK, length(ord)))]

    data.frame(
      neigh_id     = as.integer(active[keep]),
      neigh_is_FR  = as.integer(active[keep] %in% FR_IDs),
      s_time       = as.numeric(s_times[keep]),
      dt_eff       = as.numeric(dt_eff[keep]),
      w_dist       = as.numeric(w_dist[keep]),
      time_w       = as.numeric(time_w[keep]),
      speed_fac    = as.numeric(speed_fac[keep]),
      contrib      = as.numeric(contrib[keep]),
      rank_contrib = as.integer(seq_along(keep)),
      stringsAsFactors = FALSE
    )
  }

  stalled <- 0L
  prev_n  <- sum(!is.na(startle_times))

  for (t in seq_len(EXP$t_max)) {

    # ----------------------------
    # Inputs: (1) all responders, (2) FR-only baseline, (3) marginal extra
    # ----------------------------
    D_all <- get_inputs_fast(
      t = t,
      startle_times = startle_times,
      W = W,
      include_speed = spec$include_speed,
      tau = p$tau,
      c = p$c,
      speed0 = speed0,
      detection_delay = p$detection_delay,
      speed_decay = p$speed_decay
    )

    # FR-only baseline: keep only initiator(s) in startle_times
    startle_FRonly <- rep(NA_real_, N)
    startle_FRonly[FR_IDs] <- startle_times[FR_IDs]  # typically 0

    D_FRonly <- get_inputs_fast(
      t = t,
      startle_times = startle_FRonly,
      W = W,
      include_speed = spec$include_speed,
      tau = p$tau,
      c = p$c,
      speed0 = speed0,
      detection_delay = p$detection_delay,
      speed_decay = p$speed_decay
    )

    D_delta <- D_all - D_FRonly

    # (keep name compatible with rest of code)
    D_i_t <- D_all

    # hazard (same as core)
    lambda <- 1e-5 + p$b1 * pmax(D_i_t, 0)
    if (spec$include_density == 1) lambda <- lambda * exp(p$b0 * LD)

    p_vec <- 1 - exp(-lambda)
    p_vec <- pmin(pmax(p_vec, 1e-12), 1 - 1e-12)

    # --- record traces every TRACE$record_every frames ---
    if ((t %% as.integer(TRACE$record_every)) == 0L) {
      st <- startle_times[target_ids]

      n_resp_total <- sum(!is.na(startle_times) & startle_times <= t)
      n_resp_nonFR <- max(0L, n_resp_total - length(FR_IDs))

      trace_rows[[ridx]] <- data.frame(
        Rep_ID    = as.integer(rep_id),
        Area      = as.numeric(Area),
        N         = as.integer(N),
        N_req     = as.integer(group_size_requested),

        Model     = TRACE$model_name,
        Network   = as.character(EXP$network_source),
        N_FRs     = as.integer(N_FRs),

        FR_ID     = as.integer(FR_primary),
        FR_speed  = as.numeric(fr_speed0),
        FR_label  = as.character(fr_speed_label),

        t         = as.integer(t),

        k         = as.integer(target_k),
        Fish_ID   = as.integer(target_ids),
        d_to_FR   = as.numeric(target_d),

        Input_D_all    = as.numeric(D_all[target_ids]),
        Input_D_FRonly = as.numeric(D_FRonly[target_ids]),
        Input_D_delta  = as.numeric(D_delta[target_ids]),

        Lambda    = as.numeric(lambda[target_ids]),
        p         = as.numeric(p_vec[target_ids]),

        Responded = as.integer(!is.na(st) & st <= t),
        Startle_t = as.numeric(st),

        N_resp_total = as.integer(n_resp_total),
        N_resp_nonFR = as.integer(n_resp_nonFR),

        stringsAsFactors = FALSE
      )
      ridx <- ridx + 1L
    }

    # update responses
    new_ids <- which(is.na(startle_times) & (runif(N) < p_vec))

    if (length(new_ids)) {
      startle_times[new_ids] <- t

      # if any traced target fish responds NOW, store top-K contributors at this moment
      if (STORE_CONTRIB_TOPK) {
        new_targets <- intersect(new_ids, target_ids)

        if (length(new_targets)) {
          for (focal_i in new_targets) {

            topdf <- compute_topk_contrib(
              focal_i = focal_i,
              t_now = t,
              startle_times = startle_times,
              W = W,
              speed0 = speed0,
              FR_IDs = FR_IDs,
              p = p,
              include_speed = spec$include_speed,
              TOPK = TOPK_CONTRIB
            )

            if (!is.null(topdf) && nrow(topdf)) {
              topdf$Rep_ID    <- as.integer(rep_id)
              topdf$Area      <- as.numeric(Area)
              topdf$N         <- as.integer(N)
              topdf$N_req     <- as.integer(group_size_requested)
              topdf$Network   <- as.character(EXP$network_source)
              topdf$Model     <- TRACE$model_name
              topdf$N_FRs     <- as.integer(N_FRs)

              topdf$FR_ID     <- as.integer(FR_primary)
              topdf$FR_label  <- as.character(fr_speed_label)
              topdf$FR_speed  <- as.numeric(fr_speed0)

              topdf$Fish_ID   <- as.integer(focal_i)
              topdf$startle_t <- as.numeric(t)

              topdf$d_to_FR   <- as.numeric(D[focal_i, FR_primary])

              contrib_rows[[cidx]] <- topdf
              cidx <- cidx + 1L
            }
          }
        }
      }
    }

    # stall logic
    now_n <- sum(!is.na(startle_times))
    if (now_n == prev_n) {
      stalled <- stalled + 1L
      if (stalled >= EXP$stall_limit) break
    } else {
      stalled <- 0L
      prev_n  <- now_n
    }
  }

  # summary per target fish 
  trace_df <- bind_rows(trace_rows)

  summ_df <- trace_df %>%
    group_by(Rep_ID, Area, N, N_req, Model, Network, N_FRs, FR_ID, FR_speed, FR_label, k, Fish_ID, d_to_FR) %>%
    summarise(
      responded   = any(Responded == 1),
      startle_t   = suppressWarnings(min(Startle_t, na.rm = TRUE)),

      max_input_all    = max(Input_D_all, na.rm = TRUE),
      max_input_FRonly = max(Input_D_FRonly, na.rm = TRUE),
      max_input_delta  = max(Input_D_delta, na.rm = TRUE),

      max_p      = max(p, na.rm = TRUE),
      max_lambda = max(Lambda, na.rm = TRUE),

      # input components at the moment this fish responds (if it responds)
      input_all_at_startle = {
        tt <- suppressWarnings(min(Startle_t, na.rm = TRUE))
        if (is.finite(tt)) {
          x <- Input_D_all[which(t == tt)]
          if (length(x)) x[1] else NA_real_
        } else NA_real_
      },
      input_FRonly_at_startle = {
        tt <- suppressWarnings(min(Startle_t, na.rm = TRUE))
        if (is.finite(tt)) {
          x <- Input_D_FRonly[which(t == tt)]
          if (length(x)) x[1] else NA_real_
        } else NA_real_
      },
      input_delta_at_startle = {
        tt <- suppressWarnings(min(Startle_t, na.rm = TRUE))
        if (is.finite(tt)) {
          x <- Input_D_delta[which(t == tt)]
          if (length(x)) x[1] else NA_real_
        } else NA_real_
      },

      .groups = "drop"
    ) %>%
    mutate(startle_t = ifelse(is.infinite(startle_t), NA_real_, startle_t)) %>%
    group_by(Rep_ID, Area, N, N_req, Model, Network, N_FRs, FR_ID, FR_speed, FR_label) %>%
    mutate(
      order_among_targets = ifelse(is.na(startle_t), NA_integer_,
                                   rank(startle_t, ties.method = "first")),
      is_first_secondary_target = !is.na(startle_t) & order_among_targets == 1L
    ) %>%
    ungroup()

  contrib_df <- if (STORE_CONTRIB_TOPK && length(contrib_rows)) bind_rows(contrib_rows) else data.frame()

  list(trace = trace_df, summary = summ_df, contrib = contrib_df)
}

# -------------------------
# 6) RUN GRID (weak/strong) × group sizes × reps
# -------------------------
fr_cond_labels <- if (EXP$FR_speed_mode == "binary") EXP$FR_speed_labels else NA_character_
stopifnot(EXP$FR_speed_mode == "binary")

worker_one <- function(A, Nreq, rep_id, fr_lab) {
  tryCatch({
    set.seed(seed_for(A, Nreq, rep_id))

    if (EXP$network_source == "gmm") {
      gmm <- get_gmm_for_size(as.integer(Nreq))
      positions <- sample_from_gmm_3d(as.integer(Nreq), gmm)
      metaN <- as.integer(Nreq)
    } else {
      geom <- sample_real_positions_3d(group_size = as.integer(Nreq), Area = A, rep_id = rep_id)
      positions <- geom$positions
      metaN <- nrow(positions)
    }

    out <- simulate_trace_on_positions(
      positions = positions,
      Area = A,
      group_size_requested = as.integer(Nreq),
      N_FRs = as.integer(TRACE$N_FRs),
      rep_id = as.integer(rep_id),
      fr_speed_label = fr_lab,
      fr_speed_value = NA_real_
    )

    out$trace$N_used <- as.integer(metaN)
    out$summary$N_used <- as.integer(metaN)
    if (!is.null(out$contrib) && nrow(out$contrib)) {
      out$contrib$N_used <- as.integer(metaN)
    }

    out
  }, error = function(e) {
    list(
      trace = data.frame(),
      summary = data.frame(
        Rep_ID = as.integer(rep_id),
        Area = as.numeric(A),
        N_req = as.integer(Nreq),
        FR_label = as.character(fr_lab),
        Error = conditionMessage(e),
        stringsAsFactors = FALSE
      ),
      contrib = data.frame()
    )
  })
}

all_trace <- list()
all_summ  <- list()
all_contrib <- list()
idx <- 1L

for (A in TRACE$areas) {
  for (Nreq in TRACE$group_sizes) {
    log_msg("Tracing: A=%.1f N=%d", A, Nreq)

    jobs <- expand.grid(
      rep_id = seq_len(TRACE$num_reps),
      fr_lab = fr_cond_labels,
      stringsAsFactors = FALSE
    )

    res_list <- mclapply(
      seq_len(nrow(jobs)),
      function(i) worker_one(A, Nreq, jobs$rep_id[i], jobs$fr_lab[i]),
      mc.cores = TRACE$mc_cores
    )

    tr <- bind_rows(lapply(res_list, `[[`, "trace"))
    su <- bind_rows(lapply(res_list, `[[`, "summary"))
    co <- bind_rows(lapply(res_list, `[[`, "contrib"))

    # Save per-condition chunks
    f_tr <- file.path(safe_path, sprintf("trace_A%.1f_N%d.csv", A, Nreq))
    f_su <- file.path(safe_path, sprintf("summary_A%.1f_N%d.csv", A, Nreq))
    data.table::fwrite(tr, f_tr)
    data.table::fwrite(su, f_su)
    log_msg("  saved: %s (%d rows)", basename(f_tr), nrow(tr))
    log_msg("  saved: %s (%d rows)", basename(f_su), nrow(su))

    if (STORE_CONTRIB_TOPK) {
      f_co <- file.path(safe_path, sprintf("contrib_A%.1f_N%d.csv", A, Nreq))
      data.table::fwrite(co, f_co)
      log_msg("  saved: %s (%d rows)", basename(f_co), nrow(co))
    }

    all_trace[[idx]] <- tr
    all_summ[[idx]]  <- su
    all_contrib[[idx]] <- co
    idx <- idx + 1L

    gc(FALSE)
  }
}

trace_df <- bind_rows(all_trace)
summ_df  <- bind_rows(all_summ)
contrib_df <- bind_rows(all_contrib)

data.table::fwrite(trace_df, file.path(safe_path, "TRACE_ALL.csv"))
data.table::fwrite(summ_df,  file.path(safe_path, "SUMMARY_ALL.csv"))
if (STORE_CONTRIB_TOPK) {
  data.table::fwrite(contrib_df, file.path(safe_path, "CONTRIB_TOPK_ALL.csv"))
}

log_msg("DONE. Combined files:")
log_msg("  %s", file.path(safe_path, "TRACE_ALL.csv"))
log_msg("  %s", file.path(safe_path, "SUMMARY_ALL.csv"))
if (STORE_CONTRIB_TOPK) {
  log_msg("  %s", file.path(safe_path, "CONTRIB_TOPK_ALL.csv"))
}
