# S01_model_fitting_time_dependent.R
# Supplementary: fits models using time-dependent distances (distance measured at
# each neighbour's response frame). Confirms main model results are robust.
#
# Input:  data/derived/model_input_time_dependent.csv
#         (from S01_time_dependent_distances.ipynb)
# Output: outputs/model_fitting/model_comparisons_time_dependent.csv

library(here)

library(dplyr)

# ---------- PATHS ----------
IN_PATH  <- here::here("data", "derived")
OUT_PATH <- here::here("outputs", "model_fitting")

# ---------- Helpers ----------
clean_numeric_column <- function(column) {
  column <- as.numeric(gsub(",", "", column))
  column
}

# ---------- Load + clean ----------
df <- read.csv(file.path(IN_PATH, "model_input_time_dependent.csv"),
               sep = ",", dec = ".", header = TRUE, stringsAsFactors = FALSE) %>%
  filter(Focal_Response_Category %in% c("SRNL", "NRNL_r")) %>%
  mutate(
    Input_Responseframe_LF = clean_numeric_column(Input_Responseframe_cam3),
    focal_Response_frame   = clean_numeric_column(Focal_Responseframe_cam3),
    Input_Speed            = clean_numeric_column(Input_Speed),
    Input_Distance         = clean_numeric_column(Input_Distance_at_InputResponse),
    Focal_Speed            = clean_numeric_column(Focal_Speed),
    Focal_Distance_Coral   = clean_numeric_column(Focal_Distance_Coral),
    Focal_Neighbour_Proximity    = clean_numeric_column(Focal_Neighbour_Proximity)   # precomputed density from CSV
  )

# ---------- Frame window per video ----------
s_times <- df %>%
  group_by(Video_ID) %>%
  summarize(
    min = min(Input_Responseframe_LF, na.rm = TRUE),
    max = max(Input_Responseframe_LF, na.rm = TRUE) + 30,
    .groups = "drop"
  )

# ---------- Core per-frame builder ----------
calculate_decision_and_response <- function(df, s_times, tau = 5, b0, b1, c,
                                            detection_delay = 5, include_speed = 1,
                                            include_density = 1, double_lobe = 0) {
  unique_fish <- unique(df[, c('Video_ID', 'Focal_Fish_ID')])
  D_i_ts <- responses <- local_densities <- focal_speeds <- focal_coral_ms <- list()

  for (i in seq_len(nrow(unique_fish))) {
    video_id      <- unique_fish$Video_ID[i]
    focal_fish_id <- unique_fish$Focal_Fish_ID[i]

    start_frame <- s_times$min[s_times$Video_ID == video_id]
    end_frame <- if (df$Focal_Response_Category[df$Video_ID == video_id & df$Focal_Fish_ID == focal_fish_id][1] %in% c("SRNL", "SRL")) {
      df$focal_Response_frame[df$Video_ID == video_id & df$Focal_Fish_ID == focal_fish_id][1]
    } else {
      s_times$max[s_times$Video_ID == video_id]
    }

    input_data <- df[df$Video_ID == video_id & df$Focal_Fish_ID == focal_fish_id, ]
    input_data <- input_data[input_data$Input_Response_Category %in% c("SRNL", "SRL", "FR"), ]

    speeds          <- input_data$Input_Speed
    distances       <- input_data$Input_Distance
    response_frames <- input_data$Input_Responseframe_LF + detection_delay

    neighbour_proximity   <- input_data$Focal_Neighbour_Proximity[1]                 # from CSV
    focal_speed     <- input_data$Focal_Speed[1]
    focal_coral_m   <- input_data$Focal_Distance_Coral[1] / 1000   # meters

    frame_list    <- start_frame:end_frame
    D_i_t_list    <- numeric(length(frame_list))
    response_list <- rep(0, length(frame_list))

    if (df$Focal_Response_Category[df$Video_ID == video_id & df$Focal_Fish_ID == focal_fish_id][1] %in% c("SRNL", "SRL")) {
      response_list[length(frame_list)] <- 1
    }

    for (t in seq_along(frame_list)) {
      current_frame <- frame_list[t]
      valid <- which(current_frame >= response_frames)
      if (length(valid) > 0) {
        times <- (current_frame - response_frames[valid])
        d <- distances[valid] / 1000
        if (include_speed == 1) {
          s <- speeds[valid] / 1000
          D_i_t <- sum(exp(-times / tau) * s^c / (d^2))
        } else {
          D_i_t <- sum(exp(-times / tau) / (d^2))
        }
        D_i_t_list[t] <- ifelse(is.finite(D_i_t), D_i_t, 0)
      } else {
        D_i_t_list[t] <- 0
      }
    }

    D_i_ts[[i]]           <- D_i_t_list
    responses[[i]]        <- response_list
    local_densities[[i]]  <- rep(neighbour_proximity, length(frame_list))
    focal_speeds[[i]]     <- rep(focal_speed, length(frame_list))
    focal_coral_ms[[i]]   <- rep(focal_coral_m, length(frame_list))
  }

  list(
    D_i_t          = unlist(D_i_ts),
    response       = unlist(responses),
    neighbour_proximity  = unlist(local_densities),
    focal_speed    = unlist(focal_speeds),
    focal_coral_m  = unlist(focal_coral_ms)
  )
}

# ---------- Objective (log coral when included) ----------
objective_function <- function(params, df, s_times, settings) {
  idx <- 1
  if (settings$include_time == 1) { tau <- max(1e-6, params[idx]); idx <- idx + 1 } else { tau <- 1e6 }
  if (settings$include_density == 1) { b0 <- params[idx]; idx <- idx + 1 } else { b0 <- 0 }
  b1 <- params[idx]; idx <- idx + 1
  if (settings$include_speed == 1) { c <- params[idx]; idx <- idx + 1 } else { c <- 1 }
  if (settings$include_coral == 1) { bc <- params[idx]; idx <- idx + 1 } else { bc <- 0 }

  res <- calculate_decision_and_response(
    df, s_times, tau = tau, b0 = b0, b1 = b1, c = c,
    detection_delay = settings$detection_delay,
    include_speed   = settings$include_speed,
    include_density = settings$include_density
  )

  D_i_t         <- pmax(res$D_i_t, 0)
  responses     <- res$response
  neighbour_proximity <- res$neighbour_proximity
  coral_m       <- res$focal_coral_m

  # Always log when coral is included; guard against zero
  coral_log <- log(pmax(coral_m, 1e-9))

  lambda <- 1e-5 + b1 * D_i_t
  if (settings$include_density == 1) lambda <- lambda * exp(b0 * neighbour_proximity)
  if (settings$include_coral   == 1) lambda <- lambda * exp(bc * coral_log)

  lambda <- pmax(lambda, 0)
  p <- 1 - exp(-lambda)

  eps <- 1e-12
  p <- pmin(pmax(p, eps), 1 - eps)

  -2 * sum(responses * log(p) + (1 - responses) * log(1 - p))
}

# ---------- Model grid: all 16 combinations (coral is log when included) ----------
grid <- expand.grid(
  include_time    = c(0, 1),
  include_density = c(0, 1),
  include_speed   = c(0, 1),
  include_coral   = c(0, 1)
)

name_model <- function(row) {
  parts <- c("Di")
  if (row$include_time == 1)    parts <- c(parts, "T")
  if (row$include_density == 1) parts <- c(parts, "De")
  if (row$include_speed == 1)   parts <- c(parts, "S")
  name <- paste0(parts, collapse = "")
  if (row$include_coral == 1) name <- paste0(name, "C")
  name
}

model_variants <- setNames(
  lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ])),
  sapply(split(grid, seq_len(nrow(grid))), name_model)
)

# ---------- Null deviance ----------
tmp <- calculate_decision_and_response(
  df, s_times,
  tau = 10, b0 = 0, b1 = 0, c = 0,
  detection_delay = 5, include_speed = 1, include_density = 1
)
resp_series   <- tmp$response
p0            <- min(max(mean(resp_series), 1e-12), 1 - 1e-12)
null_deviance <- -2 * sum(resp_series * log(p0) + (1 - resp_series) * log(1 - p0))

# ---------- Fit loop ----------
results <- data.frame()

for (model_name in names(model_variants)) {
  cat("\n--- Running model:", model_name, "---\n")
  settings <- model_variants[[model_name]]
  settings$detection_delay <- 5

  # Initial params in the SAME order objective_function unpacks them
  init <- numeric(0)
  if (settings$include_time == 1)     init <- c(init, tau = 12.3)
  if (settings$include_density == 1)  init <- c(init, b0  = -5.2e-3)
  init <- c(init, b1 = 3.61e-3)
  if (settings$include_speed == 1)    init <- c(init, c   = 1.78)
  if (settings$include_coral == 1)    init <- c(init, bc  = 0.0)

  # name parameters so we can read them by name later
  names(init) <- c(
    if (settings$include_time == 1)     "tau",
    if (settings$include_density == 1)  "b0",
    "b1",
    if (settings$include_speed == 1)    "c",
    if (settings$include_coral == 1)    "bc"
  )

  fit <- tryCatch({
    optim_result <- optim(
      par = init,
      fn = objective_function,
      method = "Nelder-Mead",
      df = df,
      s_times = s_times,
      settings = settings,
      control = list(maxit = 2000)
    )

    k <- (settings$include_time == 1) +
         (settings$include_density == 1) +
         1 +  # b1
         (settings$include_speed == 1) +
         (settings$include_coral == 1)   # bc

    dev <- optim_result$value
    aic <- dev + 2 * k
    r2  <- 1 - (dev / null_deviance)

    pars <- optim_result$par
    getp <- function(nm) if (nm %in% names(pars)) as.numeric(pars[[nm]]) else NA_real_

    data.frame(
      model = model_name,
      tau   = getp("tau"),
      b0    = getp("b0"),
      b1    = getp("b1"),
      c     = getp("c"),
      bc    = getp("bc"),
      deviance = dev,
      AIC      = aic,
      R2       = r2,
      Include_Time    = settings$include_time,
      Include_Density = settings$include_density,
      Include_Speed   = settings$include_speed,
      Include_Coral   = settings$include_coral,   # coral is log when included
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning(paste("Model", model_name, "failed:", e$message))
    data.frame(
      model = model_name, tau = NA, b0 = NA, b1 = NA, c = NA, bc = NA,
      deviance = NA, AIC = NA, R2 = NA,
      Include_Time    = settings$include_time,
      Include_Density = settings$include_density,
      Include_Speed   = settings$include_speed,
      Include_Coral   = settings$include_coral,
      stringsAsFactors = FALSE
    )
  })

  results <- bind_rows(results, fit)
}

# ---------- Save ----------
outfile <- file.path(OUT_PATH, "model_comparisons_time_dependent.csv")
write.csv(results, outfile, row.names = FALSE)
message("Wrote model_comparisons_time_dependent.csv (full factorial; coral as log): ", outfile)

print(results)
