# 02_load_model_fits.R

load_model_fits <- function(fits_csv, EXP) {
  model_fits <- read.csv(fits_csv, stringsAsFactors = FALSE) %>%
    mutate(
      model           = as.character(model),
      tau             = numify(tau),
      b0              = numify(b0),
      b1              = numify(b1),
      c               = numify(c),
      Include_Speed   = as.integer(Include_Speed),
      Include_Density = as.integer(Include_Density),
      Include_Time    = as.integer(Include_Time)
    )

  params_by_model <- model_fits %>%
    rowwise() %>%
    mutate(
      params = list(list(
        tau = if (Include_Time == 1) tau else 1e6,
        b0  = if (Include_Density == 1) b0 else 0,
        b1  = b1,
        c   = if (Include_Speed == 1) c else 0,
        detection_delay      = EXP$detection_delay,
        speed_decay          = EXP$speed_decay,
        max_interaction_dist = EXP$max_interaction_dist
      ))
    ) %>%
    ungroup() %>%
    select(model, Include_Speed, Include_Density, Include_Time, params)

  get_model_spec <- function(model_name) {
    row <- params_by_model %>% filter(model == model_name)
    if (nrow(row) == 0) stop(sprintf("No fitted row found for model: %s", model_name))
    list(
      include_speed   = row$Include_Speed[1],
      include_density = row$Include_Density[1],
      include_time    = row$Include_Time[1],
      params          = row$params[[1]]
    )
  }

  list(params_by_model = params_by_model, get_model_spec = get_model_spec)
}
