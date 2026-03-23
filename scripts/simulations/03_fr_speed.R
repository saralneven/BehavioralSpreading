# 03_fr_speed.R
# FR speed utilities: fit distribution, draw speeds, define fixed/binary speeds

fit_fr_speed_lognormal <- function(path) {
  dt <- data.table::fread(path, showProgress = FALSE)

  v <- dt$Speed
  v <- v[!is.na(v) & is.finite(v) & v > 0]
  v <- v / 1000  # mm/s -> m/s

  list(mu = mean(log(v)), sigma = sd(log(v)))
}

draw_speed0 <- function(n, mu, sigma) {
  rlnorm(n, meanlog = mu, sdlog = sigma)
}

get_weak_strong_FR_speeds <- function(path, probs = c(0.25, 0.75)) {
  dt <- data.table::fread(path, showProgress = FALSE)

  v <- dt$Speed
  v <- v[!is.na(v) & is.finite(v) & v > 0]
  v <- v / 1000  # mm/s -> m/s

  qs <- as.numeric(quantile(v, probs = probs, na.rm = TRUE))

  data.frame(
    label = c("weak", "strong"),
    quantile = probs,
    FR_speed_m_per_s = qs
  )
}

# 03_fr_speed.R

make_speed0 <- function(N, fr_index, mu, sigma, EXP,
                        fr_speed_label = NA_character_,
                        fr_speed_value = NA_real_) {

  mode <- EXP$FR_speed_mode

  if (mode == "lognormal") {
    speed0 <- draw_speed0(N, mu, sigma)
    fr_speed0 <- speed0[fr_index]

  } else if (mode == "fixed") {
    fr_speed0 <- as.numeric(EXP$FR_speed_fixed_value)
    speed0 <- rep(fr_speed0, N)

  } else if (mode == "binary") {
    if (is.na(fr_speed_label)) stop("binary mode requires fr_speed_label")
    i <- match(fr_speed_label, EXP$FR_speed_labels)
    if (is.na(i)) stop("Unknown fr_speed_label: ", fr_speed_label)
    fr_speed0 <- as.numeric(EXP$FR_speed_binary_values[i])
    speed0 <- rep(fr_speed0, N)

  } else if (mode == "grid") {
    if (!is.finite(fr_speed_value)) stop("grid mode requires fr_speed_value (numeric)")
    fr_speed0 <- as.numeric(fr_speed_value)
    speed0 <- rep(fr_speed0, N)

  } else {
    stop("Unknown FR_speed_mode: ", mode)
  }

  # your rule: everyone has FR speed
  speed0[] <- fr_speed0
  list(speed0 = speed0, fr_speed0 = fr_speed0)
}

get_fr_speed_conditions <- function(EXP) {
  if (EXP$FR_speed_mode == "lognormal") return(list(type="label", values=""))
  if (EXP$FR_speed_mode == "fixed")     return(list(type="label", values="fixed"))
  if (EXP$FR_speed_mode == "binary")    return(list(type="label", values=EXP$FR_speed_labels))
  if (EXP$FR_speed_mode == "grid")      return(list(type="value", values=as.numeric(EXP$FR_speed_grid_values)))
  stop("Unknown FR_speed_mode: ", EXP$FR_speed_mode)
}