# 01_io_and_helpers.R

numify <- function(x) suppressWarnings(as.numeric(x))

fmt_num <- function(x, digits = 3) {
  s <- format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE)
  gsub("\\.", "p", s)
}

log_msg <- function(..., .ts = TRUE) {
  pref <- if (.ts) sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")) else ""
  cat(pref, sprintf(...), "\n")
}

make_run_tag <- function(EXP) {
  paste(
    format(Sys.time(), "%Y%m%d-%H%M%S"),
    paste0("NFR", paste(EXP$N_FRs, collapse = "-")),
    paste0("FRspd", EXP$FR_speed_mode),
    paste0("dec", fmt_num(EXP$speed_decay, 4)),
    paste0("tmax", EXP$t_max),
    sep = "_"
  )
}

seed_for <- function(A, D, rep_id) {
  Ai <- as.integer(round(A * 1000))
  Di <- as.integer(round(D * 1000))
  as.integer((Ai * 1e8 + Di * 1e4 + rep_id * 13L) %% .Machine$integer.max)
}

make_safe_folder <- function(out_base, run_tag) {
  p <- file.path(out_base, paste0("experimental_run_", run_tag))
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  p
}

make_partial_filename <- function(folder, model_name, A, N, nfr, EXP, run_tag, spd_lab = NA_character_) {
  parts <- c(
    paste0("Model", model_name),
    paste0("A", fmt_num(A, 1)),
    paste0("N", N),
    paste0("NFR", nfr),
    paste0("FRspd", EXP$FR_speed_mode),
    if (EXP$FR_speed_mode == "fixed")  paste0("FRv", fmt_num(EXP$FR_speed_fixed_value, 3)) else NULL,
    if (EXP$FR_speed_mode == "binary") paste0("FR", spd_lab) else NULL,
    paste0("dec", fmt_num(EXP$speed_decay, 4)),
    run_tag
  )
  file.path(folder, paste0(paste(parts, collapse = "_"), ".csv"))
}
