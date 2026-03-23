# 10_calculate_Rt.R
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

library(here)
BASE_PATH <- here::here("outputs", "simulations", "GMM_3D", "diagnostics")

IN_DIR  <- file.path(BASE_PATH, "TRACEINPUTS_20260302-154039_DiTDeS_NFR1_gmm_dec0")  # update with actual trace folder name
OUT_DIR <- IN_DIR

SUMMARY_FILE <- file.path(IN_DIR, "SUMMARY_ALL.csv")
CONTRIB_FILE <- file.path(IN_DIR, "CONTRIB_TOPK_ALL.csv")

OUT_RT <- file.path(OUT_DIR, "Rt_long.csv")
OUT_RI <- file.path(OUT_DIR, "Ri_long.csv")

summ   <- fread(SUMMARY_FILE)
contrib <- fread(CONTRIB_FILE)

# Key columns that exist in your files
KEY <- c("Rep_ID","Area","N","N_req","Model","Network","N_FRs","FR_ID","FR_label","FR_speed")

# --- Fix responded: your file uses 0/1 ---
summ$responded <- as.integer(summ$responded)
summ$responded <- summ$responded == 1

# --- Startle times for responders (targets) ---
startle_map <- summ %>%
  filter(responded, is.finite(startle_t)) %>%
  select(all_of(KEY), Fish_ID, startle_t) %>%
  distinct()

# --- Add FR row at t=0 ---
fr_rows <- startle_map %>%
  select(all_of(KEY)) %>% distinct() %>%
  mutate(Fish_ID = as.integer(FR_ID), startle_t = 0)

startle_map <- bind_rows(startle_map, fr_rows) %>%
  distinct(across(all_of(KEY)), Fish_ID, .keep_all = TRUE)

# Fish that startled
startled_fish <- startle_map %>%
  filter(is.finite(startle_t)) %>%
  select(all_of(KEY), Fish_ID, startle_t)

# --- Compute p_ji from contrib weights per recipient event ---
contrib2 <- contrib %>%
  mutate(
    contrib   = as.numeric(contrib),
    Fish_ID   = as.integer(Fish_ID),
    neigh_id  = as.integer(neigh_id),
    startle_t = as.numeric(startle_t)
  ) %>%
  filter(is.finite(startle_t), is.finite(contrib), contrib > 0) %>%
  group_by(across(all_of(KEY)), Fish_ID, startle_t) %>%
  mutate(
    w_sum = sum(contrib, na.rm = TRUE),
    p_ji  = contrib / w_sum
  ) %>%
  ungroup()

# --- R_i: sum of fractional attributions ---
Ri_raw <- contrib2 %>%
  group_by(across(all_of(KEY)), neigh_id) %>%
  summarise(R_i = sum(p_ji, na.rm = TRUE), .groups = "drop") %>%
  rename(Fish_ID = neigh_id)

Ri_full <- startled_fish %>%
  left_join(Ri_raw, by = c(KEY, "Fish_ID")) %>%
  mutate(R_i = ifelse(is.na(R_i), 0, R_i))

# --- R(t): mean R_i for fish that startled at time t ---
Rt_long <- Ri_full %>%
  group_by(across(all_of(KEY)), time = startle_t) %>%
  summarise(
    R      = mean(R_i, na.rm = TRUE),
    n_fish = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(Rep_ID, Model, Network, N_req, FR_label, time)

fwrite(Rt_long, OUT_RT)
fwrite(Ri_full, OUT_RI)

cat("Wrote:\n", OUT_RT, "\n", OUT_RI, "\n")
cat("Rt_long rows:", nrow(Rt_long), "\n")