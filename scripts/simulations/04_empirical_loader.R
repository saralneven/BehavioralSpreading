# 04_empirical_loader.R
# Load + sample empirical 3D networks from OverallPositions.csv
# Returns: loader object with:
#   - real_dt_unique (data.table)
#   - available_real_sizes (int)
#   - videos_by_size (list)
#   - sample_real_positions_3d(group_size, Area, rep_id)

load_empirical_positions_3d <- function(real_positions_csv, EXP) {
  stopifnot(file.exists(real_positions_csv))

  scale_factor <- 1 / 1000  # mm -> m

  dt_raw <- data.table::fread(real_positions_csv, showProgress = FALSE)

  keep_cols <- intersect(names(dt_raw), c("Video_ID", "Fish_ID", "X", "Y", "Z"))
  if (!all(c("Video_ID", "X", "Y", "Z") %in% keep_cols)) {
    stop("OverallPositions.csv must include at least Video_ID, X, Y, Z (and optionally Fish_ID).")
  }

  dt <- dt_raw[, ..keep_cols]

  dt[, `:=`(
    Video_ID = as.integer(Video_ID),
    Fish_ID  = if ("Fish_ID" %in% names(dt)) as.integer(Fish_ID) else NA_integer_,
    X = as.numeric(X) * scale_factor,
    Y = as.numeric(Y) * scale_factor,
    Z = as.numeric(Z) * scale_factor
  )]

  # drop missing coords
  dt <- dt[is.finite(X) & is.finite(Y) & is.finite(Z)]

  # fallback Fish_ID if missing
  if (!("Fish_ID" %in% names(dt)) || all(is.na(dt$Fish_ID))) {
    warning("Fish_ID missing or all NA in OverallPositions.csv; creating fallback Fish_ID per row.")
    dt[, Fish_ID := seq_len(.N)]
  }

  # robust default snapshot strategy: one row per (Video_ID, Fish_ID)
  dt_unique <- dt[, .SD[1], by = .(Video_ID, Fish_ID)]

  video_sizes <- dt_unique[, .(Group_Size = uniqueN(Fish_ID)), by = Video_ID]
  videos_by_size_raw <- split(video_sizes$Video_ID, video_sizes$Group_Size)

  # drop NA names + empties
  videos_by_size <- videos_by_size_raw
  videos_by_size <- videos_by_size[!is.na(names(videos_by_size))]
  videos_by_size <- videos_by_size[vapply(videos_by_size, length, integer(1)) > 0]

  available_sizes <- sort(as.integer(names(videos_by_size)))
  if (length(available_sizes) == 0) stop("No valid group sizes found in OverallPositions.csv after filtering.")

  nearest_nonempty_size <- function(gs) {
    gs <- as.integer(gs)
    available_sizes[which.min(abs(available_sizes - gs))]
  }

  sample_real_positions_3d <- function(group_size, Area, rep_id) {
    requested <- as.integer(group_size)

    gs <- requested
    if (!(as.character(gs) %in% names(videos_by_size))) {
      if (isTRUE(EXP$empirical_use_nearest_size)) {
        gs <- nearest_nonempty_size(gs)
      } else {
        stop("Requested group size not available in empirical networks: ", requested)
      }
    }

    vids <- videos_by_size[[as.character(gs)]]
    if (is.null(vids) || length(vids) == 0) {
      stop("No empirical videos available for mapped size: ", gs)
    }

    # try multiple times in case a sampled video has too few usable fish
    for (attempt in seq_len(as.integer(EXP$empirical_max_tries))) {
      vid <- sample(vids, 1)

      pos_dt <- dt_unique[Video_ID == vid, .(Fish_ID, X, Y, Z)]
      pos_dt <- pos_dt[!is.na(Fish_ID)]
      pos_dt <- pos_dt[!duplicated(Fish_ID)]

      if (nrow(pos_dt) < as.integer(EXP$empirical_min_fish)) next

      # if more fish than target size, subsample down to gs
      if (nrow(pos_dt) > gs) {
        take_ids <- sample(pos_dt$Fish_ID, gs, replace = FALSE)
        pos_dt <- pos_dt[Fish_ID %in% take_ids]
      }

      return(list(
        positions            = as.data.frame(pos_dt[, .(X, Y, Z)]),
        Video_ID             = as.integer(vid),
        Group_Size_Used      = as.integer(nrow(pos_dt)),
        Group_Size_Target    = as.integer(gs),
        Group_Size_Requested = as.integer(requested)
      ))
    }

    stop("Failed to sample a usable empirical video for requested size ", requested, " (mapped to ", gs, ").")
  }

  list(
    real_dt_unique         = dt_unique,
    videos_by_size         = videos_by_size,
    available_real_sizes   = available_sizes,
    sample_real_positions_3d = sample_real_positions_3d
  )
}
