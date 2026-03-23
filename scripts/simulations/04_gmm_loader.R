# 04_gmm_loader.R

load_gmms_by_size <- function(gmm_csv) {
  gmm_params_df <- as.data.frame(data.table::fread(gmm_csv))

  scale_factor <- 1/1000  # mm -> m

  gmm_by_size <- lapply(
    split(gmm_params_df, gmm_params_df$Group_Size),
    function(df_size) {

      comps <- split(df_size, df_size$Component)

      weights <- sapply(comps, function(x) x$Weight[1])
      mus <- lapply(comps, function(x) c(x$Mean_X[1], x$Mean_Y[1], x$Mean_Z[1]) * scale_factor)

      Sigmas <- lapply(comps, function(x) {
        S_mm <- matrix(c(
          x$Cov_xx[1], x$Cov_xy[1], x$Cov_xz[1],
          x$Cov_xy[1], x$Cov_yy[1], x$Cov_yz[1],
          x$Cov_xz[1], x$Cov_yz[1], x$Cov_zz[1]
        ), nrow = 3, byrow = TRUE)

        S_m <- S_mm * (scale_factor^2)
        S_m + diag(1e-9, 3)
      })

      list(weights = as.numeric(weights), mus = mus, Sigmas = Sigmas)
    }
  )

  available_sizes <- sort(as.integer(names(gmm_by_size)))

  get_gmm_for_size <- function(group_size) {
    gs <- as.numeric(group_size)
    nm <- as.character(gs)

    if (nm %in% names(gmm_by_size)) return(gmm_by_size[[nm]])

    if (gs <= min(available_sizes)) return(gmm_by_size[[as.character(min(available_sizes))]])
    if (gs >= max(available_sizes)) return(gmm_by_size[[as.character(max(available_sizes))]])

    lower_sizes <- available_sizes[available_sizes < gs]
    upper_sizes <- available_sizes[available_sizes > gs]
    s_lo <- max(lower_sizes)
    s_hi <- min(upper_sizes)

    g_lo <- gmm_by_size[[as.character(s_lo)]]
    g_hi <- gmm_by_size[[as.character(s_hi)]]
    alpha <- (gs - s_lo) / (s_hi - s_lo)

    w <- (1 - alpha) * g_lo$weights + alpha * g_hi$weights
    w[w < 0] <- 0
    if (sum(w) == 0) w <- rep(1 / length(w), length(w)) else w <- w / sum(w)

    mus <- Map(function(m0, m1) (1 - alpha) * m0 + alpha * m1, g_lo$mus, g_hi$mus)

    Sigmas <- Map(function(S0, S1) {
      S <- (1 - alpha) * S0 + alpha * S1
      S + diag(1e-9, 3)
    }, g_lo$Sigmas, g_hi$Sigmas)

    list(weights = w, mus = mus, Sigmas = Sigmas)
  }

  sample_from_gmm_3d <- function(N, gmm) {
    K <- length(gmm$weights)
    if (K == 0) stop("GMM has zero components.")
    w <- gmm$weights / sum(gmm$weights)

    comp_ids <- sample.int(K, size = N, replace = TRUE, prob = w)

    out <- matrix(NA_real_, nrow = N, ncol = 3L)
    for (k in seq_len(K)) {
      idx <- which(comp_ids == k)
      if (!length(idx)) next
      mu <- gmm$mus[[k]]
      S  <- gmm$Sigmas[[k]]
      out[idx, ] <- MASS::mvrnorm(n = length(idx), mu = mu, Sigma = S)
    }
    data.frame(X = out[,1], Y = out[,2], Z = out[,3])
  }

  list(
    gmm_by_size = gmm_by_size,
    available_sizes = available_sizes,
    get_gmm_for_size = get_gmm_for_size,
    sample_from_gmm_3d = sample_from_gmm_3d
  )
}
