# 05_geometry.R

make_distance_weights <- function(positions,
                                  max_interaction_dist = Inf,  # set Inf to match fitted (unlimited)
                                  min_sep = 0.02               # meters; pick based on body length
                                  ) {

  D_full <- as.matrix(dist(positions))
  p <- nrow(D_full)

  if (p == 0) {
    W <- Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0),
                              dims=c(0,0), symmetric=FALSE)
    return(list(W=W, LD=numeric(0), D=D_full))
  }

  diag(D_full) <- NA_real_

  # --- biological realism: prevent 1/d^2 blow-ups from near-overlaps ---
  D_full[is.finite(D_full) & D_full < min_sep] <- min_sep

  # --- Local density: ALL neighbours (matches your data calculation idea) ---
  LD <- apply(D_full, 1, function(row) sum(1 / (row^2), na.rm = TRUE))
  LD <- as.numeric(LD)

  # --- Social influence weights W: either unlimited (Inf) or cut at max_interaction_dist ---
  D_cut <- D_full
  if (is.finite(max_interaction_dist)) {
    D_cut[D_cut > max_interaction_dist] <- NA_real_
  }

  nz <- which(is.finite(D_cut), arr.ind = TRUE)

  if (nrow(nz) == 0) {
    W <- Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0),
                              dims=c(p,p), symmetric=FALSE)
    return(list(W=W, LD=LD, D=D_full))
  }

  vals <- 1 / (D_cut[nz]^2)
  W <- Matrix::sparseMatrix(i = nz[,1], j = nz[,2], x = vals,
                            dims = c(p, p), symmetric = FALSE)

  list(W = W, LD = LD, D = D_full)
}
