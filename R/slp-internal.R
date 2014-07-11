#  slp
#
#  Wrapper function for creation of Discrete Prolate Spheroidal (Slepian) Sequence
#  Regression Smoothers. 
#

slp <- function(x, W = NA, K = NA, deltat = 1, naive = FALSE, intercept = FALSE, forceC = TRUE) {

  # logical checks (note: W is assumed to be in the same units as deltat, i.e., days)
  stopifnot(is.numeric(x), (is.na(W) | is.numeric(W) & W > 0 & W < 1/(2*deltat)), deltat %in% c(1, 6),
            (is.numeric(K) & K <= length(x) | is.na(K)))

  namesx <- names(x)
  x <- as.vector(x)

  # x should be a contiguous time array; if not, convert to one, back-convert at the end
  if(is.na(max(x[-1L] - x[-length(x)])) | max(x[-1L] - x[-length(x)]) > 1) {
    warning("slp: Input time array is not contiguous.")
    minT <- min(x, na.rm = TRUE); maxT <- max(x, na.rm = TRUE)
    stopifnot(is.numeric(minT), is.numeric(maxT))
    wx <- seq(minT, maxT, deltat)

    if(deltat == 6) {
      wx <- seq(minT, maxT, deltat)
      if(wx[length(wx)] != maxT) { stop("slp: Input time array not properly time aligned to 6-day samples per deltat = 6.") }
    }

          cat(" Input Time Array: \n")
    cat(paste(" *        samples: ", length(x), "\n", sep = ""))
    cat(paste(" *      max - min: ", maxT - minT + 1, "\n", sep = ""))
    cat(paste(" *         deltat: ", deltat, "\n", sep = ""))
    cat(paste(" * number missing: ", length(which(is.na(x))), "\n", sep =""))
  } else {
    wx <- x 
  }

  N <- length(wx)

  # user must set _one_ of K or W
  if(is.na(W) & is.na(K)) { stop("Must set one of K or W for family selection.") }
  
  # if user does not set K, default to 2NW-1
  if(is.na(K)) { K <- floor(2 * length(wx) * W - 1) }
  # otherwise, user has selected K, default W to appropriate choice by above
  if(is.na(W)) { 
    if(K != round(K)) { 
      K <- floor(K)
      warning(paste("slp: K choice not integer. Truncated to K = ", K, sep = ""))
    }
    W <- as.numeric((K + 1)) / as.numeric((2 * length(wx)))
  }

  # Logical check: if N, K and W match one of the saved objects, simply
  # return that object; otherwise, run through the generation process
  if(checkSaved(N, W, K) & !forceC) {
    Wn <- round(W * 365.2425)
    data(paste0("basis_N_", N[j], "_W_", W[k], "_K_", K, ".RData"))

    if(!intercept) { basis <- basis[, -1] }
    
    # need to convert back to original time array -- the NAs in the original
    if(any(is.na(x))) {
      basis[which(!(wx %in% x[!is.na(x)])), ] <- rep(NA, ncol(basis))
    }

  } else {

      # start by generating baseline Slepian basis vectors
      v <- .dpss(n = N, nw = N * W, k = K)
      if(naive) {
        basis <- v
      } else {
    
        # Equations follow from Thomson (2001), see help file for full citation
        alpha <- t(rep(1, N)) %*% v %*% t(v) %*% rep(1, N) 
        R <- rep(1/sqrt(alpha), N)
        U <- t(v) %*% R
        sRaw <- (R %*% t(U) + v %*% (diag(K) - U %*% t(U))) %*% t(v) 
         
        # Find an orthonormal matrix A such that A %*% t(A) \approx sRaw -- impossible to get equality
        basis <- svd(sRaw, nu = K, nv = K)$u  # very, very slow -- actually computes all N :(
        # *** find a speed-up for computing _only_ the K highest eval/evect pairs
        # ** doesn't actually give what we want, because U != t(V)
    
        # idempotent ==> eigenvalues are 0 or 1
        # rank = K ==> K 1 eval, N-K 0 eval
        # want SVD's U or equivalent -- the top K eigenvectors corresponding to the K 1 evals
      }
    
      dimnames(basis) <- list(names(wx), 1L:ncol(basis))
    
      # two options: either capture all constant value (in which case the smooth conflicts with 
      # the intercept term, and the first column is NA'd; so use with "-1") -- "intercept = TRUE";
      # or ignore all constant value (in which case the first column should be dropped) -- "intercept = FALSE"
      if(!intercept) { basis <- basis[, -1] }
    
      # need to convert back to original time array -- the NAs in the original
      if(any(is.na(x))) {
        basis[which(!(wx %in% x[!is.na(x)])), ] <- rep(NA, ncol(basis))
      }

      a <- list(K = K, W = W, N = N, naive = naive)
      attributes(basis) <- c(attributes(basis), a)
      class(basis) <- c("slp", "basis", "matrix")
  }
  basis
}

