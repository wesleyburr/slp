#  slp
#
#  Wrapper function for creation of Discrete Prolate Spheroidal (Slepian) Sequence
#  Regression Smoothers. 
#


slp <- function(x, W = NA, K = NA, deltat = 1, naive = FALSE) {

  # logical checks (note: W is assumed to be in the same units as deltat, i.e., days)
  stopifnot(is.numeric(x), is.numeric(W), W > 0, W < 1/(2*deltat), deltat %in% c(1, 6))

  # x should be a contiguous time array; if not, convert to one, back-convert at the end
  if(is.na(max(x[-1] - x[-length(x)])) | max(x[-1] - x[-length(x)]) > 1) {
    warning("slp: Input time array is not contiguous.")
    minT <- min(x, na.rm = TRUE); maxT <- max(x, na.rm = TRUE)
    stopifnot(is.integer(minT), is.integer(maxT))
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
    if(!is.integer(K)) { 
      K <- floor(K)
      warning(paste("slp: K choice not integer. Truncated to K = ", K, sep = ""))
    }
    W <- as.numeric((K + 1)) / as.numeric((2 * length(wx)))
  }

  # start by generating baseline Slepian basis vectors
  v <- dpss(n = N, nw = N * W, k = K)
  if(naive) {
  } else {

    # Equations follow from Thomson (2001), see help file for full citation
    alpha <- t(rep(1, N)) %*% v %*% t(v) %*% rep(1, N) 
    R <- rep(1/sqrt(alpha), N)
    U <- t(v) %*% R
    sRaw <- (R %*% t(U) + v %*% (diag(K) - U %*% t(U))) %*% t(v) 
      
    # reform into positive-definite (idempotent + symmetric) basis vector form via SVD
    Xmat <- svd(sRaw, nu = K, nv = 0)$u  # very, very slow
  }

  # need to convert back to original time stamps
  # need to make similar to ns() object

  if(naive) {
    return(v) 
  } else {
    return(Xmat) 
  }
}

