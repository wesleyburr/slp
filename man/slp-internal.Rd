\name{slp-internal}
\alias{slp}
\title{Generate a Basis Matrix for Discrete Prolate Spheroidal (Slepian) Sequences}
\description{
  Generate the basis matrix for a particular \code{N, W} Slepian sequence
  family member, with the additional property that the smoother passes constants
  without distortion. Quite slow execution due to the latter property.

  Based on \code{\link[splines]{ns}} for implementation with \code{\link[gam]{gam}}
  and on \code{\link[spsmooth]{smooth.construct.sp.smooth.spec}} for implementation with
  \code{\link[mgcv]{gam}}. 
}
\usage{
slp(x, W = NA, K = NA, deltat = 1, naive = FALSE)
}
\arguments{
  \item{x}{the predictor variable.  Missing values are allowed. Assumed to be contiguous;
      if not, then converted to a contiguous series to determine appropriate \code{N},
      \code{K} and \code{W}, then the basis vectors are back-converted at the termination
      of the routine.}
  \item{W}{the time bandwidth. Computed as the frequency domain analogue of the maximum period of
      interest for a time series-regression problem using ``smooth functions of time''. For example,
      a period choice of 2 months converts to 60 days and \code{W = 1/60} cycles per day. Alternatively,
      if the interest is in a period of \code{7} cycles per year, then \code{W = 7 / 365.2425 = 0.0192}
      cycles per day.}
  \item{K}{the number of basis vectors requested. If not provided, then \code{W} must be, and
      \code{K} is automatically computed as \code{floor(2 * N * W - 1)}. Approximately equivalent
      to \code{df} for \code{\link[splines]{ns}} with fixed dimension.}
  \item{deltat}{the time step for the input \code{x}. Restricted to \code{1} and \code{6} days
      for ease of logic checking, as these are the most traditional choices. Assumes that 
      \code{W}  is in the same units, and has no real impact beyond this, so it is trivial to 
      make \code{deltat} symbolically equal an arbitrary choice and convert \code{W} to match.}
  \item{naive}{a flag for returning either the naive (default) Slepian basis vectors \code{v}. If set
      to \code{FALSE}, returns the constant-reproducing variant.}
}
\details{
  \code{slp} is based on \code{\link{dpss}}, which generates a family of Discrete
  Prolate Spheroidal (Slepian) Sequences. These vectors are orthonormal, have alternating
  even/odd parity, and form the optimally concentrated basis set for the subspace of 
  \code{R^N} corresponding to the bandwidth \code{W}. Full details are given
  in Slepian (1978).  These basis functions have natural boundary conditions, and lack any form of 
  knot structure. 

  The \code{\link{dpss}} basis vectors are adapted to provide the additional
  useful property of passing constants perfectly. That is, the smoother matrix
  \code{S} formed from the returned rectangular matrix will reproduce constants
  at near round-off precision, i.e., \code{S \%/\% rep(1, N) = rep(1, N)}. 

  The primary use is in modeling formula to directly specify a Slepian time-based smoothing
  term in a model: see the examples.
}
\value{
  A matrix of dimension \code{length(x) * K} where either \code{K} was
  supplied, or \code{W} was supplied and \code{K} converted. Note that the
  basis vectors are computed on a contiguous grid based on \code{x}, and then
  back-converted to the time structure of \code{x}.

  Attributes are returned that correspond to the arguments to \code{ns},
  and explicitly give \code{K}, \code{W} etc for use by \code{predict.slp()}.
}
\seealso{
  \code{\link[splines]{ns}} 
}
\references{
}
\examples{
}
\keyword{smooth}
