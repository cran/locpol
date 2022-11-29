%%	OK


\name{PRDenEstC}
\alias{PRDenEstC}
\title{ Parzen--Rosenblatt denstiy estimator. }
\description{
  Parzen--Rosenblat univariate density estimator.
}
\usage{
PRDenEstC(x, xeval, bw, kernel, weig = rep(1, length(x)))
}
\arguments{
  \item{x}{ vector with data points. }
  \item{xeval}{ Vector of evaluation points. }
  \item{bw}{ Smoothing parameter, bandwidth. }
  \item{kernel}{ Kernel used to perform the estimation, see \code{\link{Kernels}}}
  \item{weig}{ Vector of weights for observations. }
}
\details{
  Simple Parzen--Rosenblat univariate density estimation, computed using
  definition.
}
\value{
  Returns an \code{(x,den)} data frame.

  \item{x }{Evaluation points.}
  \item{den }{Density at each \code{x} point.}
}
\references{
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman & Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).
}
\author{
  Jorge Luis Ojeda Cabrera.
}
\seealso{
  \code{\link{density}}, that uses FT to compute a kernel density
estimator,  \code{\link{bkde}} from package \pkg{KernSmooth} for a
binned version, and \code{\link{bw.nrd0}}, \code{\link{dpik}}, \code{\link{denCVBwSelC}} for bandwidth selection.
}
\examples{
	N <- 100
	x <-  runif(N)
	xeval <- 0:10/10
	b0.125 <- PRDenEstC(x, xeval, 0.125, EpaK)
	b0.05 <- PRDenEstC(x, xeval, 0.05, EpaK)
	cbind(x = xeval, fx = 1, b0.125 = b0.125$den, b0.05 = b0.05$den)
}
\keyword{ nonparametric }
\keyword{ smooth }
