%%	OK

\name{locCteWeights}
\alias{locCteWeightsC}
\alias{locLinWeightsC}
\alias{locWeightsEvalC}
\alias{locWeightsEval}
\alias{locPolWeights}
\title{ Local Polynomial Weights }
\description{
  Local Constant and local Linear estimator with weight.
}
\usage{
  locCteWeightsC(x, xeval, bw, kernel, weig = rep(1, length(x)))
  locLinWeightsC(x, xeval, bw, kernel, weig = rep(1, length(x)))
  locPolWeights(x, xeval, deg, bw, kernel, weig = rep(1, length(x)))
  locWeightsEval(lpweig, y)
  locWeightsEvalC(lpweig, y)
}
\arguments{
	\item{x}{ x covariate data values. }
	\item{y}{ y response data values. }
	\item{xeval}{ Vector with evaluation points. }
	\item{bw}{ Smoothing parameter, bandwidth. }
	\item{deg}{  Local polynomial estimation degree (\eqn{p}). }
	\item{kernel}{ Kernel used to perform the estimation, see \code{\link{Kernels}}}
	\item{weig}{ Vector of weights for observations. }
	\item{lpweig}{ Local polynomial weights \eqn{(X^TWX)^{-1}X^TW} evaluated at \code{xeval} matrix. }
}
\details{
	\code{locCteWeightsC} and \code{locLinWeightsC} computes local
	constant and local linear weights, say any of the entries of
	the vector \eqn{(X^TWX)^{-1}X^TW} for \eqn{p=0} and \eqn{p=1} resp.
	\code{locWeightsEvalC} and \code{locWeightsEval} computes local
	the estimator for a given vector of responses \code{y}

% \code{locWeightsEvalC} and \code{locWeightsEval} returns the
%	vector with the estimation. It performs the matrix multiplication
%	between \code{locWeig} and \code{y} to obtain the estimation at
%	given \code{xeval} points.
}
\value{
	\code{locCteWeightsC} and \code{locLinWeightsC} returns a list
	with two components
	\item{den}{Estimation of \eqn{(n*h*f(x))^{p+1}} being \eqn{h} the
        bandwidth \code{bw}.}
	\item{locWeig}{\eqn{(X^TWX)^{-1}X^TW} evaluated at \code{xeval} Matrix.}
}
\references{
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman & Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).
}
\author{ Jorge Luis Ojeda Cabrera.}
\seealso{
\code{\link{Kernels}}, \link{locpol}.
}
\examples{
	size <- 200
	sigma <- 0.25
	deg <- 1
	kernel <- EpaK
	bw <- .25
	xeval <- 0:100/100
	regFun <- function(x) x^3
	x <- runif(size)
	y <- regFun(x) + rnorm(x, sd = sigma)
	d <- data.frame(x, y)
	lcw <- locCteWeightsC(d$x, xeval, bw, kernel)$locWeig
	lce <- locWeightsEval(lcw, y)
	lceB <- locCteSmootherC(d$x, d$y, xeval, bw, kernel)$beta0
	mean((lce-lceB)^2)
    llw <- locLinWeightsC(d$x, xeval, bw, kernel)$locWeig
	lle <- locWeightsEval(llw, y)
	lleB <- locLinSmootherC(d$x, d$y, xeval, bw, kernel)$beta0
	mean((lle-lleB)^2)
}
\keyword{ nonparametric }
\keyword{ smooth }
