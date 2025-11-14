%%	OK


\name{locpolSmoothers}
\alias{locCteSmootherC}
\alias{locLinSmootherC}
\alias{locCuadSmootherC}
\alias{locPolSmootherC}
\alias{looLocPolSmootherC}
\title{ Local Polynomial estimation. }
\description{
  Computes the local polynomial estimation of the regression function.
}
\usage{
locCteSmootherC(x, y, xeval, bw, kernel, weig = rep(1, length(y)))
locLinSmootherC(x, y, xeval, bw, kernel, weig = rep(1, length(y)))
locCuadSmootherC(x, y, xeval, bw, kernel, weig = rep(1, length(y)))
locPolSmootherC(x, y, xeval, bw, deg, kernel, DET = FALSE,
	weig = rep(1, length(y)))
looLocPolSmootherC(x, y, bw, deg, kernel, weig = rep(1, length(y)),
        DET = FALSE)
}
\arguments{
	\item{x}{ x covariate data values. }
	\item{y}{ y response data values. }
	\item{xeval}{ Vector of evaluation points. }
	\item{bw}{ Smoothing parameter, bandwidth. }
	\item{kernel}{ Kernel used to perform the estimation, see \code{\link{Kernels}}}
	\item{weig}{ Vector of weights for observations. }
	\item{deg}{ Local polynomial estimation degree (\eqn{p}). }
	\item{DET}{ Boolean to ask for the computation of the determinant if the matrix \eqn{X^TWX}. }
}
\details{
	All these function perform the estimation of the regression function
for different degrees. While \code{locCteSmootherC}, \code{locLinSmootherC},
and \code{locCuadSmootherC} uses direct computations for the degrees 0,1
and 2 respectively, \code{locPolSmootherC} implements a general method for any degree.
	Particularly useful can be \code{looLocPolSmootherC}(Leave one out) which computes the local polynomial estimator for any degree as \code{locPolSmootherC} does, but estimating \eqn{m(x_i)} without using \eqn{i}--th observation on the computation.
}
\value{
	A data frame whose components gives the evaluation points, the estimator
for the regression function \eqn{m(x)} and its derivatives at each point, and
the estimation of the marginal density for \code{x} to the \eqn{p+1} power.
These components are given by:
	\item{x}{Evaluation points.}
	\item{beta0, beta1, beta2,...}{Estimation of the \eqn{i}-th derivative
of the regression function (\eqn{m^{(i)}(x)}) for \eqn{i=0,1,...}. }
	\item{den}{Estimation of \eqn{(n*h*f(x))^{p+1}}, being \eqn{h} the
bandwidth \code{bw}.}
}
\references{
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman & Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).
}
\author{ Jorge Luis Ojeda Cabrera. }
\examples{
N <- 100
xeval <- 0:10/10
d <- data.frame(x = runif(N))
bw <- 0.125
fx <- xeval^2 - xeval + 1
##	Non random
d$y <- d$x^2 - d$x + 1
cuest <- locCuadSmootherC(d$x, d$y ,xeval, bw, EpaK)
lpest2 <- locPolSmootherC(d$x, d$y , xeval, bw, 2, EpaK)
print(cbind(x = xeval, fx, cuad0 = cuest$beta0,
lp0 = lpest2$beta0, cuad1 = cuest$beta1, lp1 = lpest2$beta1))
##	Random
d$y <- d$x^2 - d$x + 1 + rnorm(d$x, sd = 0.1)
cuest <- locCuadSmootherC(d$x,d$y , xeval, bw, EpaK)
lpest2 <- locPolSmootherC(d$x,d$y , xeval, bw, 2, EpaK)
lpest3 <- locPolSmootherC(d$x,d$y , xeval, bw, 3, EpaK)
cbind(x = xeval, fx, cuad0 = cuest$beta0, lp20 = lpest2$beta0,
lp30 = lpest3$beta0, cuad1 = cuest$beta1, lp21 = lpest2$beta1,
lp31 = lpest3$beta1)
}
\keyword{ nonparametric }
\keyword{ smooth }
