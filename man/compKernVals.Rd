\name{compKernVals}
\alias{computeRK}
\alias{computeK4}
\alias{computeMu}
\alias{computeMu0}
\alias{computeMu}
\alias{Kconvol}
\alias{Compute kernel values}
\title{Compute kernel values.}
\description{
  Some \R code provided to compute kernel related values.
}
\usage{
	computeRK(kernel, lower=dom(kernel)[[1]], upper=dom(kernel)[[2]], 
			subdivisions = 25)
	computeK4(kernel, lower=dom(kernel)[[1]], upper=dom(kernel)[[2]], 
			subdivisions = 25)
	computeMu(i, kernel, lower=dom(kernel)[[1]], upper=dom(kernel)[[2]], 
			subdivisions = 25)
	computeMu0(kernel, lower=dom(kernel)[[1]], upper=dom(kernel)[[2]], 
			subdivisions = 25)
	Kconvol(kernel,lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],
			subdivisions = 25)
}
\arguments{
  \item{kernel}{ Kernel used to perform the estimation, see \code{\link{Kernels}}}
  \item{i}{ Order of kernel moment to compute }
  \item{lower, upper}{ Integration limits. }
  \item{subdivisions}{ the maximum number of subintervals. }
}
\details{
  These functions uses function \code{\link{integrate}}.
}
\value{
	A numeric value returning:
\item{computeK4}{The fourth order autoconvolution of \code{K}.}
\item{computeRK}{The second order autoconvolution of \code{K}.}
\item{computeMu0}{The integral of \code{K}.}
\item{computeMu2}{The second order moment of \code{K}.}
\item{computeMu}{The \eqn{i}-th order moment of \code{K}.}
\item{Kconvol}{The autoconvolution of \code{K}.}
These functions are implemented by means of \code{\link{integrate}}.
}
\references{ 
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman \& Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).
}
\author{ Jorge Luis Ojeda Cabrera. }
\seealso{ \code{\link{RK}}, \link{Kernel characteristics}, \code{\link{integrate}}.
}
\examples{
	##	Note that lower and upper params are set in the definition to
	##	use 'dom()' function.
	g <- function(kernels)
	{
		mu0 <- sapply(kernels,function(x) computeMu0(x,))
		mu0.ok <- sapply(kernels,mu0K)
		mu2 <- sapply(kernels,function(x) computeMu(2,x))
		mu2.ok <- sapply(kernels,mu2K)
		Rk.ok <- sapply(kernels,RK)
		RK <- sapply(kernels,function(x) computeRK(x))
		K4 <- sapply(kernels,function(x) computeK4(x))
		res <- data.frame(mu0,mu0.ok,mu2,mu2.ok,RK,Rk.ok,K4)
		res
	}
	g(kernels=c(EpaK,gaussK,TriweigK,TrianK))
}
\keyword{ nonparametric }
\keyword{ smooth }
