\name{KernelChars}
\alias{K4}
\alias{RK}
\alias{RdK}
\alias{mu0K}
\alias{mu2K}
\alias{dom}
\alias{Kernel characteristics}
\title{Kernel characteristics}
\description{
	For a given kernel returns some of the most commonly used numeric values 
for functionals related to them.
}
\usage{
	RK(K)
	RdK(K)
	mu2K(K)
	mu0K(K)
	K4(K)
	dom(K)
}
\arguments{
  \item{K}{ A kernel as given in \code{\link{Kernels}} }
}
\details{
	These functions returns some of the most commonly used numeric 
values related to the kernels. They are implemented as an attribute 
of every kernel. For the computations of the numeric value for these 
quantities, see references.
}
\value{
	A numeric value returning:
\itemize{
\item{RK}{The $L_2$ norm of \code{K}.}
\item{RdK}{The $L_2$ norm of the derivative of \code{K}.}
\item{mu2K}{The second order moment of \code{K}.}
\item{mu2K}{The zeroth order moment of \code{K}.}
\item{dom}{The support of \code{K}.}
\item{K4}{The fourth order autoconvolution of \code{K} at \eqn{x=0}{x=0}.}
}
These functions are implemented via kernel attributes.
}
\references{ 
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman \& Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).
}
\author{ Jorge Luis Ojeda Cabrera.}
\seealso{ 
\code{\link{Kernels}}, \link{Compute kernel values}. 
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
