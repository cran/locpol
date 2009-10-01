\name{equivKernel}
\alias{equivKernel}
\title{ Equivalent Kernel. }
\description{
  Computes the Equivalent kernel for the local polynomial estimation.
}
\usage{
equivKernel(kernel,nu,deg,lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],
	subdivisions=25)
}
\arguments{
  \item{nu}{ Orders of derivative to estimate. }
  \item{deg}{ Degree of Local polynomial estimator. }
  \item{kernel}{ Kernel used to perform the estimation, see \code{\link{Kernels}}}
  \item{lower, upper}{ Integration limits. }
  \item{subdivisions}{ the maximum number of subintervals. }
}
\details{
	The definition of the Equivalent kernel for the local polynomial 
estimation can be found in page 64 in Fan and Gijbels(1996). The 
implementation uses \code{computeMu} to compute matrix \eqn{S} and then
returns a function object 
}
\value{
	Returns a vector whose components are the equivalent kernel used to
compute the local polynomial estimator for the derivatives in \code{nu}.
}
\references{ 
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman \& Hall, London (1996).
}
\author{ Jorge Luis Ojeda Cabrera. 
}
\seealso{ \code{\link{cteNuK}}, \code{\link{adjNuK}}.
}
\examples{
##	Some kernels and equiv. for higher order
##	compare with p=1 
curve(EpaK(x),-3,3,ylim=c(-.5,1))
f <- equivKernel(EpaK,0,3)
curve(f(x),-3,3,add=TRUE,col="blue")
curve(gaussK(x),-3,3,add=TRUE)
f <- equivKernel(gaussK,0,3)
curve(f(x),-3,3,add=TRUE,col="blue")
##	Draw several Equivalent locl polynomial kernels
curve(EpaK(x),-3,3,ylim=c(-.5,1))
for(p in 1:5){
	curve(equivKernel(gaussK,0,p)(x),-3,3,add=TRUE)
    }
}
\keyword{ nonparametric }
\keyword{ smooth }
