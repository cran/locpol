\name{kernelCte}
\alias{cteNuK}
\alias{adjNuK}
\title{ Kernel Constants used in Bandwidth Selection. }
\description{
  These are values depending on the kernel and the local polynomial degrees that are used in bandwidth selection, as proposed in Fan and Gijbels(1996).
}
\usage{
cteNuK(nu,p,kernel,lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],
	subdivisions= 25)
adjNuK(nu,p,kernel,lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],
	subdivisions= 25)
}
\arguments{
  \item{nu}{ Order of derivative to estimate. }
  \item{p}{ Degree of Local polynomial estimator. }
  \item{kernel}{ Kernel used to perform the estimation, see \code{\link{Kernels}}}
  \item{lower, upper}{ Integration limits. }
  \item{subdivisions}{ the maximum number of subintervals. }
}
\details{
	\code{cteNuK} is computed using \link{Compute kernel values} and \code{link{equivKernel}} jointly with the numerical integration utility
\code{\link{integrate}}. \code{adjNuK} is implemented using quotients
of previous functions. See Fan and Gijbels(1996) pages 67 and 119.
}
\value{
	Both functions returns numeric values.
}
\references{ 
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman \& Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).
}
\author{ Jorge Luis Ojeda Cabrera. 
}
\seealso{ \code{\link{regCVBwSelC}}, \code{\link{pluginBw}}, \code{\link{integrate}}.
}
% \examples{}
\keyword{ nonparametric }
\keyword{ smooth }
