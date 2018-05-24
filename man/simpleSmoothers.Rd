%%	OK

\name{simpleSmoothers}
\alias{simpleSmootherC}
\alias{simpleSqSmootherC}
\title{ Simple smoother }
\description{
  Computes simple kernel smoothing  
}
\usage{
  simpleSmootherC(x, y, xeval, bw, kernel, weig = rep(1, length(y)))
  simpleSqSmootherC(x, y, xeval, bw, kernel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ x covariate data values. }
  \item{y}{ y response data values. }
  \item{xeval}{ Vector with evaluation points. }
  \item{bw}{ Smoothing parameter, bandwidth. }
  \item{kernel}{ Kernel used to perform the estimation, see \code{\link{Kernels}}}
  \item{weig}{ weights if they are required. }
}
\details{
  Computes simple smoothing, that is to say: it averages \code{y} 
values times kernel evaluated on \code{x} values. \code{simpleSqSmootherC}
does the average with the square of such values.
}
\value{
  Both functions returns a \code{data.frame} with
  \item{ x }{ \eqn{x} evaluation points. }
  \item{ reg }{ the smoothed values at \code{x} points. }
  ...
}
% \references{ }
\author{
  Jorge Luis Ojeda Cabrera. 
}
% \note{ }

\seealso{ 
  \code{\link{PRDenEstC}}, \code{\link{Kernel characteristics}}
}
\examples{
	size <- 1000
	x <- runif(100)
	bw <- 0.125
	kernel <- EpaK
	xeval <- 1:9/10
	y <- rep(1,100)	
	##	x kern. aver. should give density f(x)
	prDen <- PRDenEstC(x,xeval,bw,kernel)$den
	ssDen <- simpleSmootherC(x,y,xeval,bw,kernel)$reg
	all(abs(prDen-ssDen)<1e-15)
	##	x kern. aver. should be f(x)*R2(K) aprox.
	s2Den <- simpleSqSmootherC(x,y,xeval,bw,kernel)$reg
	summary( abs(prDen*RK(kernel)-s2Den) )
	summary( abs(1*RK(kernel)-s2Den) )
	##	x kern. aver. should be f(x)*R2(K) aprox.
	for(n in c(1000,1e4,1e5))
	{
		s2D <- simpleSqSmootherC(runif(n),rep(1,n),xeval,bw,kernel)$reg
		cat("\n",n,"\n")
		print( summary( abs(1*RK(kernel)-s2D) ) )
	}
}
\keyword{ nonparametric }
\keyword{ smooth }
