%%	OK

\name{denCVBwSelC}
\alias{denCVBwSelC}
\title{ CV bandwidth selector for density }
\description{
	Computes Cross Validation bandwidth selector for the 
Parzen--Rosenblatt density estimator...
}
\usage{
denCVBwSelC(x, kernel = gaussK, weig = rep(1, length(x)), 
            interval = .lokestOptInt)
}
\arguments{
	\item{x}{ vector with data points. }
	\item{kernel}{ Kernel used to perform the estimation, see \code{\link{Kernels}}. }
	\item{weig}{ Vector of weights for observations. }
	\item{interval}{ A range of values where to look for the bandwidth 
parameter. }
}
\details{
	The selector is implemented using its definition.
}
\value{
	A numeric value with the bandwidth.
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
\seealso{ \code{\link{bw.nrd0}}, \code{\link{dpik}}.
}
\examples{
stdy <- function(size=100,rVar=rnorm,dVar=dnorm,kernel=gaussK,x=NULL)
{
	if( is.null(x) ) x <- rVar(size)
	Tc <- system.time( dbwc <- denCVBwSelC(x,kernel) )[3]
	ucvT <- system.time( ucvBw <- bw.ucv(x,lower=0.00001,upper=2.0) )[3]
	nrdT <- system.time( nrdBw <- bw.nrd(x) )[3]
	{
		xeval <- seq( min(x)+dbwc , max(x)-dbwc ,length=50)
		hist(x,probability=TRUE)
		lines(xeval,trueDen <- dVar(xeval),col="red")
		lines(density(x),col="cyan")
		xevalDenc <- PRDenEstC(x,xeval,dbwc,kernel)
		dencMSE <- mean( (trueDen-xevalDenc)^2 )
		xevalucvDen <- PRDenEstC(x,xeval,ucvBw,kernel)
		ucvMSE <- mean( (trueDen-xevalucvDen)^2 )
		xevalDenNrd <- PRDenEstC(x,xeval,nrdBw,kernel)
		nrdMSE <- mean( (trueDen-xevalDenNrd)^2 )
		lines(xevalDenc,col="green")
		lines(xevalucvDen,col="blue")
		lines(xevalDenNrd,col="grey")
	}
	return( cbind(  bwVal=c(evalC=dbwc,ucvBw=ucvBw,nrdBw=nrdBw),
				mse=c(dencMSE,ucvMSE,nrdMSE), 
				time=c(Tc,ucvT,nrdT) ) )
}
stdy(100,kernel=gaussK) 
stdy(100,rVar=rexp,dVar=dexp,kernel=gaussK)
##	check stdy with other kernel, distributions	
}
\keyword{ nonparametric }
\keyword{ smooth }

