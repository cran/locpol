%% 	OK


\name{pluginBw}
\alias{pluginBw}
\title{ Plugin Bandwidth selector. }
\description{
  Implements a plugin bandwidth selector for the regression function.
}
\usage{
pluginBw(x, y, deg, kernel, weig = rep(1,length(y)))
}
\arguments{
	\item{x}{ x covariate values. }
	\item{y}{ y response values. }
	\item{deg}{ degree of the local polynomial. }
	\item{kernel}{ Kernel used to perform the estimation, see \code{\link{Kernels}}. }
	\item{weig}{ Vector of weights for observations.}
}
\details{
    Computes the plug-in bandwidth selector as shown in Fan and Gijbels(1996) book using pilots estimates as given in page 110-112 (Rule of thumb for bandwidth selection).  Currently, only even values of p are can be used.
}
\value{
	A numeric value.
}
\references{
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman & Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).
}
\author{ Jorge Luis Ojeda Cabrera.
}
\seealso{ \code{\link{thumbBw}}, \code{\link{regCVBwSelC}}.
}
\note{Currently, only even values of p are can be used.
}
\examples{
	size <- 200
	sigma <- 0.25
	deg <- 1
	kernel <- EpaK
	xeval <- 0:100/100
	regFun <- function(x) x^3
	x <- runif(size)
	y <- regFun(x) + rnorm(x, sd = sigma)
	d <- data.frame(x, y)
	cvBwSel <- regCVBwSelC(d$x,d$y, deg, kernel, interval = c(0, 0.25))
	thBwSel <- thumbBw(d$x, d$y, deg, kernel)
	piBwSel <- pluginBw(d$x, d$y, deg, kernel)
	est <- function(bw, dat, x) return(locPolSmootherC(dat$x,dat$y, x, bw, deg,
					kernel)$beta0)
	ise <- function(val, est) return(sum((val - est)^2 * xeval[[2]]))
	plot(d$x, d$y)
	trueVal <- regFun(xeval)
	lines(xeval, trueVal, col = "red")
	xevalRes <- est(cvBwSel, d, xeval)
	cvIse <- ise(trueVal, xevalRes)
	lines(xeval, xevalRes, col = "blue")
	xevalRes <- est(thBwSel, d, xeval)
	thIse <- ise(trueVal, xevalRes)
	xevalRes <- est(piBwSel, d, xeval)
	piIse <- ise(trueVal, xevalRes)
	lines(xeval, xevalRes, col = "blue", lty = "dashed")
	res <- rbind(	bw = c(cvBwSel, thBwSel, piBwSel),
					ise = c(cvIse, thIse, piIse) )
	colnames(res) <- c("CV", "th", "PI")
	res
}
\keyword{ nonparametric }
\keyword{ smooth }
