%%	OK

\name{thumbBw}
\alias{thumbBw}
\alias{compDerEst}
\title{ Rule of thumb for bandwidth selection.}
\description{
  Implements Fan and Gijbels(1996)'s Rule of thumb for bandwidth selection
}
\usage{
thumbBw(x, y, deg, kernel, weig = rep(1, length(y)))
compDerEst(x, y, p, weig = rep(1, length(y)))
}
\arguments{
  \item{x}{ x covariate data values. }
  \item{y}{ y response data values. }
  \item{p}{ order of local polynomial estimator. }
  \item{deg}{ Local polynomial estimation degree($p$). }
  \item{kernel}{ Kernel used to perform the estimation.}
  \item{weig}{ weights if they are required. }
}
\details{
	See Fan and Gijbels(1996) book, Section 4.2. This implementation is 
also considering weights. \code{compDerEst} computes the \eqn{p+1} 
derivative of the regression function in a simple manner, assuming it
is a polynomial in \eqn{x}. \code{thumbBw} gives a bandwidth selector
by means of pilot estimator given by \code{compDerEst} and the mean of
residuals. 
}
\value{
	 \code{thumbBw} returns a single numeric value, while \code{compDerEst} returns a data frame whose components are:
	\item{x}{x values.}
	\item{y}{y values.}
	\item{res}{residuals for the parametric estimation.}
	\item{der}{derivative estimation at x values.}
}
\references{ 
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman \& Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).
}
\author{ 
	Jorge Luis Ojeda Cabrera. 
}
\seealso{ 
	\code{\link{regCVBwSelC}}, \code{\link{pluginBw}}.
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
