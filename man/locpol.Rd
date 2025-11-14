%%	OK
\name{locpol}
\alias{locpol}
\alias{confInterval}
\alias{residuals.locpol}
\alias{fitted.locpol}
\alias{summary.locpol}
\alias{print.locpol}
\alias{plot.locpol}
\title{ Local Polynomial estimation. }
\description{
  Formula interface for the local polynomial estimation.
}
\usage{
    locpol(formula,data,weig=rep(1,nrow(data)),bw=NULL,kernel=EpaK,deg=1,
            xeval=NULL,xevalLen=100)
    confInterval(x)
    \method{residuals}{locpol}(object,...)
    \method{fitted}{locpol}(object,deg=0,...)
    \method{summary}{locpol}(object,...)
    \method{print}{locpol}(x,...)
    \method{plot}{locpol}(x,...)
}
\arguments{
	\item{formula}{ formula as in \code{lm}, only first covariate is used. }
	\item{data}{ data frame with data. }
	\item{weig}{ Vector of weights for each observations. }
	\item{bw}{ Smoothing parameter, bandwidth. }
	\item{kernel}{ Kernel used to perform the estimation, see
        \code{\link{Kernels}}}
	\item{deg}{ Local polynomial estimation degree (\eqn{p}). }
	\item{xeval}{ Vector of evaluation points. By default \code{xevalLen} points between the min. and the max. of the regressors.}
	\item{xevalLen}{ Length of \code{xeval} if it is \code{NULL} }
    \item{x}{ A \code{locpol} object. }
    \item{object}{ A \code{locpol} object. }
    \item{...}{ Any other required argument. }
}
\details{
	This is an interface to the local polynomial estimation function
    that provides basic \code{lm} functionality. \code{summary} and
    \code{print} methods shows very basic information about the fit,
    \code{fitted} return the estimation of the derivatives if \code{deg}
    is larger than 0, and \code{plot} provides a plot of data, local
    polynomial estimation and the variance estimation.

    Variance estimation is carried out by means of the local constant
    regression estimation of the squared residuals.

    \code{confInterval} provides confidence intervals for all points
    in \code{x$lpFit$[,x$X]}, say those in \code{xeval}.
}
\value{
	A list containing among other components:
    \item{mf}{ Model frame for \code{data} and \code{formula}. }
    \item{data}{ data frame with data. }
    \item{weig}{ Vector of weight for each observations.}
    \item{xeval}{ Vector of evaluation points. }
    \item{bw}{ Smoothing parameter, bandwidth. }
    \item{kernel}{ Kernel used, see \code{\link{Kernels}}}
    \item{KName}{ Kernel name, a string with the name of kernel.}
    \item{deg}{ Local polynomial estimation degree (\eqn{p}).}
    \item{X,Y}{ Names in \code{data} of the response and covariate. They are also used in \code{lpFit} to name the fitted data.}
    \item{residuals}{ Residuals of the local polynomial fit.}
    \item{lpFit}{ Data frame with the local polynomial fit. It contains covariate, response, derivatives estimation, \eqn{X} density estimation, and variance estimation.}
}
\references{
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman & Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).

Crist'obal, J. A. and Alcal\'a, J. T. (2000).
\emph{ Nonparametric regression estimators for length biased data\/}.
J. Statist. Plann. Inference, 89, pp. 145-168.

Ahmad, Ibrahim A. (1995)
\emph{On multivariate kernel estimation for samples from weighted
distributions\/}.
Statistics & Probability Letters, 22, num. 2, pp. 121-129
}
\author{ Jorge Luis Ojeda Cabrera.
}
\seealso{ \code{\link[KernSmooth]{locpoly}} from package \pkg{KernSmooth},
\code{\link{ksmooth}} and \code{\link{loess}} in \pkg{stats} (but from earlier package \code{modreg}).
}
\examples{
    N <- 250
    xeval <- 0:100/100
    ##  ex1
    d <- data.frame(x = runif(N))
    d$y <- d$x^2 - d$x + 1 + rnorm(N, sd = 0.1)
    r <- locpol(y~x,d)
    plot(r)
    ##  ex2
    d <- data.frame(x = runif(N))
    d$y <- d$x^2 - d$x + 1 + (1+d$x)*rnorm(N, sd = 0.1)
    r <- locpol(y~x,d)
    plot(r)
    ## notice:
    rr <- locpol(y~x,d,xeval=runif(50,-1,1))
    ## notice x has null dens. outside (0,1)
    ## plot(rr) raises an error, no conf. bands outside (0,1).
    ## length biased data !!
    d <- data.frame(x = runif(10*N))
    d$y <- d$x^2 - d$x + 1 + (rexp(10*N,rate=4)-.25)
    posy <- d$y[ whichYPos <- which(d$y>0) ];
    d <- d[sample(whichYPos, N,prob=posy,replace=FALSE),]
    rBiased <- locpol(y~x,d)
    r <- locpol(y~x,d,weig=1/d$y)
    plot(d)
    points(r$lpFit[,r$X],r$lpFit[,r$Y],type="l",col="blue")
    points(rBiased$lpFit[,rBiased$X],rBiased$lpFit[,rBiased$Y],type="l")
    curve(x^2 - x + 1,add=TRUE,col="red")

}
\keyword{ nonparametric }
\keyword{ smooth }
