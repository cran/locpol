%%	no OK
\name{bivNPest}
\alias{bivDens}
\alias{bivReg}
\alias{plot.bivNpEst}
\alias{predict.bivNpEst}
\alias{epaK2d} 
\alias{gauK2d}
\alias{plotBivNpEstOpts}
\alias{mayBeBwSel}
\title{ Bivariate Local estimation. }
\description{
  Simple bivariate Local density and regression estimation with weights.
}
\usage{
    bivDens(X,weig,K,H)
    bivReg(X,Y,weig,K,H)
    \method{predict}{bivNpEst}(object,newdata,...)
    \method{plot}{bivNpEst}(x,...)
}
\arguments{
    \item{X}{ Covariate or independent data, should be a \code{data.frame} or \code{matrix}, whose two first two columns are used. }
    \item{Y}{ Response data, a \code{vector}}. 	
    \item{weig}{ Vector of weigths for each observations. }	
    \item{K}{ Bivariate kernel function as \code{bivDens} and \code{bivReg}. }	
    \item{H}{ Bandwidth matrix. Its default value is determined by \code{mayBeBwSel}.}
    \item{object, x}{\code{bivNpEst} class objects, those returned by  \code{bivDens} and \code{bivReg} functions.}
    \item{newdata}{Data, should be a \code{data.frame} where the density or regressions is going to be predicted.}
    \item{...}{Further graphical parameters. These parameters should agree with those in \code{\link{persp}}.}
}
\details{
    The functions \code{bivDens} and \code{bivReg} provide a very basic 
interface that allows bivariate local estimation with weights. It implements 
basic kernel density estimator and Nadaraya--Watson estimator for bivariate 
data. Very simple interface methods allow the prediction and plotting of 
these estimators.

    The only bivariate kernels provided are \code{epaK2d} and \code{gauK2d}.
New ones can be added in the same way as functions with a vector of length 2.

    The default bandwidth selector (see \code{mayBeBwSel}) that has been 
provided \emph{is not optimal or good in any sense}. It has been added as 
a simple way to provide an easy, fast and simple way to be able to use the 
estimators. 

    The graphical parameters allowed for \code{...} in \code{plot(x,...)} are those that appears in the function \link{persp}. The list \code{plotBivNpEstOpts} provide a default for some of these graphical parameters.
}
\value{
    A list containing:
    \item{X}{ Covariate data. }
    \item{Y}{ Response data }
    \item{H}{ Bandwidth matrix }
    \item{estFun}{Estimator function. }
}
% \references{ }
\author{ Jorge Luis Ojeda Cabrera. 
}
\examples{
  n <- 100
  d <- data.frame(x=rexp(n,rate=1/2),y=rnorm(n))
  ## x is a length-biased version of an exp. dist. with rate 1.
  dDen <- bivDens(d,weig=1/d$x)
  plot(dDen,r=5)
  d <- data.frame(X1=runif(n),X2=runif(n))
  d$Y <-  exp(10*d$X1+d$X2^2) 
  dDen <- bivDens(d[,c("X1","X2")])
  plot(dDen,r=5)
  dReg <- bivReg(d[,c("X1","X2")],d$Y)
  plot(dReg,r=5)
  plot(dReg,r=5,phi=20,theta=40)
}
\keyword{ nonparametric }
\keyword{ smooth }
