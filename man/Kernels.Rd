\name{Kernels}
\alias{Kernels}
\alias{biweigK}
\alias{CosK}
\alias{EpaK}
\alias{Epa2K}
\alias{gaussK}
\alias{gaussKlf}
\alias{SqK}
\alias{QuartK}
\alias{TrianK}
\alias{tricubK}
\alias{tricubKlf}
\alias{TriweigK}
\title{ Kernels. }
\description{
  Definition of common kernels used in local polynomial estimation.
}
\usage{
	CosK(x)
	EpaK(x)
	Epa2K(x)
	gaussK(x)
}

\arguments{
  \item{x}{ Numeric vector o value. }
}
\details{
  The implementation of these kernels is done by means
functions that can operate on vectors.

  Most common referred numeric values for these kernels are provided as
attributes, see \code{\link{RK}}, \code{\link{mu0K}}, etc\dots.
}
% \value{  }
\references{
Fan, J. and Gijbels, I.
\emph{ Local polynomial modelling and its applications\/}.
Chapman & Hall, London (1996).

Wand, M.~P. and Jones, M.~C.
\emph{ Kernel smoothing\/}.
Chapman and Hall Ltd., London (1995).
}
\author{ Jorge Luis Ojeda Cabrera. }
\seealso{ \code{\link{RK}}, \code{\link{mu0K}}.
}
\keyword{ nonparametric }
\keyword{ smooth }
