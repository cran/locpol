% 	OK

\name{selKernel}
\alias{selKernel}
\title{ Kernel selection.}
\description{
  Uses kernel attributes to selects kernels. This function is mainly used for internal purposes.
}
\usage{
selKernel(kernel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernel}{ kernel to use. }
}
\details{
	Uses \code{RK(K)} to identify a kernel. 
	The integer is used in the C code part to perform 
computations with given kernel. It allows for a kernel 
selection in C routines. It is used only for internal 
purposes. 
}
\value{
  An integer that is unique for each kernel.
}
% \references{ }
\author{ Jorge Luis Ojeda Cabrera. 
}
% \note{ }
\section{Warning }{  
	Used only for internal purposes. 
}
% \seealso{}

\keyword{ nonparametric }
\keyword{ smooth }
