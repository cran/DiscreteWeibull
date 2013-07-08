\name{Edweibull}
\alias{Edweibull}
\alias{E2dweibull}
\alias{Vdweibull}
\alias{ERdweibull}
\title{
Expected values
}
\description{
First and second order moments, variance and expected value of the reciprocal for the type 1 discrete Weibull distribution
}
\usage{
Edweibull(q, beta, eps = 1e-04, nmax = 1000, zero = FALSE)
E2dweibull(q, beta, eps = 1e-04, nmax = 1000, zero = FALSE)
Vdweibull(q, beta, eps = 1e-04, nmax = 1000, zero = FALSE)
ERdweibull(q, beta, eps = 1e-04, nmax = 1000)
}
\arguments{
  \item{q}{
first parameter
}
  \item{beta}{
second parameter
}
  \item{eps}{
error threshold for the numerical computation of the expected value
}
  \item{nmax}{
maximum value considered for the numerical computation of the expected value
}
  \item{zero}{
\code{TRUE}, if the support contains \eqn{0}; \code{FALSE} otherwise
}
}
\details{
The expected values are numerically computed considering a truncated support: integer values smaller than or equal to
\eqn{2\min(nmax;F^{-1}(1-eps;q,\beta))}, where \eqn{F^{-1}} is the inverse of the cumulative distribution function (implemented by the function \code{\link{qdweibull}})
}
\value{
the (approximate) expected values of the discrete Weibull distribution:
\code{Edweibull} gives the first order moment, \code{E2dweibull} the second order moment, \code{Vdweibull} the variance, \code{ERdweibull} the expected value of the reciprocal (only if \code{zero} is \code{FALSE})
}
\author{
Alessandro Barbiero
}
\examples{
q<-0.75
beta<-1.25
Edweibull(q, beta)
E2dweibull(q, beta)
Vdweibull(q, beta)
ERdweibull(q, beta)
# if beta=1...
beta<-1
Edweibull(q, beta)
# which equals...
1/(1-q)
}
\keyword{distribution}