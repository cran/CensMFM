\name{rMSN}
\alias{rMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating from Multivariate Skew-normal and Normal Random Distributions.
}
\description{
It generates random realizations from a multivariate Skew-normal and Normal distribution.
}
\usage{
rMSN(n, mu, Sigma, shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{n}{number of observations.}
   \item{mu}{
  a numeric vector of length \eqn{p} representing the location parameter.
  }
  \item{Sigma}{
  a numeric positive definite matrix with dimension \eqn{p}x\eqn{p} representing the scale parameter.
  }
  \item{shape}{
  a numeric vector of length \eqn{p} representing the skewness parameter for Skew-normal(SN) case. If \code{shape == 0}, the SN case reduces to a normal (symmetric) distribution.
  }
}

\value{
It returns a \eqn{n} x \eqn{p} matrix containing the generated random realizations.
}

\references{
Cabral, C. R. B., Lachos, V. H., & Prates, M. O. (2012). Multivariate mixture modeling using skew-normal independent distributions. Computational Statistics & Data Analysis, 56(1), 126-142.

Prates, M. O., Lachos, V. H., & Cabral, C. (2013). mixsmsn: Fitting finite mixture of scale mixture of skew-normal distributions. Journal of Statistical Software, 54(12), 1-20.

C.E. Galarza, L.A. Matos, D.K. Dey & V.H. Lachos. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut.

F.H.C. de Alencar, C.E. Galarza, L.A. Matos & V.H. Lachos. (2019) Finite Mixture Modeling of Censored and Missing Data Using the Multivariate Skew-Normal Distribution. echnical report. ID 19-31. University of Connecticut.
}
\author{
Francisco H. C. de Alencar \email{hildemardealencar@gmail.com},
Christian E. Galarza \email{cgalarza88@gmail.com},
Victor Hugo Lachos \email{hlachos@uconn.edu} and
Larissa A. Matos \email{larissam@ime.unicamp.br}

Maintainer: Francisco H. C. de Alencar \email{hildemardealencar@gmail.com}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fit.FMMSNC}}, \code{\link{rMMSN}} and \code{\link{rMMSN.contour}}
}
\examples{
mu     <- c(-3,-4)
Sigma  <- matrix(c(3,1,1,4.5), 2,2)
shape <- c(-3,2)
rMSN(10,mu = mu,Sigma = Sigma,shape = shape)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Multivariate distributions}% use one of  RShowDoc("KEYWORDS")
\keyword{Skew_normal}% __ONLY ONE__ keyword per line
