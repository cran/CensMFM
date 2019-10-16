\name{rMMSN}
\alias{rMMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Generator of Finite Mixture of Multivariate Distributions.
}
\description{
It generates random realizations following a multivariate finite mixture of Skew-normal (\code{family == "SN"}) and normal (\code{family == "Normal"}) distributions under censoring. Censoring level can be set as a percentage and it can be adjusted per group if desired.
}
\usage{
rMMSN(n = NULL, mu = NULL, Sigma = NULL, shape = NULL, percent = NULL,
each = FALSE, pii = NULL, family = "SN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of observations
}
  \item{mu}{
a list with \eqn{g} entries, where each entry represents location parameter per group, being a vector of dimension \eqn{p}.
}
  \item{Sigma}{
a list with \eqn{g} entries, where each entry represents a scale parameter per group, a matrix with dimension \eqn{p}x\eqn{p}.
}
  \item{shape}{
a list with \eqn{g} entries, where each entry represents a skewness parameter, being a vector of dimension \eqn{p}.
}
  \item{percent}{
Percentage of censored data in each group or data as a whole (see next item).
}
  \item{each}{
If \code{each == TRUE}, the data will be censored in each group,
where \code{percent} must be a vector of dimension \eqn{p}. Besides, if \code{each == FALSE} (by default), the data will be censored in the whole set, then \code{percent} must be a vector of dimension 1.
}
  \item{pii}{
a vector of weights for the mixture of dimension \eqn{g}, the number of clusters. It must sum to one!
}
  \item{family}{
distribution family to be used for fitting. Options are "SN" for the Skew-normal and "Normal" for the Normal distribution respectively.
}
}

\value{
It returns a list that depending of the case, it returns one or more of the following objects:
  \item{y}{a \eqn{n} x \eqn{p} matrix containing the generated random realizations.}
  \item{G}{a vector of length \eqn{n} containing the group classification per subject.}
  \item{cutoff}{a vector containing the censoring cutoffs per group.}
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
\code{\link{fit.FMMSNC}}, \code{\link{rMSN}} and \code{\link{rMMSN.contour}}
}
\examples{
mu <- Sigma <- shape <- list()
mu[[1]]     <- c(-3,-4)
mu[[2]]     <- c(2,2)
shape[[1]]  <- c(-2,2)
shape[[2]]  <- c(-3,4)
Sigma[[1]]  <- matrix(c(3,1,1,4.5), 2,2)
Sigma[[2]]  <- matrix(c(2,1,1,3.5), 2,2)
pii         <- c(0.6,0.4)
percent   <- c(0.1,0.1)
family <- "SN"
n <-100

set.seed(20)
rMMSN(n = n,pii = pii, mu = mu, Sigma = Sigma, shape = shape,
percent = percent, each = TRUE, family = family)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Multivariate distributions}% use one of  RShowDoc("KEYWORDS")
\keyword{Skew_normal}% __ONLY ONE__ keyword per line
\keyword{Censored data}
