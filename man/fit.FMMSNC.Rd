\name{fit.FMMSNC}
\alias{fit.FMMSNC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Finite Mixture of Multivariate Distributions.
}
\description{
It adjusts a finite mixture of censored and/or missing multivariate distributions (FM-MC). These are the Skew-normal, normal and Student-t multivariate distributions. It uses a EM-type algorithm for iteratively computing maximum likelihood estimates of the parameters.
}
\usage{
fit.FMMSNC(cc, LI, LS, y, mu = NULL, Sigma = NULL, shape = NULL, pii = NULL,
nu = NULL, g = NULL, get.init = TRUE, criteria = TRUE, family = "SN", error = 1e-05,
iter.max = 350, uni.Gama = FALSE, kmeans.param = NULL, cal.im = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cc}{
vector of censoring indicators. For each observation it takes 0 if non-censored, 1 if censored.
}
  \item{LI}{
the matrix of lower limits of dimension \eqn{n}x\eqn{p}. See details section.
}
  \item{LS}{
the matrix of upper limits of dimension \eqn{n}x\eqn{p}. See details section.
}
  \item{y}{
the response matrix with dimension \eqn{n}x\eqn{p}.
}
  \item{mu}{
a list with \eqn{g} entries, where each entry represents location parameter per group, being a vector of dimension. \eqn{p}.
}
  \item{Sigma}{
a list with \eqn{g} entries, where each entry represents a scale parameter per group, a matrix with dimension. \eqn{p}x\eqn{p}.
}
  \item{shape}{
a list with \eqn{g} entries, where each entry represents a skewness parameter, being a vector of dimension \eqn{p}.
}
  \item{pii}{
a vector of weights for the mixture (dimension of the number \eqn{g} of clusters). Must sum to one!
}
  \item{nu}{
the degrees of freedom for the Student-t distribution case, being a vector with dimension \eqn{g}.
}
  \item{g}{
number of mixture components.
}
  \item{get.init}{
Logical, \code{TRUE} or \code{FALSE}. If (\code{get.init==TRUE}) the function computes the initial values, otherwise (\code{get.init==FALSE}) the user should enter the initial values manually.
}
  \item{criteria}{
Logical, \code{TRUE} or \code{FALSE}. It indicates if likelihood-based criteria selection methods (AIC, BIC and EDC) are computed for comparison purposes.
}
  \item{family}{
distribution family to be used. Available distributions are the Skew-normal ("SN"), normal ("Normal") or Student-t ("t") distribution.
}
  \item{error}{
relative error for stopping criterion of the algorithm. See details.
}
  \item{iter.max}{
the maximum number of iterations of the EM algorithm.
}
  \item{uni.Gama}{
Logical, \code{TRUE} or \code{FALSE}. If \code{uni.Gama==TRUE}, the scale matrices per group are considered to be equals.
}
  \item{kmeans.param}{
a list with alternative parameters for the kmeans function when generating initial values. List by default is
\code{list(iter.max = 10, n.start = 1, algorithm = "Hartigan-Wong")}.
}
  \item{cal.im}{
Logical, \code{TRUE} or \code{FALSE}. If \code{cal.im==TRUE}, the information matrix is calculated and the standard errors are reported.
}
}
\details{
The information matrix is calculated with respect to the entries of
the square root matrix of Sigma, this using the Empirical information matrix. Disclaimer: User must be careful since the inference is asymptotic, so it must be used for decent sample sizes. Stopping criterion is \code{abs((loglik/loglik-1))<epsilon}.
}
\value{
It returns a list that depending of the case, it returns one or more of the following objects:
  \item{mu}{a list with \eqn{g} components, where each component is a vector with dimension \eqn{p} containing the estimated values of the location parameter.}
  \item{Sigma}{a list with \eqn{g} components, where each component is a matrix with dimension \eqn{p}x\eqn{p} containing the estimated values of the scale matrix.}
  \item{Gamma}{a list with \eqn{g} components, where each component is a matrix with dimension \eqn{p}x\eqn{p} containing the estimated values of the \eqn{Gamma} scale matrix.}
  \item{shape}{a list with \eqn{g} components, where each component is a vector with dimension \eqn{p} containing the estimated values of the skewness parameter.}
  \item{nu}{a vector with one element containing the value of the degreees of freedom \eqn{nu} parameter.}
  \item{pii}{a vector with \eqn{g} elements containing the estimated values of the weights \eqn{pii}.}
   \item{Zij}{a \eqn{n} x \eqn{p} matrix containing the estimated weights values of the subjects for each group.}
  \item{yest}{a \eqn{n} x \eqn{p} matrix containing the estimated values of \eqn{y}.}
  \item{MI}{a list with the standard errors for all parameters.}
  \item{logLik}{the log-likelihood value for the estimated parameters.}
  \item{aic}{the AIC criterion value for the estimated parameters.}
  \item{bic}{the BIC criterion value for the estimated parameters.}
  \item{edc}{the EDC criterion value for the estimated parameters.}
  \item{iter}{number of iterations until the EM algorithm converges.}
  \item{group}{a \eqn{n} x \eqn{p} matrix containing the classification for the subjects to each group.}
  \item{time}{time in minutes until the EM algorithm converges.}
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
\note{
The \code{uni.Gama} parameter refers to the \eqn{\Gamma} matrix for the Skew-normal distribution, while for the normal and student-t distribution, this parameter refers to the \eqn{\Sigma} matrix.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rMSN}}, \code{\link{rMMSN}} and \code{\link{rMMSN.contour}}
}
\examples{
mu          <- Sigma <- shape <- list()
mu[[1]]     <- c(-3,-4)
mu[[2]]     <- c(2,2)
Sigma[[1]]  <- matrix(c(3,1,1,4.5), 2,2)
Sigma[[2]]  <- matrix(c(2,1,1,3.5), 2,2)
shape[[1]]  <- c(-2,2)
shape[[2]]  <- c(-3,4)
nu          <- c(0,0)
pii         <- c(0.6,0.4)
percen <- c(0.1,0.2)
n <- 200
g <- 2
seed <- 654678

set.seed(seed)
test = rMMSN(n = n, pii = pii,mu = mu,Sigma = Sigma,shape = shape,
percen = percen, each = TRUE, family = "SN")

Zij <- test$G
cc <- test$cc
y <- test$y

## left censure ##
LI <-cc
LS <-cc
LI[cc==1]<- -Inf
LS[cc==1]<- y[cc==1]

test_fit.cc0 = fit.FMMSNC(cc, LI, LS, y, mu=mu,
Sigma = Sigma, shape=shape, pii = pii, g = 2, get.init = FALSE,
criteria = TRUE, family = "Normal", error = 0.0001,
iter.max = 200, uni.Gama = FALSE, cal.im = FALSE)

#full analysis may take a few seconds more...

\donttest{
test_fit.cc = fit.FMMSNC(cc, LI, LS, y, mu=mu,
Sigma = Sigma, shape=shape, pii = pii, g = 2, get.init = FALSE,
criteria = TRUE, family = "SN", error = 0.00001,
iter.max = 350, uni.Gama = FALSE, cal.im = TRUE)

## missing data ##
pctmiss <- 0.2 # 20\% of missing data in the whole data
missing <- matrix(runif(n*g), nrow = n) < pctmiss
y[missing] <- NA

cc <- matrix(nrow = n,ncol = g)
cc[missing] <- 1
cc[!missing] <- 0

LI <- cc
LS <-cc
LI[cc==1]<- -Inf
LS[cc==1]<- +Inf

test_fit.mis = fit.FMMSNC(cc, LI, LS, y, mu=mu,
Sigma = Sigma, shape=shape, pii = pii, g = 2, get.init = FALSE,
criteria = TRUE, family = "SN", error = 0.00001,
iter.max = 350, uni.Gama = FALSE, cal.im = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Multivariate distributions}% use one of  RShowDoc("KEYWORDS")
\keyword{Skew_normal}% __ONLY ONE__ keyword per line
\keyword{Censored data}
\keyword{Missing data}
\keyword{Finite Mixture}
