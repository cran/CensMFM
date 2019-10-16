\name{rMMSN.contour}
\alias{rMMSN.contour}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise Scatter Plots and Histograms for Finite Mixture of Multivariate Distributions.
}
\description{
It plots the scatter plots with density contours for different multivariate distributions. Possible options are the Skew-normal (\code{family == "SN"}), Normal (\code{family == "Normal"}) and Student-t (\code{family == "t"}) distribution. Different colors are used by groups. Histograms are shown in the diagonal.
}
\usage{
rMMSN.contour(model = NULL, y = NULL, mu = NULL, Sigma = NULL, shape = NULL,
nu = NULL, pii = NULL, Zij = NULL, contour = FALSE, hist.Bin = 30, slice = 100,
col.names = NULL, length.x = c(0.5, 0.5), length.y = c(0.5, 0.5), family = "SN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
 is an object resultant from the code{\link{fit.FMMSNC}} function.
}
  \item{y}{
the response matrix with dimension \eqn{n}x\eqn{p}.
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
  \item{nu}{
the degrees of freedom for the Student-t distribution case, being a vector with dimension \eqn{g}.
}
  \item{pii}{
a vector of weights for the mixture of dimension \eqn{g}, the number of clusters. It must sum to one!
}
  \item{Zij}{
a matrix of dimension \eqn{n}x\eqn{p} indicating the group for each observation.
}
  \item{contour}{
If \code{contour == TRUE} the density contour will be shown, if \code{contour == FALSE} the density contour must be not returned.
}
  \item{hist.Bin}{
number of bins in the histograms. Default is 30.
}
  \item{slice}{
number of slices for the contours.
}
  \item{col.names}{
names passed to the data matrix \eqn{y} of dimension \eqn{p}.
}
  \item{length.x}{
a vector of dimension 2 with the value to be subtracted and added from the minimum and maximum observation in the x-axis respectively. Default is \code{c(0.5,0.5)}.
}
  \item{length.y}{
a vector of dimension 2 with the value to be subtracted and added from the minimum and maximum observation in the y-axis respectively. Default is \code{c(0.5,0.5)}.
}
  \item{family}{
distribution family to be used. Available distributions are the Skew-normal ("SN"), normal ("Normal") or Student-t ("t") distribution.
}
}
\details{
If the \code{model} object is used, the user still has the option to choose the \code{family}. If the \code{model} object is not used, the user must input all other parameters. User may use the \code{rMMSN} function to generate data.
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
This functions works well for any length of \eqn{g} and \eqn{p}, but contour densities are only shown for \eqn{p = 2}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fit.FMMSNC}}, \code{\link{rMMSN}} and \code{\link{fit.FMMSNC}}
}
\examples{
mu          <- Sigma <- shape <- list()
mu[[1]]     <- c(-3,-4)
mu[[2]]     <- c(2,2)
Sigma[[1]]  <- matrix(c(3,1,1,4.5), 2,2)
Sigma[[2]]  <- matrix(c(2,1,1,3.5), 2,2)
shape[[1]]  <- c(-2,2)
shape[[2]]  <- c(-3,4)
nu          <- 0
pii         <- c(0.6,0.4)
percent     <- c(0.1,0.2)
n <- 100
seed <- 654678

set.seed(seed)
test = rMMSN(n = n, pii = pii,mu = mu,Sigma = Sigma,shape = shape,
percent = percent, each = TRUE, family = "SN")


## SN ##
SN.contour = rMMSN.contour(model = NULL, y = test$y, Zij = test$G
,mu = mu, Sigma = Sigma, shape = shape, pii = pii, family = "SN")

#Plotting contours may take some time...
\donttest{
## SN ##
SN.contour = rMMSN.contour(model = NULL, y = test$y, Zij = test$G
,mu = mu, Sigma = Sigma, shape = shape, pii = pii, contour = TRUE,
family = "SN")

## Normal ##
N.contour = rMMSN.contour(model = NULL,y = test$y, Zij = test$G
,mu = mu, Sigma = Sigma, shape = shape, pii = pii, contour = TRUE,
family = "Normal")

## t ##
t.contour = rMMSN.contour(model = NULL,y = test$y, Zij = test$G
,mu = mu, Sigma = Sigma, shape = shape, pii = pii, nu = c(4,3),
contour = TRUE, family = "t")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Multivariate distributions}% use one of  RShowDoc("KEYWORDS")
\keyword{Skew_normal}% __ONLY ONE__ keyword per line
\keyword{Censored data}
\keyword{Scatter plots}
