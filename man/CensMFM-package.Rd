\name{CensMFM-package}
\alias{CensMFM-package}
\alias{CensMFM}
\docType{package}
\title{
\packageTitle{CensMFM}
}
\description{
\packageDescription{CensMFM}
}
\details{

The DESCRIPTION file:

\packageIndices{CensMFM}

  The \pkg{CensMFM} package provides comprehensive tools for fitting and analyzing finite mixture models on censored and/or missing data using several multivariate distributions. This package supports the normal, Student-t, and skew-normal distributions, facilitating point estimation and asymptotic inference through the empirical information matrix. Additionally, it allows for the generation of censored data.

  Key functions include:
  \itemize{
    \item \code{\link{fit.FMMSNC}}: Fits finite mixtures of censored and/or missing multivariate distributions using an EM-type algorithm. This function supports skew-normal, normal, and Student-t distributions.
    \item \code{\link{rMMSN.contour}}: Generates pairwise scatter plots and contour plots for analyzing the relationships within the fitted models.
    \item \code{\link{rMMSN}}: Provides functionality to generate random realizations from a finite mixture of multivariate distributions, particularly useful for simulation studies involving censored data.
    \item \code{\link{rMSN}}: Focuses on generating random realizations from multivariate Skew-normal and Normal distributions.
  }

  This package serves as an extension and complement to the methodologies presented in the paper by Lachos, V. H., Moreno, E. J. L., Chen, K. & Cabral, C. R. B. (2017) <doi:10.1016/j.jmva.2017.05.005>, specifically for the multivariate skew-normal case.


}
\author{
\packageAuthor{CensMFM}

Maintainer: \packageMaintainer{CensMFM}
}
\references{
Cabral, C. R. B., Lachos, V. H., & Prates, M. O. (2012). Multivariate mixture modeling using skew-normal independent distributions. Computational Statistics & Data Analysis, 56(1), 126-142.

Prates, M. O., Lachos, V. H., & Cabral, C. (2013). mixsmsn: Fitting finite mixture of scale mixture of skew-normal distributions. Journal of Statistical Software, 54(12), 1-20.

C.E. Galarza, L.A. Matos, D.K. Dey & V.H. Lachos. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut.

F.H.C. de Alencar, C.E. Galarza, L.A. Matos & V.H. Lachos. (2019) Finite Mixture Modeling of Censored and Missing Data Using the Multivariate Skew-Normal Distribution. echnical report. ID 19-31. University of Connecticut.
}

\keyword{Multivariate distributions}% use one of  RShowDoc("KEYWORDS")
\keyword{Skew_normal}% __ONLY ONE__ keyword per line
\keyword{Censored data}
\keyword{Missing data}
\keyword{Finite Mixture}
\seealso{
\code{\link{fit.FMMSNC}}, \code{\link{rMSN}}, \code{\link{rMMSN}} and \code{\link{rMMSN.contour}}
}
