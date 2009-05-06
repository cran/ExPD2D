\name{pdsample1}
\alias{pdsample1}
\title{Projection depth of a point with respect to a given bivariate data}
\description{
\code{pdsample1} function calculates projection depth of a given bivariate data points. This function is almost the same as \code{pdsample}, but only the object output is the projection depth of the given arbitrary point (not necessary sample point), with respect to the given data. It will be called by \code{pdpoint} function.
}
\usage{pdsample1(x,y)}

\arguments{
    \item{x}{A numeric vector (the x coordinates of the data)}
    \item{y}{A numeric vector (the y coordinates of the data)}
  }
\author{Yijun Zuo, Xiangyang Ye}
\references{
Yijun Zuo, Projection-based Depth Functions and Associated Medians. \emph{The Annals of Statistics} 2003; Vol 31, No. 5, 1460-1490\cr

Yijun zuo, Exactly Compute Bivariate Projection Depth and Stahel-Donoho Estimator. Manuscript, 2009}

\seealso{
  \code{\link{pdsample}}
  }
\examples{
data(pds)
pdsample1(pds$x,pds$y)
}

