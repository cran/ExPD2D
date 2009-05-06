\name{pdsample}
\alias{pdsample}
\title{Exact projection depth of points in a given bivariate data }
\description{
 \code{pdsample} function calculates projection depth of a given bivariate data points.
 }

\usage{pdsample(x,y)}

\arguments{
    \item{x}{A numeric vector (the x coordinates of the data)}
    \item{y}{A numeric vector (the y coordinates of the data)}
  }
\author{Yijun Zuo, Xiangyang Ye}
\references{
Yijun Zuo, Projection-based Depth Functions and Associated Medians. \emph{The Annals of Statistics} 2003; Vol 31, No. 5, 1460-1490\cr

Yijun zuo, Exactly Compute Bivariate Projection Depth and Stahel-Donoho Estimator. Manuscript, 2009}

\seealso{
  \code{\link{pdpoint}}
  }
\examples{
data(pds)
pdsample(pds$x,pds$y)
}

