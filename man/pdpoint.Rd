\name{pdpoint}
\alias{pdpoint}
\title{Projection depth of a point with respect to a given bivariate data}
\description{
 \code{pdpoint} function calculates projection depth of a given point with respect to the given  bivariate data.
}
\usage{pdpoint(xpt,x,y)}

\arguments{
    \item{xpt}{An arbitrary bivariate data set for its projection depth with respect to x and y}
    \item{x}{A numeric vector (x coordinates of the given data)}
    \item{y}{A numeric vector (y coordinates of the given data)}
    }
\author{Yijun Zuo, Xiangyang Ye}

\references{
Yijun Zuo, Projection-based Depth Functions and Associated Medians. \emph{The Annals of Statistics} 2003; Vol 31, No. 5, 1460-1490\cr

Yijun zuo, Exactly Compute Bivariate Projection Depth and Stahel-Donoho Estimator. Manuscript, 2009
}



\seealso{
  \code{\link{pdsample}}
  }
\examples{
data(pds)
xpt <- pds[1,]

#one point prediction
pdpoint(xpt,pds$x,pds$y)

#multiple points prediction
pdpoint(pds,pds$x,pds$y)
}

