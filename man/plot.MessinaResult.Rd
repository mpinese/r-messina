\name{plot.MessinaResult}
\alias{plot.MessinaResult}
\title{Plot the results of a Messina fit.}
\usage{
  plot.MessinaResult(result, ...)
}
\arguments{
  \item{result}{the result of a Messina fit, either
  returned by messina or messinasurv.}

  \item{i}{an optional index of which feature (as ordered
  by the rows of the x matrix supplied to messina) to plot.
  If missing, the feature with the largest fit margin is
  used.}

  \item{type}{the plotting engine to use, either "base" for
  base graphics, or "ggplot2" for ggplot2 (the ggplot2
  package must be available to use this engine).}
}
\description{
  TODO
}
\author{
  Mark Pinese \email{m.pinese@garvan.org.au}
}
\references{
  TODO
}
\seealso{
  \code{\link{MessinaResult-class}}

  \code{\link{messina}}

  \code{\link{messinaSurv}}

  \code{\link{ggplot2-package}}
}

